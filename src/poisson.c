#include "poisson.h"

static inline REAL sqr(REAL x)
{
    return x*x;
}

void applyPboundaryCond(REAL **P, lattice *grid, char **FLAG)
{
    /* Apply boundary conditions for the pressure field */
    short flag = 0;
    /* First set values on the actual boundary of the region */
    if (grid->edges & LEFT)
        for (int j = 0; j <= grid->delj; j++)
            P[0][j] = P[1][j];
    if (grid->edges & RIGHT)
        for (int j = 0; j <= grid->delj; j++)
            P[grid->deli+1][j] = P[grid->deli][j];
    if (grid->edges & BOTTOM)
        for (int i = 0; i <= grid->deli; i++)
            P[i][0] = P[i][1];
    if (grid->edges & TOP)
        for (int i = 0; i <= grid->deli; i++)
            P[i][grid->delj+1] = P[i][grid->delj];
    REAL dxSqrd = sqr(grid->delx);
    REAL dySqrd = sqr(grid->dely);
    for (int i = 1; i <= grid->deli; i++)
        for (int j = 1; j <= grid->delj; j++)
        {
            flag = FLAG[i][j];
            if (flag == C_F)
                continue;
            else if (flag == C_B)
                continue;
            switch (flag ^ C_B)
            {
            case B_N:
                P[i][j] = P[i][j+1];
                break;
            case B_O:
                P[i][j] = P[i+1][j];
                break;
            case B_S:
                P[i][j] = P[i][j-1];
                break;
            case B_W:
                P[i][j] = P[i-1][j];
                break;
            case (B_N | B_O):
                P[i][j] = (dxSqrd*P[i][j+1] + dySqrd*P[i+1][j])/(dxSqrd + dySqrd);
                break;
            case (B_N | B_W):
                P[i][j] = (dxSqrd*P[i][j+1] + dySqrd*P[i-1][j])/(dxSqrd + dySqrd);
                break;
            case (B_S | B_O):
                P[i][j] = (dxSqrd*P[i][j-1] + dySqrd*P[i+1][j])/(dxSqrd + dySqrd);
                break;
            case (B_S | B_W):
                P[i][j] = (dxSqrd*P[i][j-1] + dySqrd*P[i-1][j])/(dxSqrd + dySqrd);
                break;
            default:
                P[i][j] = 0;
            }

        }
    return;
}

int solveSORforPoisson(REAL **p, REAL **rhs, char **FLAG,
                       fluidSim *sim, lattice *grid, MPI_Comm Region)
{
    /* Use a SOR algorithm to solve the poisson equation */
    int  i, j, it = 0;
    REAL temporary, error, eps = sqr(sim->eps);
    REAL invXWidthSqrd = 1/sqr(grid->delx);
    REAL invYWidthSqrd = 1/sqr(grid->dely);
    REAL scale = sim->omega/(2*(invXWidthSqrd+invYWidthSqrd));
    REAL buf[grid->deli + grid->delj + 2];
    /* Count the number of fluid cells */
    int numberOfCells = 0;
    for (i = 1; i <= grid->deli; i++)
        for (j = 1; j <= grid->delj; j++)
        {
            if (FLAG[i][j] == C_F)
                numberOfCells++;
        }
    temporary = eps*numberOfCells;
    MPI_Allreduce(&temporary,&eps,1,MPI_DOUBLE,MPI_SUM,Region);
    do {
        /* Apply the boundary condition */
        applyPboundaryCond(p,grid,FLAG);
        exchangeMat(p,1,1,buf,grid,Region);
        /* Compute the new coefficients iteratively */
        for (i = 1; i <= grid->deli; i++)
            for (j = 1; j <= grid->delj; j++)
            {
                if (FLAG[i][j] != C_F)
                    continue;
                temporary = invXWidthSqrd*(p[i+1][j] + p[i-1][j])
                          + invYWidthSqrd*(p[i][j+1] + p[i][j-1])
                          - rhs[i-1][j-1];
                p[i][j] = scale*temporary + (1-sim->omega)*p[i][j];
            }
        error = 0;
        /* Calculate the residue with respect to the rhs field */
        for (i = 1; i <= grid->deli; i++)
            for (j = 1; j <= grid->delj; j++)
            {
                if (FLAG[i][j] != C_F)
                    continue;
                temporary = invXWidthSqrd*(p[i+1][j] - 2*p[i][j] + p[i-1][j])
                          + invYWidthSqrd*(p[i][j+1] - 2*p[i][j] + p[i][j-1])
                          - rhs[i-1][j-1];
                error += sqr(temporary);
            }
        if (++it == sim->itmax)
            break;
        temporary = error;
        MPI_Allreduce(&temporary,&error,1,MPI_DOUBLE,MPI_SUM,Region);
    } while (error > eps);
    return it;
}

REAL compDelt(lattice *grid, REAL **U, REAL **V, fluidSim *sim)
{
    /* Find the optimal step width in time */
    if (sim->tau <= 0)
        return sim->dt;
    REAL dt = sim->Re*(sqr(grid->delx) + sqr(grid->dely))/2;
    REAL utime, vtime;
    for (int i = 1; i <= grid->deli; i++)
        for (int j = 1; j <= grid->delj; j++)
        {
            utime = grid->delx/U[i+1][j];
            if (utime < 0)
                utime = -utime;
            /* One of these inequalities will always be trivial */
            if (utime < dt)
                dt = utime;
            vtime = (grid->dely)/V[i][j+1];
            if (vtime < 0)
                vtime = -vtime;
            if (vtime < dt)
                dt = vtime;
        }
    dt = sim->tau*dt;
    if (dt > sim->dt)
        return sim->dt;
    return dt;
}

void compRHS(REAL **F, REAL **G, REAL **RHS, char **FLAG, lattice *grid, REAL delt)
{
    for (int i = 0; i < grid->deli; i++)
        for (int j = 0; j < grid->delj; j++)
        {
            if (FLAG[i][j] != C_F)
                continue;
            RHS[i][j] = (F[i+1][j+1] - F[i][j+1])/(delt*(grid->delx));
            RHS[i][j] += (G[i+1][j+1] - G[i+1][j])/(delt*(grid->dely));
        }
    return;
}

void adaptUV(REAL **U, REAL **V, REAL **P, REAL **F, REAL **G,
             REAL delt, char **FLAG, lattice *grid)
{
    REAL facX = delt/(grid->delx);
    REAL facY = delt/(grid->dely);
    for (int i = 1; i <= grid->deli; i++)
        for (int j = 1; j <= grid->delj; j++)
        {
            if (FLAG[i][j] != C_F)
                continue;
            U[i+1][j] = F[i][j] - facX*(P[i+1][j] - P[i][j]);
            V[i][j+1] = G[i][j] - facY*(P[i][j+1] - P[i][j]);
        }
    return;
}

REAL delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz)
{
    REAL duvdz = (U[i+1][j] + U[i+1][j+1])*(V[i][j+1] + V[i+1][j+1]);
    REAL correctionTerm = 0;
    delz *= 4;
    if (z == DERIVE_BY_X)
    {
        duvdz -= (U[i][j] + U[i][j+1])*(V[i-1][j+1] + V[i][j+1]);
        if (alpha == 0)
            return duvdz/delz;
        correctionTerm = abs(U[i+1][j] + U[i+1][j+1]) * (V[i][j+1] - V[i+1][j+1]);
        correctionTerm -= abs(U[i][j]+U[i][j+1]) * (V[i-1][j+1] - V[i][j+1]);
    }
    else
    {
        duvdz -= (U[i+1][j-1] + U[i+1][j])*(V[i][j] + V[i+1][j]);
        if (alpha == 0)
            return duvdz/delz;
        correctionTerm = abs(V[i][j+1] + V[i+1][j+1]) * (U[i+1][j] - U[i+1][j+1]);
        correctionTerm -= abs(V[i][j] + V[i+1][j]) * (U[i+1][j-1] - U[i+1][j]);
    }
    return (duvdz + alpha*correctionTerm)/delz;
}

REAL delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz)
{
    int dx = (z == DERIVE_BY_X) ? 1 : 0;
    int dy = (z == DERIVE_BY_Y) ? 1 : 0;
    delz *= 4;
    REAL df2dz = sqr(F[i][j] + F[i+dx][j+dy]) - sqr(F[i-dx][j-dy] + F[i][j]);
    if (alpha == 0)
        return df2dz/delz;
    REAL correctionTerm = sqr(F[i][j]) - sqr(F[i+dx][j+dy]);
    if (F[i+dx][j+dy] < -F[i][j])
        correctionTerm *= -1;
    df2dz += alpha*correctionTerm;
    correctionTerm = sqr(F[i-dx][j-dy]) - sqr(F[i][j]);
    if (F[i-dx][j-dy] < - F[i][j])
        correctionTerm *= -1;
    df2dz -= alpha*correctionTerm;
    return df2dz/delz;
}

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, char **FLAG, REAL delt,
                lattice *grid, fluidSim *simulation)
{
    /* Compute the auxiliary arrays F and G */
    REAL d2ux, d2uy, d2vx, d2vy;
    REAL du2x, dv2y;
    REAL duvx, duvy;
    short flag;
    int i,j;
    for (i = 1; i <= grid->deli; i++)
    {
        for (j = 1; j <= grid->delj; j++)
        {
            flag = FLAG[i][j];
            if (flag == C_B)      /* Boundary cells with no neighboring fluid cells */
                continue;
            else if (flag == C_F) /* Pure fluid cells */
            {
                d2ux = (U[i+2][j] - 2*U[i+1][j] + U[i][j])/sqr(grid->delx);
                d2uy = (U[i+1][j+1] - 2*U[i+1][j] + U[i+1][j-1])/sqr(grid->dely);

                duvy = delUVbyDelZ(U,V,i,j,DERIVE_BY_Y,simulation->alpha,grid->dely);
                du2x = delFSqrdByDelZ(U,i+1,j,DERIVE_BY_X,simulation->alpha,grid->delx);

                F[i][j] = U[i+1][j] + delt*((d2ux+d2uy)/simulation->Re - du2x - duvy + simulation->GX);

                d2vx = (V[i+1][j+1] - 2*V[i][j+1] + V[i-1][j+1])/sqr(grid->delx);
                d2vy = (V[i][j] - 2*V[i][j+1] + V[i][j])/sqr(grid->dely);

                dv2y = delFSqrdByDelZ(V,i,j+1,DERIVE_BY_Y,simulation->alpha,grid->dely);
                duvx = delUVbyDelZ(U,V,i,j,DERIVE_BY_X,simulation->alpha,grid->delx);

                G[i][j] = V[i][j+1] + delt*((d2vx+d2vy)/simulation->Re - dv2y - duvx + simulation->GY);
            }
            if (flag & B_N)       /* North */
                G[i][j] = V[i][j+1];
            else if (flag & B_S)  /* South */
                G[i-1][j] = V[i-1][j+1];
            if (flag & B_O)       /* East */
                F[i][j] = U[i+1][j];
            else if (flag & B_W)  /* West */
                F[i-1][j] = U[i][j];
        }
    }
    if (grid->edges & LEFT)
    {
        for (j = 1; j<= grid->delj; j++)
            F[0][j] = U[1][j];
    }
    if (grid->edges & RIGHT)
    {
        for (j = 1; j<= grid->delj; j++)
            F[grid->deli][j] = U[grid->deli+1][j];
    }
    if (grid->edges & BOTTOM)
    {
        for (i = 1; i <= grid->deli; i++)
            G[i][0] = V[i][1];
    }
    if (grid->edges & TOP)
    {
        for (i = 1; i <= grid->deli; i++)
            G[i][grid->delj] = V[i][grid->delj+1];
    }
    return;
}

int simulateFluid (REAL **U, REAL **V, REAL **P,
                   boundaryCond* bCond, lattice *grid, fluidSim *sim, MPI_Comm Region,
                   REAL t_end, const char *problem, int opt)
{
    /* Simulate the fluid characterized by sim, grid and bCond */
    /* Error checking */
    if (!U || !V|| !P)
        return 0;
    if (!bCond || !grid || !sim)
        return 0;
    /* Auxiliary Grids: */
    /* RHS is used for the Poisson-Solver so no ghost cells are neccessary */
    REAL **F = create2Dfield(grid->deli+1,grid->delj+1);
    REAL **G = create2Dfield(grid->deli+1,grid->delj+1);
    REAL **RHS = create2Dfield(grid->deli,grid->delj);
    if (!F || !G || !RHS)
        return 0;
    int partcount = 5000, n = 1, rank;
    MPI_Comm_rank(Region,&rank);
    REAL del_vec, dt = 0, delt = sim->dt;
    REAL buf[grid->deli+grid->delj];
    if (opt >= OUTPUT)
        del_vec = t_end/(opt/OUTPUT);
    else
        del_vec = t_end*2;
    /* Create Particles
    particle *parts = createParticleArray(partcount);
    if (!parts)
    {
        printf("Creating particles failed. Proceed? (y/n)");
        if (getchar() == 'n')
            return 0;
    } */
    /* Begin the simulation */
    if (rank == 0)
        printf("Computing Reynoldsnumber %lg.\n",sim->Re);
    for (REAL time = 0; time <= t_end; time += delt)
    {
        if ((rank == 0) && (opt & PRINT))
            printf("Time is at %lg seconds\n",time);
        /* Update all parameters and fields for the iteration */
        setBCond(U,V,grid,bCond);
        setSpecBCond(U,V,grid,problem);
        compFG(U,V,F,G,bCond->FLAG,delt,grid,sim);
        compRHS(F,G,RHS,bCond->FLAG,grid,delt);

        /* Solve the Poisson Equation */
        solveSORforPoisson(P,RHS,bCond->FLAG,sim,grid, Region);
        /* Update U and V through F,G and P */
        adaptUV(U,V,P,F,G,delt,bCond->FLAG,grid);
        exchangeMat(U,2,1,buf,grid,Region);
        exchangeMat(V,1,2,buf,grid,Region);
        dt = compDelt(grid,U,V,sim);

        if (time > del_vec*n)
        {
            dumpFields(Region,U,V,P,grid,n++);
            /*WriteParticle(parts,partcount,n);
            ParticleSeed(parts,0.1,0.9,0.1,0.9,partcount,50);*/
        }
        /*ParticleVelocity(U,V,parts,grid,bCond->FLAG,partcount);
        ParticleTransport(parts,partcount,delt);*/
        MPI_Allreduce(&dt,&delt,1,MPI_DOUBLE,MPI_MIN,Region);
    }
    if (rank == 0)
        printf("[Simulation complete!]\n");

    /* Destroy non-simulated grids */
    destroy2Dfield(F,grid->deli+1);
    destroy2Dfield(G,grid->deli+1);
    destroy2Dfield(RHS,grid->deli);
    /* Destroy structures
    destroyParticleArray(parts); */
    return n;
}
