#include "poisson.h"

void applyPboundaryCond(REAL **P, lattice *grid, short **FLAG)
{
    /* Apply boundary conditions for the pressure field */
    if (P == NULL || FLAG == NULL)
        return;
    short flag = 0;
    /* First set values on the actual boundary of the region */
    if (grid->il == 0)
        for (int j = grid->jb; j <= grid->jt; j++)
            P[0][j] = P[1][j];
    if (grid->ir == grid->imax)
        for (int j = grid->jb; j <= grid->jt; j++)
            P[grid->ir-grid->il+1][j] = P[grid->ir-grid->il][j];
    if (grid->jb == 0)
        for (int i = grid->il; i <= grid->ir; i++)
            P[i][0] = P[i][1];
    if (grid->jt == grid->jmax)
        for (int i = grid->il; i <= grid->ir; i++)
            P[i][grid->jt-grid->jb+1] = P[i][grid->jt-grid->jb];
    REAL dxSqrd = grid->delx*grid->delx;
    REAL dySqrd = grid->dely*grid->dely;
    for (int i = grid->il+1; i <= grid->ir; i++)
        for (int j = grid->jb+1; j <= grid->jt; j++)
        {
            flag = FLAG[i-1][j-1];
            if (flag == C_F)
                continue;
            else if (flag == C_B)
            {
                continue;
            }
            switch (flag - C_B)
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

static inline REAL sqr(REAL x)
{
    return x*x;
}

int solveSORforPoisson(REAL **p, REAL **rhs, short **FLAG,
                       fluidSim *sim, lattice *grid)
{
    /* Use a SOR algorithm to solve the poisson equation */
    if (p == NULL || rhs == NULL || grid == NULL || FLAG == NULL)
        return 0;
    int it = 0;
    REAL temporary, error;
    REAL invXWidthSqrd = 1/(grid->delx*grid->delx);
    REAL invYWidthSqrd = 1/(grid->dely*grid->dely);
    REAL scale = sim->omega/(2*(sqr(1/grid->delx)+sqr(1/grid->dely)));
    REAL eps = sim->eps*sim->eps;
    int i,j;
    /* Count the number of fluid cells */
    int numberOfCells = 0;
    for (i = grid->il+1; i <= grid->ir; i++)
        for (j = grid->jb+1; j <= grid->jt; j++)
        {
            if (FLAG[i-1][j-1] == C_F)
                numberOfCells++;
        }
    do {
        /* Apply the boundary condition */
        applyPboundaryCond(p,grid,FLAG);
        /* Compute the new coefficients iteratively */
        for (i = grid->il+1; i <= grid->ir; i++)
            for (j = grid->jb+1; j <= grid->jt; j++)
            {
                if (FLAG[i-1][j-1] != C_F)
                    continue;
                temporary = invXWidthSqrd*(p[i+1][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] + p[i][j-1]) - rhs[i-1][j-1];
                p[i][j] = scale*temporary + (1-sim->omega)*p[i][j];
            }
        error = 0;
        /* Calculate the residue with respect to the rhs field */
        for (i = grid->il+1; i <= grid->ir; i++)
            for (j = grid->jb+1; j <= grid->jt; j++)
            {
                if (FLAG[i-1][j-1] != C_F)
                    continue;
                temporary = invXWidthSqrd*(p[i+1][j] - 2*p[i][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] - 2*p[i][j] + p[i][j-1]);
                temporary -= rhs[i-1][j-1];
                error += temporary*temporary;
            }
        if (++it >= sim->itmax)
            break;
    } while (error/numberOfCells > eps);
    if (it != 0)
        printf("Remaining error after %i iterations: %e/%i vs. %e\n",it,error,numberOfCells,eps);
    return it;
}

void compDelt(REAL *delt, lattice *grid, REAL **U, REAL **V, fluidSim *sim)
{
    /* Find the optimal step width in time */
    if (sim->tau <= 0)
        return;
    REAL dt = sim->Re/(2/(sqr(grid->delx) + sqr(grid->dely)));
    for (int i = grid->il+1; i <= grid->ir; i++)
        for (int j = grid->jb+1; j <= grid->jt; j++)
        {
            /* Once again, one of the inequalities is trivial */
            if (U[i][j] > 0 && grid->delx/U[i][j] < dt)
                dt = grid->delx/U[i][j];
            else if (U[i][j] < 0 && grid->delx/U[i][j] > -dt)
                dt = -(grid->delx)/U[i][j];
            if (V[i][j] > 0 && grid->dely/V[i][j] < dt)
                dt = (grid->dely)/V[i][j];
            else if (V[i][j] < 0 && grid->dely/V[i][j] > -dt)
                dt = -(grid->dely)/V[i][j];
        }
    *delt = sim->tau*dt;
    if (*delt > sim->dt)
        *delt = sim->dt;
    return;
}

void compRHS(REAL **F, REAL **G, REAL **RHS, short **FLAG, lattice *grid, REAL delt)
{
    if (F == NULL || G == NULL || RHS == NULL)
    {
        fill2Dfield(INFINITY,RHS,grid->imax,grid->jmax);
        return;
    }
    for (int i = grid->il+1; i <= grid->ir; i++)
        for (int j = grid->jb+1; j <= grid->jt; j++)
        {
            if (FLAG[i-1][j-1] != C_F)
                continue;
            RHS[i-1][j-1] = (F[i][j] - F[i-1][j])/(delt*(grid->delx));
            RHS[i-1][j-1] += (G[i][j] - G[i][j-1])/(delt*(grid->dely));
        }
    return;
}

void adaptUV(REAL **U, REAL **V, REAL **P, REAL **F, REAL **G,
             REAL delt, short **FLAG, lattice *grid)
{
    for (int i = grid->il+1; i <= grid->ir; i++)
        for (int j = grid->jb+1; j <= grid->jt; j++)
        {
            if (FLAG[i-1][j-1] != C_F)
                continue;
            U[i][j] = F[i][j] - delt/grid->delx*(P[i+1][j] - P[i][j]);
            V[i][j] = G[i][j] - delt/grid->dely*(P[i][j+1] - P[i][j]);
        }
    return;
}

REAL delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz)
{
    REAL duvdz = (U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j]);
    REAL correctionTerm = 0;
    if (z == DERIVE_BY_X)
    {
        duvdz -= (U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j]);
        if (alpha == 0)
            return duvdz/(4*delz);
        correctionTerm = abs(U[i][j] + U[i][j+1]) * (V[i][j] - V[i+1][j]);
        correctionTerm -= abs(U[i-1][j]+U[i-1][j+1]) * (V[i-1][j] - V[i][j]);
    }
    else
    {
        duvdz -= (U[i][j-1] + U[i][j])*(V[i][j-1] + V[i+1][j-1]);
        if (alpha == 0)
            return duvdz/(4*delz);
        correctionTerm = abs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]);
        correctionTerm -= abs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]);
    }
    return (duvdz + alpha*correctionTerm)/(4*delz);
}

REAL delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz)
{
    int dx = (z == DERIVE_BY_X) ? 1 : 0;
    int dy = (z == DERIVE_BY_Y) ? 1 : 0;
    REAL df2dz = sqr(F[i][j] + F[i+dx][j+dy]) - sqr(F[i-dx][j-dy] + F[i][j]);
    if (alpha == 0)
        return df2dz/(4*delz);
    REAL correctionTerm = sqr(F[i][j]) - sqr(F[i+dx][j+dy]);
    if (F[i+dx][j+dy] < -F[i][j])
        correctionTerm *= -1;
    df2dz += alpha*correctionTerm;
    correctionTerm = sqr(F[i-dx][j-dy]) - sqr(F[i][j]);
    if (F[i-dx][j-dy] < - F[i][j])
        correctionTerm *= -1;
    df2dz -= alpha*correctionTerm;
    return df2dz/(4*delz);
}

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, short **FLAG, REAL delt,
                lattice *grid, fluidSim *simulation)
{
    /* Compute the auxiliary arrays F and G */
    REAL d2ux, d2uy, d2vx, d2vy;
    REAL du2x, dv2y;
    REAL duvx, duvy;
    short flag;
    int i,j;
    for (i = grid->il+1; i <= grid->ir; i++)
    {
        for (j = grid->jb+1; j <= grid->jt; j++)
        {
            flag = FLAG[i-1][j-1];
            if (flag == C_B) /* Boundary cells with no neighboring fluid cells */
                continue;
            if (flag & B_N)     /* North */
                G[i][j] = V[i][j];
            if (flag & B_O)     /* East */
                F[i][j] = U[i][j];
            if (flag & B_S)     /* South */
                G[i-1][j] = V[i-1][j];
            if (flag & B_W)     /* West */
                F[i-1][j] = U[i-1][j];
            if (flag == C_F)    /* Pure fluid cells */
            {
                d2ux = (U[i+1][j] - 2*U[i][j] + U[i-1][j])/sqr(grid->delx);
                d2uy = (U[i][j+1] - 2*U[i][j] + U[i][j-1])/sqr(grid->dely);

                duvy = delUVbyDelZ(U,V,i,j,DERIVE_BY_Y,simulation->alpha,grid->dely);
                du2x = delFSqrdByDelZ(U,i,j,DERIVE_BY_X,simulation->alpha,grid->delx);

                F[i][j] = U[i][j] + delt*((d2ux+d2uy)/simulation->Re - du2x - duvy + simulation->GX);

                d2vx = (V[i+1][j] - 2*V[i][j] + V[i-1][j])/sqr(grid->delx);
                d2vy = (V[i][j+1] - 2*V[i][j] + V[i][j-1])/sqr(grid->dely);

                dv2y = delFSqrdByDelZ(V,i,j,DERIVE_BY_Y,simulation->alpha,grid->dely);
                duvx = delUVbyDelZ(U,V,i,j,DERIVE_BY_X,simulation->alpha,grid->delx);

                G[i][j] = V[i][j] + delt*((d2vx+d2vy)/simulation->Re - dv2y - duvx + simulation->GY);
            }
        }
    }
    if (grid->il == 0)
    {
        for (j = grid->jb+1; j<= grid->jt; j++)
            F[0][j] = U[0][j];
    }
    if (grid->ir == grid->imax)
    {
        for (j = grid->jb+1; j<= grid->jt; j++)
            F[grid->ir-grid->il][j] = U[grid->ir-grid->il][j];
    }
    if (grid->jb == 0)
    {
        for (i = grid->il+1; i <= grid->ir; i++)
            G[i][0] = V[i][0];
    }
    if (grid->jt == grid->jmax)
    {
        for (i = grid->il+1; i <= grid->ir; i++)
            G[i][grid->jt-grid->jb] = V[i][grid->jt-grid->jb];
    }
    return;
}

int simulateFluid (REAL **U, REAL **V, REAL **P,
                        boundaryCond* bCond, lattice *grid, fluidSim *sim,
                        REAL delt, REAL t_end, const char *problem, int opt)
{
    /* Simulate the fluid characterized by sim, grid and bCond */
    /* Error checking */
    if (!U || !V|| !P)
        return 0;
    /* Auxiliary Grids: */
    REAL **F = create2Dfield(grid->imax+1,grid->jmax+1);
    REAL **G = create2Dfield(grid->imax+1,grid->jmax+1);
    /* RHS is used for the Poisson-Solver so no ghost cells are neccessary */
    REAL **RHS = create2Dfield(grid->imax,grid->jmax);
    if (!F || !G || !RHS)
        return 0;
    int partcount = 5000, n = 0;
    REAL del_vec;
    if (opt >= OUTPUT)
        del_vec = t_end/(opt/OUTPUT);
    else
        del_vec = t_end*2;
    /* Create Particles */
    particle *parts = createParticleArray(partcount);
    if (!parts)
    {
        printf("Creating particles failed. Proceed? (y/n)");
        if (getchar() == 'n')
            return 0;
    }
    /* Begin the simulation */
    printf("Computing Reynoldsnumber %lg.\n",sim->Re);
    for (REAL time = 0; time <= t_end; time += delt)
    {
        if (opt & PRINT)
            printf("Time is at %lg seconds\n",time);
        /* Update all parameters and fields for the iteration */
        setBCond(U,V,grid,bCond);
        setSpecBCond(U,V,grid,problem);
        compFG(U,V,F,G,bCond->FLAG,delt,grid,sim);
        compRHS(F,G,RHS,bCond->FLAG,grid,delt);

        /* Solve the Poisson Equation */
        solveSORforPoisson(P,RHS,bCond->FLAG,sim,grid);
        /* Update U and V through F,G and P */
        adaptUV(U,V,P,F,G,delt,bCond->FLAG,grid);
        compDelt(&delt,grid,U,V,sim);

        if (time > del_vec*n)
        {
            outputVec(U,V,P,grid,++n);
            WriteParticle(parts,partcount,n);
            ParticleSeed(parts,0.1,0.9,0.1,0.9,partcount,50);
        }
        ParticleVelocity(U,V,parts,grid,bCond->FLAG,partcount);
        ParticleTransport(parts,partcount,delt);
    }
    printf("[Simulation complete!]\n");
    writeVTKfileFor2DintegerField("GeometryField.vtk","geometryfield",bCond->FLAG,grid);

    /* Destroy non-simulated grids */
    destroy2Dfield(F,grid->imax+1);
    destroy2Dfield(G,grid->imax+1);
    destroy2Dfield(RHS,grid->imax);
    /* Destroy structures */
    destroyParticleArray(parts);
    return 1;
}
