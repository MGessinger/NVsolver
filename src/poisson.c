#include "poisson.h"
/*
REAL** create2DpoissonMatrix(REAL ilength, REAL jlength, int imax, int jmax)
{
    int size = imax*jmax;
    REAL DX = (imax+1)/ilength;
    REAL DY = (jmax+1)/jlength;
    REAL **A = createMatrix(size,size);
    REAL **B = createMatrix(size,size);
    if (A == NULL || B == NULL)
        return NULL;
    fill2Dfield(0,A,size,size);
    fill2Dfield(0,B,size,size);
    for (int i = 0; i < size; i++)
    {
        A[i][i] = 2;
        if (i < size-jmax)
        {
            A[i][i+jmax] = -1;
            A[i+jmax][i] = -1;
        }
        B[i][i] = 2;
        if ((i+1)%jmax != 0)
        {
            B[i][i+1] = -1;
            B[i+1][i] = -1;
        }
    }
    scal2Dfield(DY*DY,B,size,size);
    axpy2Dfield(DX*DX,A,B,size,size);
    destroyMatrix(A,size);
    return B;
}

int solveSOR(REAL **A, REAL *x, REAL *b, int rows, int cols, REAL omega, REAL epsilon, int itermax)
{
    if (rows != cols)
    {
        return 0;
    }
    if (A == NULL || x == NULL || b == NULL)
    {
        return 0;
    }
    REAL DeltaX, error;
    int it = 0;
    do {
        error = 0;
        for (int i = 0; i < cols; i++)
        {
            DeltaX = omega*(b[i] - dot(A[i],x,cols))/A[i][i];
            x[i] += DeltaX;
            error += DeltaX*DeltaX;
            if (error >= 1e100)
            {
                fill1Dfield(0,x,cols);
                return solveSOR(A,x,b,rows,cols,1,epsilon,itermax);
            }
        }
        it++;
    } while (sqrt(error) > epsilon && it < itermax);
    return it;
}*/

void applyPboundaryCond(REAL **P, int imax, int jmax, REAL dxSqrd, REAL dySqrd, short **FLAG)
{
    if (P == NULL || FLAG == NULL)
        return;
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            if (FLAG[i-1][j-1] == C_F)
                continue;
            else if (FLAG[i-1][j-1] == C_B)
            {
                P[i][j] = 0;
                continue;
            }
            switch (FLAG[i-1][j-1] - C_B)
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

int solveSORforPoisson(REAL **p, REAL **rhs, short **FLAG, REAL omega, REAL epsilon,
                        int itermax, lattice *grid)
{
    if (p == NULL || rhs == NULL || grid == NULL || FLAG == NULL)
        return 0;
    int it = 0;
    REAL temporary, error;
    REAL invXWidthSqrd = 1/(grid->delx*grid->delx);
    REAL invYWidthSqrd = 1/(grid->dely*grid->dely);
    REAL eps = epsilon*epsilon;
    int numberOfCells;
    do {
        numberOfCells = grid->imax*grid->jmax;
        /* Apply the boundary condition */
        applyHomogeneousNeumannBC(p,grid->imax,grid->jmax);
        applyPboundaryCond(p,grid->imax,grid->jmax,invXWidthSqrd,invYWidthSqrd,FLAG);
        /* Compute the new coefficients iteratively */
        for (int i = 1; i <= grid->imax; i++)
            for (int j = 1; j <= grid->jmax; j++)
            {
                if (FLAG[i-1][j-1] != C_F)
                {
                    numberOfCells--;
                    continue;
                }
                temporary = invXWidthSqrd*(p[i+1][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] + p[i][j-1]) - rhs[i-1][j-1];
                p[i][j] = omega*temporary/(2*invXWidthSqrd+2*invYWidthSqrd) + (1-omega)*p[i][j];
            }
        error = 0;
        /* Calculate the residue with respect to the rhs field */
        for (int i = 1; i <= grid->imax; i++)
            for (int j = 1; j <= grid->jmax; j++)
            {
                if (FLAG[i-1][j-1] != C_F)
                    continue;
                temporary = invXWidthSqrd*(p[i+1][j] - 2*p[i][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] - 2*p[i][j] + p[i][j-1]);
                temporary -= rhs[i-1][j-1];
                error += temporary*temporary;
            }
        if (++it >= itermax)
            break;
    } while (error/numberOfCells > eps);
    if (it != 0)
        printf("Remaining error after %i iterations: %e/%i vs. %e\n",it,error,numberOfCells,eps);
    return it;
}

REAL sqr(REAL x)
{
    return x*x;
}

void compDelt(REAL *delt, lattice *grid, REAL **U, REAL **V, REAL Re, REAL tau)
{
    if (tau <= 0)
        return;
    REAL dt = Re/(2/(sqr(grid->delx) + sqr(grid->dely)));
    for (int i = 1; i <= grid->imax; i++)
        for (int j = 1; j <= grid->jmax; j++)
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
    *delt = tau*dt;
    return;
}

void compRHS(REAL **F, REAL **G, REAL **RHS, short **FLAG, lattice *grid, REAL delt)
{
    if (F == NULL || G == NULL || RHS == NULL)
    {
        fill2Dfield(INFINITY,RHS,grid->imax,grid->jmax);
        return;
    }
    for (int i = 1; i <= grid->imax; i++)
        for (int j = 1; j <= grid->jmax; j++)
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
    for (int i = 1; i <= grid->imax; i++)
        for (int j = 1; j <= grid->jmax; j++)
        {
            if (FLAG[i-1][j-1] != C_F)
            {
                U[i][j] = V[i][j] = 0;
                continue;
            }
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
    REAL d2ux, d2uy, d2vx, d2vy;
    REAL du2x, dv2y;
    REAL duvx, duvy;
    short flag;
    for (int i = 1; i <= grid->imax; i++)
    {
        for (int j = 1; j <= grid->jmax; j++)
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
            else    /* Pure fluid cells */
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
    for (int j = 1; j<= grid->jmax; j++)
    {
        F[0][j] = U[0][j];
        F[grid->imax][j] = U[grid->imax][j];
    }
    for (int i = 1; i <= grid->imax; i++)
    {
        G[i][0] = V[i][0];
        G[i][grid->jmax] = V[i][grid->jmax];
    }
    return;
}

lattice* simulateFluid (REAL ***U, REAL ***V, REAL ***P, const char *fileName, int opt)
{
    lattice *grid = malloc(sizeof(lattice));
    if (grid == NULL)
        return NULL;
    grid->imax = 0;
    grid->jmax = 0;
    boundaryCond *bCond;

    int n = 0;
    int partcount = 1000;
    fluidSim simulation;
    REAL delt, t_end;
    char problem[128];
    /* Read prameters from a file */
    if (readParameters(fileName,U,V,P,grid,&simulation,&bCond,&delt,&t_end,problem) < 17)
    {
        destroyBoundCond(bCond,grid->imax);
        free(grid);
        return NULL;
    }
    if (bCond->FLAG == NULL)
    {
        bCond->FLAG = create2DIntegerField(grid->imax,grid->jmax);
        initFlags(problem,bCond->FLAG,grid->imax,grid->jmax);
    }
    REAL del_vec;
    if (opt >= OUTPUT)
        del_vec = t_end/(opt/OUTPUT);
    else
        del_vec = t_end*2;
    printf("Computing Reynoldsnumber %lg.\n",simulation.Re);
    particle *parts = createParticleArray(partcount);

    /* Helping Grids: */
    REAL **F = createMatrix(grid->imax+1,grid->jmax+1);
    REAL **G = createMatrix(grid->imax+1,grid->jmax+1);
    /* RHS is used for the Poisson-Solver so no ghost cells are neccessary */
    REAL **RHS = createMatrix(grid->imax,grid->jmax);
    /* Create Particles */

    /* Error checking */
    if (U == NULL || V == NULL || P == NULL)
        return grid;
    if (F == NULL || G == NULL || RHS == NULL)
        return grid;
    if (parts == NULL)
        return grid;

    /* Begin the simulation */
    for (REAL time = 0; time <= t_end; time += delt)
    {
        if (opt & PRINT)
            printf("Time is at %lg seconds\n",time);
        /* Update all parameters and fields for the iteration */
        setBCond(*U,*V,grid->imax,grid->jmax,bCond);
        setSpecBCond(*U,*V,grid->imax,grid->jmax,problem);
        compFG(*U,*V,F,G,bCond->FLAG,delt,grid,&simulation);
        compRHS(F,G,RHS,bCond->FLAG,grid,delt);

        /* Solve the Poisson Equation */
        solveSORforPoisson(*P,RHS,bCond->FLAG,simulation.omega,simulation.eps,simulation.itmax,grid);
        /* Update U and V through F,G and P */
        adaptUV(*U,*V,*P,F,G,delt,bCond->FLAG,grid);
        compDelt(&delt,grid,*U,*V,simulation.Re,simulation.tau);

        if (time > del_vec*n)
        {
            outputVec(*U,*V,*P,parts,grid,partcount,++n);
            ParticleSeed(parts,0.1,0.1,0.7,1.3,partcount,50);
        }
        ParticleVelocity(*U,*V,parts,grid,bCond->FLAG,partcount);
        ParticleTransport(parts,partcount,delt);
    }
    printf("[Simulation complete!]\n");
    writeVTKfileFor2DintegerField("GeometryField.vtk","geometryfield",bCond->FLAG,grid->imax,grid->jmax,grid->delx,grid->dely);

    /* Destroy non-simulated grids */
    destroyMatrix(F,grid->imax+1);
    destroyMatrix(G,grid->imax+1);
    destroyMatrix(RHS,grid->imax);
    /* Destroy structures */
    destroyParticleArray(parts);
    destroyBoundCond(bCond,grid->imax);
    return grid;
}
