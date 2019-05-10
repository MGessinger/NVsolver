#include "poisson.h"

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
    /* Explicitly set-up the Poisson-matrix for the dicrete lattice */
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
    /* Use the relaxed Gauss-Seidel-Algorithm to solve the equation Ax=b */
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
}

int solveSORforPoisson(REAL **p, REAL **rhs, REAL omega, REAL epsilon,
                        int itermax, int useNeumann, lattice *grid)
{
    if (p == NULL || rhs == NULL)
    {
        return 0;
    }
    int it = 0;
    REAL temporary, error;
    REAL invXWidthSqrd = (grid->imax*grid->imax)/(grid->xlength*grid->xlength); /* invXWidthSqrd = 1/(dx)^2 */
    REAL invYWidthSqrd = (grid->jmax*grid->jmax)/(grid->ylength*grid->ylength); /* invYWidthSqrd = 1/(dy)^2 */
    REAL eps2 = epsilon*epsilon;
    do {
        /* Apply the boundary condition */
        if (useNeumann == 1)
            applyHomogeneousNeumannBC(p,grid->imax,grid->jmax);
        else
            applyHomogeneousDirichletBC(p,grid->imax,grid->jmax);
        /* Compute the new coefficients iteratively */
        for (int i = 1; i <= grid->imax; i++)
            for (int j = 1; j <= grid->jmax; j++)
            {
                temporary = invXWidthSqrd*(p[i+1][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] + p[i][j-1]) - rhs[i-1][j-1];
                temporary *= omega/(2*invXWidthSqrd+2*invYWidthSqrd);
                p[i][j] = temporary + (1-omega)*p[i][j];
            }
        error = 0;
        /* Calculate the residue with respect to the rhs field */
        for (int i = 1; i <= grid->imax; i++)
            for (int j = 1; j <= grid->jmax; j++)
            {
                temporary = invXWidthSqrd*(p[i+1][j] - 2*p[i][j] + p[i-1][j]);
                temporary += invYWidthSqrd*(p[i][j+1] - 2*p[i][j] + p[i][j-1]);
                temporary -= rhs[i-1][j-1];
                error += temporary*temporary;
            }
        error = error/(grid->imax*grid->jmax);
        it++;
    } while (error > eps2 && it < itermax);
    return it;
}

void compDelt(REAL *delt, REAL imax, REAL jmax, REAL delx, REAL dely, REAL **U, REAL **V, REAL Re, REAL tau)
{
    if (tau <= 0)
        return;
    REAL dt = Re/(2/(delx*delx + dely*dely));
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            /* Once again, one of the inequalities is trivial */
            if (U[i][j] > 0 && delx/U[i][j] < dt)
                dt = delx/U[i][j];
            else if (U[i][j] < 0 && delx/U[i][j] > -dt)
                dt = -delx/U[i][j];
            if (V[i][j] > 0 && dely/V[i][j] < dt)
                dt = dely/V[i][j];
            else if (V[i][j] < 0 && dely/V[i][j] > -dt)
                dt = -dely/V[i][j];
        }
    *delt = tau*dt;
    return;
}

void compRHS(REAL **F, REAL **G, REAL **RHS, int imax, int jmax, REAL delx, REAL dely, REAL delt)
{
    if (F == NULL || G == NULL || RHS == NULL)
    {
        fill2Dfield(INFINITY,RHS,imax,jmax);
        return;
    }
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            RHS[i-1][j-1] = (F[i][j] - F[i-1][j])/(delt*delx);
            RHS[i-1][j-1] += (G[i][j] - G[i][j-1])/(dely*delt);
        }
    return;
}

void adaptUV(REAL **U, REAL **V, REAL **P, REAL **F, REAL **G,
             REAL delt, REAL delx, REAL dely, int imax, int jmax)
{
    for (int i = 1; i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            U[i][j] = F[i][j] - delt/delx*(P[i+1][j] - P[i][j]);
            V[i][j] = G[i][j] - delt/dely*(P[i][j+1] - P[i][j]);
        }
    return;
}

REAL sqr(REAL x)
{
    return x*x;
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

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, int imax, int jmax,
                REAL delt, REAL delx, REAL dely, fluidSim *simulation)
{
    REAL d2ux, d2uy, d2vx, d2vy;
    REAL du2x, dv2y;
    REAL duvx, duvy;
    copy(U[0],F[0],imax+1);
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            if (i != imax)
            {
                d2ux = (U[i+1][j] - 2*U[i][j] + U[i-1][j])/sqr(delx);
                d2uy = (U[i][j+1] - 2*U[i][j] + U[i][j-1])/sqr(dely);
                duvy = delUVbyDelZ(U,V,i,j,DERIVE_BY_Y,simulation->alpha,dely);

                du2x = delFSqrdByDelZ(U,i,j,DERIVE_BY_X,simulation->alpha,delx);

                F[i][j] = U[i][j] + delt*((d2ux+d2uy)/simulation->Re - du2x - duvy + simulation->GX);
            }
            else
                F[i][j] = U[i][j];
            if (j == jmax) /* G will not be updated for j == jmax */
                continue;
            dv2y = delFSqrdByDelZ(V,i,j,DERIVE_BY_Y,simulation->alpha,dely);
            d2vx = (V[i+1][j] - 2*V[i][j] + V[i-1][j])/sqr(delx);
            d2vy = (V[i][j+1] - 2*V[i][j] + V[i][j-1])/sqr(dely);
            duvx = delUVbyDelZ(U,V,i,j,DERIVE_BY_X,simulation->alpha,delx);

            G[i][j] = V[i][j] + delt*((d2vx+d2vy)/simulation->Re - dv2y - duvx + simulation->GY);
        }
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];
    }
    return;
}

lattice* simulateFluid (REAL ***U, REAL ***V, REAL ***P, const char *fileName, int opt, boundaryCond *bCond)
{
    int n = 1;
    lattice *grid = malloc(sizeof(lattice));
    if (grid == NULL)
        return NULL;
    fluidSim simulation;
    REAL delx, dely, delt;
    REAL tau, UI, VI, PI;
    char problem[128];
    /* Read prameters from a file */
    if (readParameters(fileName,grid,&simulation,&delx,&dely,&delt,&tau,&UI,&VI,&PI,problem) < 17)
        return grid;
    REAL del_vec;
    if (opt >= OUTPUT)
        del_vec = simulation.t_end/(opt/OUTPUT);
    else
        del_vec = simulation.t_end*2;
    printf("Computing Reynoldsnumber %lg.\n",simulation.Re);

    /* Simulated Grids: */
    *P = createMatrix(grid->imax+2,grid->jmax+2);
    *U = createMatrix(grid->imax+2,grid->jmax+2);
    *V = createMatrix(grid->imax+2,grid->jmax+2);
    /* Helping Grids: */
    REAL **F = createMatrix(grid->imax+1,grid->jmax+1);
    REAL **G = createMatrix(grid->imax+1,grid->jmax+1);
    /* RHS is used for the *Poisson-Solver so no ghost cells are neccessary */
    REAL **RHS = createMatrix(grid->imax,grid->jmax);

    /* Error checking */
    if (U == NULL || V == NULL || P == NULL)
        return grid;
    if (F == NULL || G == NULL || RHS == NULL)
        return grid;
    /* Initialise the algorithm */
    initUVP(*U,*V,*P,grid->imax,grid->jmax,UI,VI,PI);
    for (REAL time = 0; time <= simulation.t_end; time += delt)
    {
        /* Update all parameters and fields for the iteration */
        compDelt(&delt,grid->imax,grid->jmax,delx,dely,*U,*V,simulation.Re,tau);
        setBCond(*U,*V,grid->imax,grid->jmax,bCond);
        setSpecBCond(*U,*V,grid->imax,grid->jmax,"Driven Cavity");
        compFG(*U,*V,F,G,grid->imax,grid->jmax,delt,delx,dely,&simulation);
        compRHS(F,G,RHS,grid->imax,grid->jmax,delx,dely,delt);

        /* Solve the Poisson Equation */
        solveSORforPoisson(*P,RHS,simulation.omega,simulation.eps,simulation.itmax,1,grid);
        /* Update U and V through F,G and P */
        adaptUV(*U,*V,*P,F,G,delt,delx,dely,grid->imax,grid->jmax);

        if (time > del_vec*n)
        {
            outputVec(*U,*V,*P,grid,n);
            n++;
        }
        if (opt & PRINT)
            printf("Time is at %lg\n",time);
        if (delt <= 1e-5)
        {
            printf("Aborting loop. Please verify the output for correctness.\n");
            break;
        }
    }
    printf("[Simulation complete!]\n");
    /* Destroy non-simulated grids */
    destroyMatrix(F,grid->imax+1);
    destroyMatrix(G,grid->imax+1);
    destroyMatrix(RHS,grid->imax);
    return grid;
}
