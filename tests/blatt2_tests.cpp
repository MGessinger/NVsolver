#include <iostream>

#include "gtest/gtest.h"

extern "C" {
    // Fieldlib is a C library and requires C linkage
    #include "blas.h"
    #include "fields.h"
    #include "real.h"
    #include "poisson.h"
}

TEST(Poisson,MatrixCreation)
{
    REAL **A = create2DpoissonMatrix(0.5,1.5,3,2);
    REAL **B = createMatrix(6,6);
    for (int i = 0; i < 6; i++)
        B[i][i] = 136;
    for (int i = 0; i < 5; i++)
        if (i%2 == 0)
        {
            B[i+1][i] = -4;
            B[i][i+1] = -4;
        }
    for (int i = 0; i < 4; i++)
    {
        B[i][i+2] = -64;
        B[i+2][i] = -64;
    }
    EXPECT_TRUE(isEqual2Dfield(A,B,6,6,1e-10));
    destroyMatrix(A,6);
    destroyMatrix(B,6);
}

TEST(Poisson,GaussSolver)
{
    int size = 20;
    REAL **A = create2DpoissonMatrix(1,1,size,size);
    REAL *x = createVector(size*size);
    REAL *b = createVector(size*size);
    fill1Dfield(1,x,size*size);
    solveSOR(A,x,b,size*size,size*size,1.7,1e-10,10000);
    EXPECT_TRUE(isEqual1Dfield(b,x,size*size,1e-10*size));
    destroyMatrix(A,size*size);
    destroyVector(x);
    destroyVector(b);
    return;
}

REAL f(REAL x, REAL y)
{
    return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);
}

REAL f_act(REAL x, REAL y)
{
    return sin(2*M_PI*x)*sin(2*M_PI*y);
}

int imax = 20, jmax = 20;

TEST(Poisson,NaiveSolver)
{
    lattice grid = (lattice){ .imax = imax, .jmax = jmax, .delx = 1./21, .dely = 1./21};
    REAL **A = create2DpoissonMatrix(1,1,imax,jmax);
    REAL *x = createVector(imax*jmax);
    REAL *b = createVector(imax*jmax);
    fill1Dfield(1,x,imax*jmax);
    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
            b[i*jmax+j] = -f((i+1)*grid.delx,(j+1)*grid.dely);
    solveSOR(A,x,b,imax*jmax,imax*jmax,1.7,1e-10,5000);
    REAL **fGrid = sampleFDgridOnCellCorners(f_act,&grid);
    REAL diff = 0;

    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
        {
            REAL DX = x[i*jmax+j]-fGrid[i][j];
            diff += DX*DX;
        }
    diff = sqrt(diff);
    printf("The error for the naive one was %lg.\n",diff);
    EXPECT_LE(diff,grid.delx+grid.dely);
    destroyMatrix(A,imax*jmax);
    destroyMatrix(fGrid,imax);
    destroyVector(x);
    destroyVector(b);
    return;
}

TEST(Poisson,EfficientSolver)
{
    lattice grid = (lattice){ .imax = imax, .jmax = jmax, .delx = 0.05, .dely = 0.05};
    REAL **p = createMatrix(imax+2,jmax+2);
    REAL **rhs = sampleFDgridOnCellCenters(f,&grid);
    REAL **F = sampleFDgridOnCellCenters(f_act,&grid);
    fill2Dfield(1,p,imax+2,jmax+2);
    int it = solveSORforPoisson(p,rhs,1.7,1e-10,5000,0,&grid);
    REAL diff = 0;
    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
        {
            REAL DX = p[i+1][j+1] - F[i][j];
            diff += DX*DX;
        }
    diff = sqrt(diff);
    printf("The efficient algorithm was executed %i times barring an error of %lg.\n",it,diff);
    EXPECT_LE(diff,grid.delx+grid.dely);
    destroyMatrix(p,imax+2);
    destroyMatrix(rhs,imax);
    destroyMatrix(F,imax);
    return;
}
