#include <iostream>

#include "gtest/gtest.h"

extern "C" {
    // Fieldlib is a C library and requires C linkage
    #include "blas.h"
    #include "fields.h"
    #include "real.h"
    #include "poisson.h"
    #include "IO.h"
}

TEST(Simulation,Initialisation)
{
    fluidSim fluid;
    lattice grid;
    REAL delx, dely, delt;
    REAL tau, UI, VI, PI;
    /* Must cast to char* explicitly because C++ is a pile of shit */
    int vars = readParameters((char*)"dcavity.par",&grid,&fluid,&delx,&dely,&delt,&tau,&UI,&VI,&PI);
    ASSERT_EQ(vars,17);
    EXPECT_EQ(delx,0.2);
    EXPECT_EQ(dely,0.2);
    REAL **U = createMatrix(grid.imax+2,grid.jmax+2);
    REAL **V = createMatrix(grid.imax+2,grid.jmax+2);
    REAL **P = createMatrix(grid.imax+2,grid.jmax+2);
    initUVP(U,V,P,grid.imax,grid.jmax,UI,VI,PI);
    setBCond(U,V,grid.imax,grid.jmax,NULL);
    EXPECT_TRUE(isEqual1Dfield(V[0],V[grid.imax+1],grid.jmax+2,1e-10));
    compDelt(&delt,grid.imax,grid.jmax,delx,dely,U,V,fluid.Re,tau);
    EXPECT_LE(delt-0.2,1e-10);
    destroyMatrix(U,grid.imax+2);
    destroyMatrix(V,grid.imax+2);
    destroyMatrix(P,grid.imax+2);
    return;
}

TEST(Simulation,TrivialFluid)
{
    REAL **U, **V, **P;
    lattice *grid = simulateFluid(&U,&V,&P,"empty.par", SILENT,NULL);
    REAL **zero = createMatrix(grid->imax+2,grid->jmax+2);
    EXPECT_TRUE(isEqual2Dfield(U,zero,grid->imax+2,grid->jmax+2,1e-3));

    destroyMatrix(U,grid->imax+2);
    destroyMatrix(V,grid->imax+2);
    destroyMatrix(P,grid->imax+2);
    destroyMatrix(zero,grid->imax+2);
    free(grid);
    return;
}
