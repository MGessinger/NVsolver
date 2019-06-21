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
    REAL delt, t_end;
    /* Must cast to char* explicitly because C++ is a pile of shit */
    char problem[128];
    REAL **U, **V, **P;
    int vars = readParameters((char*)"dcavity.par",&U,&V,&P,&grid,&fluid,NULL,&delt,&t_end,problem);
    ASSERT_EQ(vars,17);
    EXPECT_EQ(grid.delx,0.2);
    EXPECT_EQ(grid.dely,0.2);
    setBCond(U,V,grid.imax,grid.jmax,NULL);
    EXPECT_TRUE(isEqual1Dfield(V[0],V[grid.imax+1],grid.jmax+2,1e-10));
    compDelt(&delt,&grid,U,V,fluid.Re,fluid.tau);
    EXPECT_LE(delt-0.2,1e-10);
    destroyMatrix(U,grid.imax+2);
    destroyMatrix(V,grid.imax+2);
    destroyMatrix(P,grid.imax+2);
    return;
}

TEST(Simulation,TrivialFluid)
{
    REAL **U, **V, **P;
    int imax = 32, jmax = 32;
    boundaryCond *bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
    lattice *grid = simulateFluid(&U,&V,&P,"empty.par",SILENT);
    outputVec(U,V,P,NULL,grid,0,0);
    REAL **zero = createMatrix(grid->imax+2,grid->jmax+2);
    EXPECT_TRUE(isEqual2Dfield(U,zero,grid->imax+2,grid->jmax+2,1e-3));

    destroyMatrix(U,grid->imax+2);
    destroyMatrix(V,grid->imax+2);
    destroyMatrix(P,grid->imax+2);
    destroyMatrix(zero,grid->imax+2);
    destroyBoundCond(bCond,grid->imax);
    free(grid);
    return;
}

TEST(Simulation,Tunnel)
{
    REAL **U, **V, **P;
    int imax = 32, jmax = 32;
    boundaryCond *bCond = createBoundCond(NOSLIP,OUTFLOW,FREESLIP,FREESLIP);
    lattice *grid = simulateFluid(&U,&V,&P,"tunnel.par",SILENT);
    outputVec(U,V,P,NULL,grid,0,0);
    REAL **expected = createMatrix(grid->imax,grid->jmax);
    fill2Dfield(1,expected,grid->imax,grid->jmax);
    EXPECT_TRUE(isEqual2Dfield(U+1,expected,grid->imax,grid->jmax,1e-3));

    destroyMatrix(U,grid->imax+2);
    destroyMatrix(V,grid->imax+2);
    destroyMatrix(P,grid->imax+2);
    destroyMatrix(expected,grid->imax);
    destroyBoundCond(bCond,grid->imax);
    free(grid);
    return;
}
