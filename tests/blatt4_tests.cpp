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

TEST(Geometrie,Stufe)
{
    fluidSim fluid;
    lattice grid;
    REAL delt, t_end;
    REAL UI, VI, PI;
    char problem[128];
    ASSERT_EQ(readParameters((char*)"step.par",&grid,&fluid,&delt,&t_end,&UI,&VI,&PI,problem),17);
    boundaryCond *bCond = createBoundCond(grid.imax,grid.jmax,NOSLIP,OUTFLOW,NOSLIP,NOSLIP);
    initFlags("Step",bCond->FLAG,grid.imax,grid.jmax);
    for (int j = grid.jmax-1; j >= 0; j--)
    {
        for (int i = 0; i < grid.imax; i++)
            printf("%x,",bCond->FLAG[i][j]);
        printf("\n");
    }
    destroyBoundCond(bCond,grid.imax);
    return;
}
