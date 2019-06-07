#include <iostream>

#include "gtest/gtest.h"

extern "C" {
    // Fieldlib is a C library and requires C linkage
    #include "blas.h"
    #include "fields.h"
    #include "real.h"
    #include "poisson.h"
    #include "IO.h"
    #include "particle.h"
}

TEST(Geometrie,Stufe)
{
    fluidSim fluid;
    lattice grid;
    REAL delt, t_end;
    char problem[128];
    ASSERT_EQ(readParameters((char*)"step.par",NULL,NULL,NULL,&grid,&fluid,&delt,&t_end,problem),17);
    boundaryCond *bCond = createBoundCond(grid.imax,grid.jmax,NOSLIP,OUTFLOW,NOSLIP,NOSLIP);
    initFlags("Step",bCond->FLAG,grid.imax,grid.jmax);
    /*for (int j = grid.jmax-1; j >= 0; j--)
    {
        for (int i = 0; i < grid.imax; i++)
            printf("%x,",bCond->FLAG[i][j]);
        printf("\n");
    }*/
    destroyBoundCond(bCond,grid.imax);
    return;
}

TEST(Particle,IO)
{
    int partcount = 100;
    particle *parts = createParticleArray(partcount);
    ParticleSeed(parts,0,10,0,10,partcount,partcount);
    WriteParticle(parts,partcount,0);
    destroyParticleArray(parts);
    return;
}
