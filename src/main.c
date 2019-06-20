#include <stdlib.h>
#include <stdio.h>

#include "fields.h"
#include "real.h"
#include "poisson.h"
#include "IO.h"
#include "boundary.h"
#include "particle.h"
#include "geometry.h"

int main (int argc, char **argv)
{
    int minimiumHeight = 32, minimumWidth = 32;
    readGeometry(argv[1],minimiumHeight,minimumWidth);
    return 0;


    REAL **U = NULL, **V = NULL, **P = NULL;
    int out = (argc >= 3) ? atoi(argv[2]) : 10;
    boundaryCond *bCond = createBoundCond(0,0,NOSLIP,OUTFLOW,NOSLIP,NOSLIP);
    lattice *grid = simulateFluid(&U,&V,&P,(argc >= 2 ? argv[1] : "dcavity.par"),PRINT | out*OUTPUT, bCond);
    outputVec(U,V,P,NULL,grid,0,0);
    /* Destroy simulated grids */
    if (grid != NULL)
    {
        destroyMatrix(U,grid->imax+2);
        destroyMatrix(V,grid->imax+2);
        destroyMatrix(P,grid->imax+2);
        destroyBoundCond(bCond,grid->imax);
        free(grid);
    }
    else
    {
        destroyBoundCond(bCond,0);
    }
    return 0;
}
