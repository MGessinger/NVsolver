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
    int height = 32, width = 32;
    short **FLAG = readGeometry(argv[1],&height,&width);
    initFlags("Image",FLAG,height,width);
    for (int j = height-1; j >= 0; j--)
    {
        for (int i = 0; i < width; i++)
            printf("%2x, ",FLAG[i][j]);
        printf("\n");
    }
    destroy2DIntegerField(FLAG,height);
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
