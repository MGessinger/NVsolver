#include <stdlib.h>
#include <stdio.h>

#include "fields.h"
#include "real.h"
#include "poisson.h"
#include "IO.h"
#include "boundary.h"
#include "particle.h"

int main (int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Usage: .\\simulator parameter_file [number_of_frames]\n");
        return 0;
    }
    REAL **U = NULL, **V = NULL, **P = NULL;
    int out = (argc >= 3) ? atoi(argv[2]) : 0;

    lattice *grid = simulateFluid(&U,&V,&P,argv[1],PRINT | out*OUTPUT);
    outputVec(U,V,P,NULL,grid,0,0);
    /* Destroy simulated grids */
    if (grid != NULL)
    {
        destroyMatrix(U,grid->imax+2);
        destroyMatrix(V,grid->imax+2);
        destroyMatrix(P,grid->imax+2);
        free(grid);
    }
    return 0;
}
