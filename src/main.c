#include <stdlib.h>
#include <stdio.h>

#include "fields.h"
#include "blas.h"
#include "real.h"
#include "poisson.h"
#include "IO.h"
#include "boundary.h"

#include <math.h>

void writeHorizontalAxis(REAL **U, int imax, int jmax)
{
    FILE *out = open_file("horizontal.dat","w");
    if (out == NULL)
        return;
    for (int j = 1; j <= jmax; j++)
        fprintf(out,"%i,%lg\n",j,U[imax/2-1][j]);
    fclose(out);
    return;
}

int main (int argc, char **argv)
{
    REAL **U, **V, **P;
    int out = (argc >= 3) ? atoi(argv[2]) : 10;
    boundaryCond *bCond = createBoundCond(0,0,NOSLIP,NOSLIP,NOSLIP,NOSLIP);
    lattice *grid = simulateFluid(&U,&V,&P,(argc >= 2 ? argv[1] : "dcavity.par"),PRINT | out*OUTPUT, bCond);
    outputVec(U,V,P,grid,0);
    writeHorizontalAxis(U,grid->imax,grid->jmax);
    /* Destroy simulated grids */
    destroyMatrix(U,grid->imax+2);
    destroyMatrix(V,grid->imax+2);
    destroyMatrix(P,grid->imax+2);
    destroyBoundCond(bCond,grid->imax+2);
    if (grid != NULL)
        free(grid);
    return 0;
}
