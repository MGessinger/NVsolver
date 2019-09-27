#include "types.h"

int main (int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Usage: .\\simulator parameter_file [number_of_frames]\n");
        return 0;
    }
    REAL **U = NULL, **V = NULL, **P = NULL;
    int out = (argc >= 3) ? (PRINT | atoi(argv[2])*OUTPUT) : SILENT;

    lattice *grid = simulateFluid(&U,&V,&P,argv[1],out);
    outputVec(U,V,P,NULL,grid,0,0);
    /* Destroy simulated grids */
    if (grid != NULL)
    {
        destroy2Dfield(U,grid->imax+2);
        destroy2Dfield(V,grid->imax+2);
        destroy2Dfield(P,grid->imax+2);
        free(grid);
    }
    return 0;
}
