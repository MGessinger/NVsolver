#include "types.h"

MPI_Comm createCommGrid (int *rank, int *dims);

void splitRegion (MPI_Comm Region, int rank, int *dims, lattice *grid);

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf1, REAL *buf2, lattice *grid, MPI_Comm Region);
void exchangeIntMat (char **mat, char *bufIn, char *bufOut, lattice *grid, MPI_Comm Region);
