#include "types.h"

MPI_Comm createCommGrid (int *rank, int *dims);

void splitRegion (MPI_Comm Region, int *dims, lattice *grid);

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf, lattice *grid, MPI_Comm Region);
void exchangeIntMat (char **mat, char *buf, lattice *grid, MPI_Comm Region);
