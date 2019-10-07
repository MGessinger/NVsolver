#include "types.h"

MPI_Comm createCommGrid(int *rank, int *dims);

void splitRegion(MPI_Comm Region, int *dims, lattice *grid, char *problem, boundaryCond *bCond);

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf, lattice *grid, MPI_Comm Region);
