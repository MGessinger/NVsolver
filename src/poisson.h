#ifndef POISSON_H_
#define POISSON_H_

#include <math.h>
#include "types.h"

#define DERIVE_BY_X (1)
#define DERIVE_BY_Y (2)

void    applyPboundaryCond(REAL **P, lattice *grid, char **FLAG);
int     solveSORforPoisson (REAL **p, REAL **rhs, char **FLAG, fluidSim *sim, lattice *grid, MPI_Comm Region);
/* Solve Ax = b for x. solveSORforPoisson uses an optimised algorithm for this application. */

REAL    compDelt(lattice *grid, REAL **U, REAL **V, fluidSim *sim);
/* Copute the next iteration of delt */

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, char **FLAG, REAL delt, lattice *grid, fluidSim *simulation);
void    compRHS (REAL **F, REAL **G, REAL **RHS, char **FLAG, lattice *grid, REAL delt);
/* Compute the RHS of the poisson equation */

void    adaptUV (REAL **U, REAL **V, REAL **P, REAL **F, REAL **G, REAL delt, char **FLAG, lattice *grid);
/* Update the values of U and V */

REAL    delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz);
REAL    delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz);
/* Compute derivatives needed in compFG */

#endif /* POISSON_H_ */
