#ifndef POISSON_H_
#define POISSON_H_

#include <math.h>
#include "types.h"

#define DERIVE_BY_X (1)
#define DERIVE_BY_Y (2)

/*REAL**  create2DpoissonMatrix (REAL ilength, REAL jlength, int imax, int jmax);
Compute the 2D-Laplacian in discrete form */

void    applyPboundaryCond(REAL **P, lattice *grid, short **FLAG);
int     solveSOR(REAL **A, REAL *x, REAL *b, int rows, int cols, REAL omega, REAL epsilon, int itermax);
int     solveSORforPoisson (REAL **p, REAL **rhs, short **FLAG, fluidSim *sim, lattice *grid);
/* Solve Ax = b for x. solveSORforPoisson uses an optimised algorithm for this application. */

void    compDelt (REAL *delt, lattice *grid, REAL **U, REAL **V, fluidSim *sim);
/* Copute the next iteration of delt */

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, short **FLAG, REAL delt, lattice *grid, fluidSim *simulation);
void    compRHS (REAL **F, REAL **G, REAL **RHS, short **FLAG, lattice *grid, REAL delt);
/* Compute the RHS of the poisson equation */

void    adaptUV (REAL **U, REAL **V, REAL **P, REAL **F, REAL **G, REAL delt, short **FLAG, lattice *grid);
/* Update the values of U and V */

REAL    delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz);
REAL    delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz);
/* Compute derivatives needed in compFG */

int simulateFluid(REAL **U, REAL **V, REAL **P, boundaryCond* bCond, lattice *grid, fluidSim *sim,
                        REAL delt, REAL t_end, const char *problem, int opt);
/* Put all functions together and simulate the fluid */

#endif /* POISSON_H_ */
