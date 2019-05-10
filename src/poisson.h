#ifndef POISSON_H_
#define POISSON_H_

#include "fields.h"
#include "blas.h"
#include "boundary.h"
#include "IO.h"

#define DERIVE_BY_X (1)
#define DERIVE_BY_Y (2)

REAL**  create2DpoissonMatrix (REAL ilength, REAL jlength, int imax, int jmax);
/* Compute the 2D-Laplacian in discrete form */

int     solveSOR(REAL **A, REAL *x, REAL *b, int rows, int cols, REAL omega, REAL epsilon, int itermax);
int     solveSORforPoisson (REAL **p, REAL **rhs, REAL omega, REAL epsilon, int itermax, int useNeumann, lattice *grid);
/* Solve Ax = b for x. solveSORforPoisson uses an optimised algorithm for this application. useNeumann
 * toggles whether homogeneous Neumann (1) or Dirichlet (0) boundary conditions are to be used for p. */

void    compDelt (REAL *delt, REAL imax, REAL jmax, REAL delx, REAL dely, REAL **U, REAL **V, REAL Re, REAL tau);
/* Copute the next iteration of delt */

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, int imax, int jmax,
                REAL delt, REAL delx, REAL dely, fluidSim *simulation);
void    compRHS (REAL **F, REAL **G, REAL **RHS, int imax, int jmax, REAL delx, REAL dely, REAL delt);
/* Compute the RHS of the poisson equation */

void    adaptUV (REAL **U, REAL **V, REAL **P, REAL **F, REAL **G, REAL delt, REAL delx, REAL dely, int imax, int jmax);
/* Update the values of U and V */

REAL    delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz);
REAL    delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz);
/* Compute derivatives needed in compFG */

lattice*    simulateFluid (REAL ***U, REAL ***V, REAL ***P, const char *fileName, int option, boundaryCond *bCond);
/* Put all functions together and simulate the fluid */

#endif /* POISSON_H_ */
