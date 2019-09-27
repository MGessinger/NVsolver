#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "types.h"

#define C_F (0)
#define C_B (1)
#define B_N (2)
#define B_O (4)
#define B_S (8)
#define B_W (16)

#define NOSLIP (0)
#define FREESLIP (1)
#define OUTFLOW (2)

boundaryCond* createBoundCond(int wl, int wr, int wt, int wb);
void    destroyBoundCond(boundaryCond *bCond, int imax);
/* Create and destroy a boundaryCond structure */

void    applyHomogeneousNeumannBC (REAL **p, int imax, int jmax);
void    applyHomogeneousDirichletBC (REAL **p, int imax, int jmax);
/* Apply homogeneous boundary conditions to the lattice */

void    setBCond (REAL **U, REAL **V, lattice *grid, boundaryCond *bCond);
void    setSpecBCond (REAL **U, REAL **V, lattice *grid, const char *problem);
/* Set the boundary conditions for the given problem */

void initFlags(const char* problem, short **FLAG, int imax, int jmax);
/* Set flags for fluid and obstacle cells */

#endif /* BOUNDARY_H_ */
