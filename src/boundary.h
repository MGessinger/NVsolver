#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <stdlib.h>
#include "real.h"

typedef struct boundaryCondition {
    int wt : 3;
    int wr : 3;
    int wb : 3;
    int wl : 3;
} boundaryCond;

#define NOSLIP (0)
#define FREESLIP (1)
#define OUTFLOW (2)

void    applyHomogeneousNeumannBC (REAL **p, int imax, int jmax);
void    applyHomogeneousDirichletBC (REAL **p, int imax, int jmax);
/* Apply homogeneous boundary conditions to the lattice */

void    setBCond (REAL **U, REAL **V, int imax, int jmax, boundaryCond *bCond);
void    setSpecBCond (REAL **U, REAL **V, int imax, int jmax, char *problem);
/* Set the boundary conditions for the given problem */

#endif /* BOUNDARY_H_ */
