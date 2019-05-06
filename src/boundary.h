#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <stdlib.h>
#include "real.h"

void    applyHomogeneousNeumannBC (REAL **p, int imax, int jmax);
void    applyHomogeneousDirichletBC (REAL **p, int imax, int jmax);
/* Apply homogeneous boundary conditions to the lattice */

void    setBCond (REAL **U, REAL **V, int imax, int jmax);

#endif /* BOUNDARY_H_ */
