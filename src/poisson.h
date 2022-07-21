#ifndef POISSON_H_
#define POISSON_H_

#include "types.h"
#include "boundary.h"

void computePoissonRHS (simulation * S, double dt);
int solvePoissonEquation (simulation * S);

#endif /* POISSON_H_ */
