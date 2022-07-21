#ifndef SETUP_H_
#define SETUP_H_

#include "types.h"

simulation * newSimulation (int imax, int jmax, double dx, double dy, double Re);
simulation * newSimulationFromFile (char * filename);
void clearSimulation (simulation * S);

#endif /* SETUP_H_ */
