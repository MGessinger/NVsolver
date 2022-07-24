#ifndef IO_H_
#define IO_H_

#include "types.h"

simulation * newSimulationFromFile (char * filename);
void outputToFile (simulation * S);

#endif /* IO_H_ */
