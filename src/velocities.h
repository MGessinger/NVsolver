#ifndef UV_H_
#define UV_H_

#include "types.h"

double delUVByDelX (simulation * S, int i, int j);
double delUVByDelY (simulation * S, int i, int j);
double delU2ByDelX (simulation * S, int i, int j);
double delV2ByDelY (simulation * S, int i, int j);

void computeAuxiliaryFields (simulation * S, double dt);
void updateVelocities (simulation * S, double dt);

#endif /* UV_H_ */
