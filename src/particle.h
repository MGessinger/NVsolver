#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "real.h"
#include "fields.h"
#include "boundary.h"

typedef struct particle{
    REAL x;
    REAL y;
    REAL u;
    REAL v;
    int onScreen : 2;
} particle;

particle* createParticleArray(int partcount);

void    destroyParticleArray (particle *parts);

int     ParticleSeed (particle *parts, REAL posx1, REAL posx2, REAL posy1, REAL posy2, int partcount, int anzahl);

void    ParticleVelocity (REAL **U, REAL **V, particle *parts, lattice *grid, int partcount);

void    ParticleTransport (particle *parts, int partcount, REAL delt);

#endif /* PARTICLE_H_ */
