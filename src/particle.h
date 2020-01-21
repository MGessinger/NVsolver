#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "types.h"

particle* createParticleArray(int partcount);
void      destroyParticleArray (particle *parts);
/* Create/Delete an array of particles */

void      ParticleSeed(particle *parts, REAL posx1, REAL posx2, REAL posy1, REAL posy2, int partcount, int anzahl);
void      ParticleVelocity (REAL **U, REAL **V, particle *parts, lattice *grid, short **FLAG, int partcount);
void      ParticleTransport (particle *parts, int partcount, REAL delt);
/* Move Particles around on the lattice */

#endif /* PARTICLE_H_ */
