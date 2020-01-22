#ifndef TYPES_H_
#define TYPES_H_

#include <stdlib.h>
#include <string.h>
#include <mpi/mpi.h>

#define REAL double

#define TOP (0x1)
#define BOTTOM (0x2)
#define LEFT (0x4)
#define RIGHT (0x8)

typedef struct lattice {
	REAL delx;
	REAL dely;
	int imax, jmax;     /* Max coords of the entire grid */
	int il, jb;         /* Bottom left of the partial grid */
	int deli, delj;     /* Max coords of the partial grid */
	char edges;
} lattice;

typedef struct fluidSimulation {
	REAL Re;
	REAL tau, dt;
	REAL GX, GY;
	REAL eps;
	REAL omega;
	REAL alpha;
	int itmax;
} fluidSim;

typedef struct bndCondition {
	unsigned wt : 2;
	unsigned wr : 2;
	unsigned wb : 2;
	unsigned wl : 2;
	char **FLAG;
} bndCond;

typedef struct particle {
	REAL x;
	REAL y;
	REAL u;
	REAL v;
	unsigned onScreen : 1;
} particle;

static inline REAL sqr (REAL x)
{
	return x*x;
}

#include "boundary.h"
#include "fields.h"
#include "IO.h"
#include "particle.h"
#include "poisson.h"
#include "parallel.h"
#include "simulation.h"

#define OUTPUT (0x10)
#define PRINT  (0x01)
#define SILENT (0x0)

#endif /* TYPES_H_ */
