#include <stdio.h>
#include "types.h"
#include "setup.h"
#include "poisson.h"
#include "boundary.h"
#include "velocities.h"
#include "IO.h"

double computeDT (simulation * S) {
	/* Find the optimal step width in time */
	double dt = S->F->Reynolds * (sqr(S->G->dx) + sqr(S->G->dy)) / 2;

	lattice * G = S->G;
	double ** U = S->F->U;
	double ** V = S->F->V;

	double utime, vtime;
	for (int i = 1; i <= G->imax; i++) {
		for (int j = 1; j <= G->jmax; j++)
		{
			utime = G->dx / absd(U[i+1][j]);
			dt = maxd(dt, utime);

			vtime = G->dy / absd(V[i][j+1]);
			dt = maxd(dt, vtime);
		}
	}
	return mind(S->tau * dt, S->dt);
}

void setDrivenCavity (simulation * S) {
	double ** U = S->F->U;

	for (int i = 1; i <= S->G->imax; i++)
		U[i][S->G->jmax + 1] = 2 - U[i][S->G->jmax];
}

void setOpenTunnel (simulation * S) {
	double ** U = S->F->U;
	double ** V = S->F->V;

	for (int j = 0; j <= S->G->jmax + 1; j++) {
		U[0][j] = U[1][j];
		U[S->G->imax + 1][j] = U[S->G->imax][j];

		V[0][j] = 0;
		V[S->G->imax + 1][j] = 0;
	}
}

int main () {
	int imax = 32, jmax = 64;
	simulation * S = newSimulation(imax, jmax, 1.0 / imax, 1.0 / jmax, 500);
	S->GX = 1;

	double dt = S->dt;
	double t = 1;
	int it;
	/* Main Loop */
	do {
		dt = computeDT(S);		/* Above */
		setVelocitiesOnBoundary(S);	/* In boundary.h */
		setOpenTunnel(S);
		computeAuxiliaryFields(S, dt);	/* In velocities.h */

		computePoissonRHS(S, dt);	/* In poisson.h */
		it = solvePoissonEquation(S);	/* In poisson.h */
		printf("Iterations: %i\n", it);

		updateVelocities(S, dt);	/* In velocities.h */

		t -= dt;
	} while (t > 0);
	outputToFile(S);			/* In IO.h */

	clearSimulation(S);
	return 0;
}
