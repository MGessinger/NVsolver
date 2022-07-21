#include <stdlib.h>
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

void readArgv (simulation ** S, double * t, int argc, char ** argv) {
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			fprintf(stderr, "Unrecognized parameter: [%s]\n", argv[i]);
			continue;
		}

		switch (argv[i][1]) {
			case 'f':
				*S = newSimulationFromFile(argv[++i]);
				break;
			case 't':
				*t = atof(argv[++i]);
				break;
			default:
				fprintf(stderr, "Unrecognized parameter: [%s]\n", argv[i]);
		}
	}
}

int main (int argc, char ** argv) {
	simulation * S = NULL;
	double dt, t = 1;

	readArgv(&S, &t, argc, argv);

	if (S == NULL) {
		fprintf(stderr, "No simulation found!\n");
		return -1;
	}

	/* Main Loop */
	do {
		dt = computeDT(S);		/* Above */
		setVelocitiesOnBoundary(S);	/* In boundary.h */
		setDrivenCavity(S);
		computeAuxiliaryFields(S, dt);	/* In velocities.h */

		computePoissonRHS(S, dt);	/* In poisson.h */
		int it = solvePoissonEquation(S);	/* In poisson.h */
		printf("Iterations: %i\n", it);

		updateVelocities(S, dt);	/* In velocities.h */

		t -= dt;
	} while (t > 0);
	outputToFile(S);			/* In IO.h */

	clearSimulation(S);
	return 0;
}
