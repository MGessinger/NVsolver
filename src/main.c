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

	lattice * G = S->G;
	double ** U = S->F->U;
	double ** V = S->F->V;

	double invDX = 1 / G->dx;
	double invDY = 1 / G->dy;

	double invDT = 2 / (S->F->Reynolds * (sqr(G->dx) + sqr(G->dy)));

	double utime, vtime;
	for (int i = 1; i <= G->imax; i++) {
		for (int j = 1; j <= G->jmax; j++)
		{
			utime = invDX * absd(U[i+1][j]);
			invDT = mind(invDT, utime);

			vtime = invDX * absd(V[i][j+1]);
			invDT = mind(invDT, vtime);
		}
	}
	if (invDT == 0)
		return S->dt;
	return mind(S->tau / invDT, S->dt);
}

void setDrivenCavity (simulation * S) {
	double ** U = S->F->U;

	for (int i = 1; i <= S->G->imax; i++)
		U[i][S->G->jmax + 1] = 2 - U[i][S->G->jmax];
}

void setTunnel (simulation * S) {
	double ** U = S->F->U;
	double ** V = S->F->V;

	lattice * G = S->G;

	for (int j = 0; j < G->jmax + 2; j++) {
		U[0][j] = 1;
		V[0][j] =  - V[1][j];
	}
}

void printUsageInfo () {
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "\t-f: Specifies the file to read simulation parameters from\n");
	fprintf(stderr, "\t-t: Specifies the total time to simulate\n");
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

	if (argc == 1) {
		fprintf(stderr, "No parameters specified! Aborting...\n\n");
		printUsageInfo();
		return -1;
	}

	readArgv(&S, &t, argc, argv);

	if (S == NULL) {
		fprintf(stderr, "No simulation found!\n");
		return -1;
	}

	/* Main Loop */
	do {
		dt = computeDT(S);		/* Above */
		setVelocitiesOnBoundary(S);	/* In boundary.h */
		//setTunnel(S);
		setDrivenCavity(S);
		computeAuxiliaryFields(S, dt);	/* In velocities.h */

		computePoissonRHS(S, dt);	/* In poisson.h */
		int it = solvePoissonEquation(S);	/* In poisson.h */
		if (it <= 0) {
			printf("Aborted with %g left on the clock.\n", t);
			break;
		}
		printf("Iterations: %i\n", it);

		updateVelocities(S, dt);	/* In velocities.h */

		t -= dt;
	} while (t > 0);
	outputToFile(S);			/* In IO.h */

	clearSimulation(S);
	return 0;
}
