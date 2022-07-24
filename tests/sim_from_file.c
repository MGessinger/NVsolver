#include <stdlib.h>
#include <stdio.h>
#include "IO.h"
#include "setup.h"

int checkLattice (lattice * G) {
	if (G->imax != 123)
		return 1;
	if (G->jmax != 321)
		return 2;

	if (G->dx != 0.125)
		return 3;
	if (G->dy != 0.75)
		return 4;

	if (G->bc_top != STICK)
		return 5;
	if (G->bc_bottom != SLIP)
		return 6;
	if (G->bc_left != OUTFLOW)
		return 7;
	if (G->bc_right != SPECIAL)
		return 9;

	return 0;
}

int checkSimulation (simulation * S) {
	if (S->GX != -0.375)
		return 1;
	if (S->GY != -0.875)
		return 2;

	if (S->tau != 0.25)
		return 3;
	if (S->dt != 2.5)
		return 4;

	if (S->alpha != 1.125)
		return 5;
	if (S->omega != 0.875)
		return 6;
	
	if (S->solver_tol != 0.0625)
		return 7;
	if (S->max_iterations != 4321)
		return 8;

	return 0;
}

int main () {
	simulation * S = newSimulationFromFile("../../data/test_parameter_file");

	if (S == NULL) {
		fprintf(stderr, "Reading from file failed big time\n");
		return EXIT_FAILURE;
	}

	int ltt = checkLattice(S->G);
	if (ltt != 0) {
		fprintf(stderr, "lattice check: %i\n", ltt);
		return EXIT_FAILURE;
	}

	int Re = S->F->Reynolds;
	if (Re != 54321) {
		fprintf(stderr, "Reynolds number: %i\n", Re);
		return EXIT_FAILURE;
	}

	int sim = checkSimulation(S);
	if (sim != 0) {
		fprintf(stderr, "Simulation check: %i\n", sim);
		return EXIT_FAILURE;
	}

	clearSimulation(S);
	return EXIT_SUCCESS;
}
