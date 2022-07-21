#include <stdlib.h>
#include "poisson.h"
#include "setup.h"

#include <math.h>

#define XY(i,j) (((i) - 0.5) * S->G->dx),(((j) - 0.5) * S->G->dy)

double cosines (double x, double y) {
	return sqr(2 * M_PI) * cos(2 * M_PI * x) * cos(2 * M_PI * y);
}

double constant (double x, double y) {
	x = x;
	y = y;
	return 0;
}

double cubic (double x, double y) {
	return 6 * (x + y - 1);
}

int test (simulation * S, double (*f)(double x, double y)) {
	lattice * G = S->G;
	for (int i = 0; i <= G->imax + 1; i++) {
		for (int j = 0; j <= G->jmax + 1; j++) {
			S->RHS[i][j] = f(XY(i,j));
			S->F->P[i][j] = 0;
		}
	}
	int iterations = solvePoissonEquation(S);
	return (iterations < S->max_iterations);
}

int main () {
	int imax = 50, jmax = 25;
	simulation * S = newSimulation(imax, jmax, 1.0 / imax, 1.0 / jmax, 1000);
	int passed = 1;

	passed &= test(S, constant);
	passed &= test(S, cubic);
	passed &= test(S, cosines);

	clearSimulation(S);
	if (!passed)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
