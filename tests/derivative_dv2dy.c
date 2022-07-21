#include <stdlib.h>
#include "types.h"
#include "setup.h"
#include "velocities.h"

#include <math.h>

#define XY_V(i,j) (((i) - 0.5) * S->G->dx),(((j) - 1) * S->G->dy)

double fv (double x, double y) {
	return cos(x) * y * y;
}

double df (double x, double y) {
	const double e = 1e-5;
	return (sqr(fv(x, y+e)) - sqr(fv(x, y-e))) / (2 * e);
}

int main () {
	int imax = 100, jmax = 100;
	simulation * S = newSimulation(imax, jmax, 1.0 / imax, 1.0 / jmax, 1000);

	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {
			S->F->V[i][j] = fv(XY_V(i,j));
		}
	}

	double error = 0;
	double dv2dy;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			dv2dy = df(XY_V(i, j));
			error += absd(dv2dy - delV2ByDelY(S, i, j));
		}
	}
	error /= (imax * jmax);

	clearSimulation(S);
	if (error > 0.1)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
