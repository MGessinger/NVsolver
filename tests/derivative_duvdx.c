#include <stdlib.h>
#include "types.h"
#include "setup.h"
#include "velocities.h"

#include <math.h>

#define XY_U(i,j) (((i) - 1) * S->G->dx),(((j) - 0.5) * S->G->dy)
#define XY_V(i,j) (((i) - 0.5) * S->G->dx),(((j) - 1) * S->G->dy)

double fu (double x, double y) {
	return cos(x) * y;
}

double fv (double x, double y) {
	return x * cos(y);
}

double df (double x, double y) {
	const double e = 1e-5;
	return (fu(x+e, y) * fv(x+e, y) - fu(x-e, y)*fv(x-e,y)) / (2 * e);
}

int main () {
	int imax = 100, jmax = 100;
	simulation * S = newSimulation(imax, jmax, 1.0 / imax, 1.0 / jmax, 1000);

	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {
			S->F->U[i][j] = fu(XY_U(i,j));
			S->F->V[i][j] = fv(XY_V(i,j));
		}
	}

	double error = 0;
	double duvdx;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			duvdx = df(XY_V(i, j));
			error += absd(duvdx - delUVByDelX(S, i, j));
		}
	}
	error /= (imax * jmax);

	clearSimulation(S);
	if (error > 0.1)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
