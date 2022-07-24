#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "setup.h"
#include "velocities.h"

#include <math.h>

#define XY_U(i,j) (((i) - 1) * S->G->dx),(((j) - 0.5) * S->G->dy)
#define XY_V(i,j) (((i) - 0.5) * S->G->dx),(((j) - 1) * S->G->dy)

double fu (double x, double y) {
	return cos(x) * cos(y);
}

double fv (double x, double y) {
	return sqr(sqr(x + y));
}

double ff (double x, double y) {
	double u = fu(x, y);
	double Du = -2 * u;
	double du2dx = 2 * u * cos(y) * -sin(x);
	double duvdy = fv(x, y) * cos(x) * -sin(y) + 4 * sqr(x + y) * (x + y) * u;
	return u + 0.01 * (Du / 1000 - du2dx - duvdy);
}

double fg (double x, double y) {
	double v = fv(x, y);
	double Dv = 12 * 2 * sqr(x + y);
	double duvdx = v * -sin(x) * cos(y) + fu(x, y) * 4 * sqr(x + y) * (x + y);
	double dv2dy = 8 * sqr(sqr(x + y) * (x + y)) * (x + y);
	return v + 0.01 * (Dv / 1000 - duvdx - dv2dy);
}

int main () {
	int imax = 100, jmax = 100;
	double dt = 0.01;
	simulation * S = newSimulation (imax, jmax, 1.0 / imax, 1.0 / jmax, 1000);
	S->GX = -0.375;
	S->GY = -0.75;

	for (int i = 0; i <= imax + 1; i++) {
		for (int j = 0; j <= jmax + 1; j++) {
			S->F->U[i][j] = fu(XY_U(i, j));
			S->F->V[i][j] = fv(XY_V(i, j));
		}
	}

	computeAuxiliaryFields(S, dt);

	double expF, expG;
	double avg_err_f = 0, avg_err_g = 0;
	for (int i = 1; i < S->G->imax; i++) {
		for (int j = 1; j < S->G->jmax; j++) {
			expF = ff(XY_U(i, j)) + dt * S->GX;
			expG = fg(XY_U(i, j)) + dt * S->GY;

			avg_err_f += absd(S->auxF[i][j] - expF);
			avg_err_g += absd(S->auxG[i][j] - expG);
		}
	}

	avg_err_f /= imax * jmax;
	avg_err_g /= imax * jmax;

	clearSimulation(S);
	printf("Average errors: %g, %g\n", avg_err_f, avg_err_g);
	if (avg_err_f > 0.01 || avg_err_g > 0.01)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
