#include "poisson.h"

int solvePoissonEquation (simulation * S) {
	lattice * G = S->G;
	double ** P = S->F->P;
	double invXWidthSqrd = sqr(1 / G->dx);
	double invYWidthSqrd = sqr(1 / G->dy);
	double scale = 2 * (invXWidthSqrd + invYWidthSqrd);
	double factor = S->omega / scale;

	double temp = 0;
	double error = 0;
	double max_error = sqr(S->solver_tol) * (G->imax * G->jmax);
	int iterations = 0;
	do {
		// Update the boundary values
		setPressureOnBoundary(S);

		// Compute the new P values
		for (int i = 1; i <= G->imax; i++) {
			for (int j = 1; j <= G->jmax; j++) {
				temp = invXWidthSqrd * (P[i + 1][j] + P[i - 1][j])
				     + invYWidthSqrd * (P[i][j + 1] + P[i][j - 1])
				     - S->RHS[i][j];
				P[i][j] = (1 - S->omega) * P[i][j] + factor * temp;
			}
		}

		// Compute the error
		error = 0;
		for (int i = 1; i <= G->imax; i++) {
			for (int j = 1; j <= G->jmax; j++) {
				temp = invXWidthSqrd * (P[i + 1][j] + P[i - 1][j])
				     + invYWidthSqrd * (P[i][j + 1] + P[i][j - 1])
				     - scale * P[i][j]
				     - S->RHS[i][j];
				error += sqr(temp);
			}
		}

		if (++iterations == S->max_iterations)
			break;
	} while (error > max_error);
	return iterations;
}

void computePoissonRHS (simulation * S, double dt) {
	double ** F = S->auxF;
	double ** G = S->auxG;
	lattice * Grid = S->G;

	double fx = 1 / (dt * Grid->dx);
	double fy = 1 / (dt * Grid->dy);

	for (int i = 1; i <= Grid->imax; i++) {
		for (int j = 1; j <= Grid->jmax; j++) {
			S->RHS[i][j] = (F[i][j] - F[i - 1][j]) * fx
				     + (G[i][j] - G[i][j - 1]) * fy;
		}
	}
}
