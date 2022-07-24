#include "poisson.h"

int solvePoissonEquation (simulation * S) {
	lattice * G = S->G;
	double ** P = S->F->P;
	double invXWidthSqrd = sqr(1 / G->dx);
	double invYWidthSqrd = sqr(1 / G->dy);
	double scale = S->omega / (2 * (invXWidthSqrd + invYWidthSqrd) );

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
				P[i][j] = (1 - S->omega) * P[i][j] + scale * temp;
			}
		}

		// Compute the error
		error = 0;
		for (int i = 1; i <= G->imax; i++) {
			for (int j = 1; j <= G->jmax; j++) {
				temp = invXWidthSqrd * (P[i + 1][j] - 2 * P[i][j] + P[i - 1][j])
				     + invYWidthSqrd * (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1])
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
			/* Us the following code, if the poisson solver produces a NaN value.
			 * Since no divisions by time-dependent values take place there, any
			 * NaN values must most likely come from here. */
			/*if (isnan(S->RHS[i][j])) {
				fprintf(stderr, "RHS computation yielded zero on input:\n");
				fprintf(stderr, "i: %i, j: %i\nF[i][j] = %g, F[i - 1][j] = %g\nG[i][j] = %g, G[i][j - 1] = %g\n", i, j, F[i][j], F[i - 1][j], G[i][j], G[i][j - 1]);
				fprintf(stderr, "dt = %g, fx = %g, fy = %g\n", dt, fx, fy);
				exit(-1);
			}*/
		}
	}
}
