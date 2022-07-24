#include "velocities.h"

double delUVByDelX (simulation * S, int i, int j) {
	double ** U = S->F->U;
	double ** V = S->F->V;

	double duvdx, correction;
	double du1, du2;

	du1 = U[i][j] + U[i][j + 1];
	du2 = U[i - 1][j] + U[i - 1][j + 1];

	duvdx = du1 * (V[i][j] + V[i + 1][j])
	      - du2 * (V[i - 1][j] + V[i][j]);
	correction = absd(du1) * (V[i][j] - V[i + 1][j])
		   - absd(du2) * (V[i - 1][j] - V[i][j]);
	return (duvdx + S->alpha * correction) / (4 * S->G->dx);
}

double delUVByDelY (simulation * S, int i, int j) {
	double ** U = S->F->U;
	double ** V = S->F->V;

	double duvdy, correction;
	double dv1, dv2;

	dv1 = V[i][j] + V[i + 1][j];
	dv2 = V[i][j - 1] + V[i + 1][j - 1];

	duvdy = dv1 * (U[i][j] + U[i][j + 1])
	      - dv2 * (U[i][j - 1] + U[i][j]);
	correction = absd(dv1) * (U[i][j] - U[i][j + 1])
		   - absd(dv2) * (U[i][j - 1] - U[i][j]);

	return (duvdy + S->alpha * correction) / (4 * S->G->dy);
}

double delU2ByDelX (simulation * S, int i, int j) {
	double ** U = S->F->U;

	double du2dx, correction;
	double du1, du2;

	du1 = U[i][j] + U[i + 1][j];
	du2 = U[i - 1][j] + U[i][j];

	du2dx = du1 * du1 - du2 * du2;
	correction = absd(du1) * (U[i][j] - U[i + 1][j])
		   - absd(du2) * (U[i - 1][j] - U[i][j]);

	return (du2dx + S->alpha * correction) / (4 * S->G->dx);
}

double delV2ByDelY (simulation * S, int i, int j) {
	double ** V = S->F->V;

	double dv2dy, correction;
	double dv1, dv2;

	dv1 = V[i][j] + V[i][j + 1];
	dv2 = V[i][j - 1] + V[i][j];

	dv2dy = dv1 * dv1 - dv2 * dv2;
	correction = absd(dv1) * (V[i][j] - V[i][j + 1])
		   - absd(dv2) * (V[i][j - 1] - V[i][j]);

	return (dv2dy + S->alpha * correction) / (4 * S->G->dy);
}

void computeAuxiliaryFields (simulation * S, double dt) {
	/* Variables for the derivatives: */
	double d2udx, d2udy, d2vdx, d2vdy;
	double du2dx, dv2dy;
	double duvdx, duvdy;

	/* Variables for fluid parameters */
	double ** U = S->F->U;
	double ** V = S->F->V;
	lattice * G = S->G;
	double invDxSqrd = sqr(1 / G->dx);
	double invDySqrd = sqr(1 / G->dy);

	for (int i = 1; i <= G->imax; i++) {
		for (int j = 1; j <= G->jmax; j++) {
			d2udx = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) * invDxSqrd;
			d2udy = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) * invDySqrd;
			du2dx = delU2ByDelX(S, i, j);
			duvdx = delUVByDelX(S, i, j);
			duvdy = delUVByDelY(S, i, j);
			dv2dy = delV2ByDelY(S, i, j);
			d2vdx = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) * invDxSqrd;
			d2vdy = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) * invDySqrd;

			S->auxF[i][j] = U[i][j] + dt * ( (d2udx + d2udy) / S->F->Reynolds - du2dx - duvdy + S->GX );
			S->auxG[i][j] = V[i][j] + dt * ( (d2vdx + d2vdy) / S->F->Reynolds - duvdx - dv2dy + S->GY );
		}

		/* Sneak in the boundary term as well */
		S->auxG[i][0] = V[i][0] + dt * S->GY;
		S->auxG[i][G->jmax] = V[i][G->jmax] + dt * S->GY;
	}

	for (int j = 1; j <= G->jmax; j++) {
		S->auxF[0][j] = U[0][j] + dt * S->GX;
		S->auxF[G->imax][j] = U[G->imax][j] + dt * S->GX;
	}
}

void updateVelocities (simulation * S, double dt) {
	double ** U = S->F->U;
	double ** V = S->F->V;
	double ** P = S->F->P;
	lattice * G = S->G;

	double fx = dt / G->dx;
	double fy = dt / G->dy;

	for (int i = 1; i <= G->imax; i++) {
		for (int j = 1; j <= G->jmax; j++) {
			U[i][j] = S->auxF[i][j] - fx * (P[i + 1][j] - P[i][j]);
			V[i][j] = S->auxG[i][j] - fy * (P[i][j + 1] - P[i][j]);
		}
	}
}
