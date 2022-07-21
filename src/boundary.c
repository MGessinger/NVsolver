#include "boundary.h"

void setPressureOnBoundary (simulation * S) {
	double ** P = S->F->P;
	lattice * G = S->G;

	for (int i = 1; i <= G->imax; i++) {
		P[i][0] = P[i][1];
		P[i][G->jmax + 1] = P[i][G->jmax];
	}
	for (int j = 1; j < G->jmax; j++) {
		P[0][j] = P[1][j];
		P[G->imax + 1][j] = P[G->imax][j];
	}
}

void setVelocitiesOnBoundary (simulation * S) {
	double ** U = S->F->U;
	double ** V = S->F->V;
	lattice * G = S->G;

	/* NOSLIP: */
	for (int i = 1; i <= G->imax; i++) {
		U[i][0] = -U[i][1];
		V[i][0] = 0;
		U[i][G->jmax + 1] = -U[i][G->jmax];
		V[i][G->jmax] = 0;
	}
	for (int j = 1; j <= G->jmax; j++) {
		U[0][j] = 0;
		U[G->imax][j] = 0;
		V[0][j] = -V[1][j];
		V[G->imax + 1][j] = V[G->imax][j];
	}
}
