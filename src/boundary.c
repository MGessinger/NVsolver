#include "boundary.h"

void setPressureOnBoundary (simulation * S) {
	double ** P = S->F->P;
	lattice * G = S->G;

	int i_ghost = G->imax + 1;
	int j_ghost = G->jmax + 1;

	for (int i = 1; i <= G->imax; i++) {
		P[i][0] = P[i][1];
		P[i][j_ghost] = P[i][j_ghost - 1];
	}
	for (int j = 1; j <= G->jmax; j++) {
		P[0][j] = P[1][j];
		P[i_ghost][j] = P[i_ghost - 1][j];
	}
}

void setVelocitiesOnBoundary (simulation * S) {
	double ** U = S->F->U;
	double ** V = S->F->V;
	lattice * G = S->G;

	int i, j;
	int i_ghost = G->imax + 1;
	int j_ghost = G->jmax + 1;

	/* TOP (j maximal) */
	switch (G->bc_top) {
		case STICK:
			for (i = 0; i <= i_ghost; i++) {
				U[i][j_ghost] = -1 * U[i][j_ghost - 1];
				V[i][j_ghost - 1] = 0;
			}
			break;
		case SLIP:
			for (i = 0; i <= i_ghost; i++) {
				U[i][j_ghost] = U[i][j_ghost - 1];
				V[i][j_ghost - 1] = 0;
			}
			break;
		case OUTFLOW:
			for (i = 0; i <= i_ghost; i++) {
				U[i][j_ghost] = U[i][j_ghost - 1];
				V[i][j_ghost - 1] = V[i][j_ghost - 2];
			}
			break;
		default: break;
	}

	/* BOTTOM (j minimal) */
	switch (G->bc_bottom) {
		case STICK:
			for (i = 0; i <= i_ghost; i++) {
				U[i][0] = -1 * U[i][1];
				V[i][0] = 0;
			}
			break;
		case SLIP:
			for (i = 0; i <= i_ghost; i++) {
				U[i][0] = U[i][1];
				V[i][0] = 0;
			}
			break;
		case OUTFLOW:
			for (i = 0; i <= i_ghost; i++) {
				U[i][0] = U[i][1];
				V[i][0] = V[i][1];
			}
			break;
		default: break;
	}

	/* LEFT (i minimal) */
	switch (G->bc_left) {
		case STICK:
			for (j = 0; j <= j_ghost; j++) {
				U[0][j] = 0;
				V[0][j] = -1 * V[1][j];
			}
			break;
		case SLIP:
			for (j = 0; j <= j_ghost; j++) {
				U[0][j] = 0;
				V[0][j] = V[1][j];
			}
			break;
		case OUTFLOW:
			for (j = 0; j <= j_ghost; j++) {
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
			break;
		default: break;
	}

	/* RIGHT (i maximal) */
	switch (G->bc_left) {
		case STICK:
			for (j = 0; j <= j_ghost; j++) {
				U[i_ghost - 1][j] = 0;
				V[i_ghost][j] = -1 * V[i_ghost - 1][j];
			}
			break;
		case SLIP:
			for (j = 0; j <= j_ghost; j++) {
				U[i_ghost - 1][j] = 0;
				V[i_ghost][j] = V[i_ghost - 1][j];
			}
			break;
		case OUTFLOW:
			for (j = 0; j <= j_ghost; j++) {
				U[i_ghost - 1][j] = U[i_ghost - 2][j];
				V[i_ghost][j] = V[i_ghost - 1][j];
			}
			break;
		default: break;
	}
}
