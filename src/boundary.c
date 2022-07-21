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

	int i, j;
	int i_ghost = G->imax + 1;
	int j_ghost = G->jmax + 1;

	/* TOP (j maximal) */
	switch (G->bc_top) {
		case STICK:
			for (i = 1; i <= G->imax; i++) {
				U[i][j_ghost] = -U[i][G->jmax];
				V[i][G->jmax] = 0;
			}
			break;
		case SLIP:
			for (i = 1; i <= G->imax; i++) {
				U[i][j_ghost] = U[i][G->jmax];
				V[i][G->jmax] = 0;
			}
			break;
		case OUTFLOW:
			for (i = 1; i <= G->imax; i++) {
				U[i][j_ghost] = U[i][G->jmax];
				V[i][G->jmax] = V[i][G->jmax - 1];
			}
			break;
		default: break;
	}

	/* BOTTOM (j minimal) */
	switch (G->bc_bottom) {
		case STICK:
			for (i = 1; i <= G->imax; i++) {
				U[i][0] = -U[i][1];
				V[i][0] = 0;
			}
			break;
		case SLIP:
			for (i = 1; i <= G->imax; i++) {
				U[i][0] = U[i][1];
				V[i][0] = 0;
			}
			break;
		case OUTFLOW:
			for (i = 1; i <= G->imax; i++) {
				U[i][0] = U[i][1];
				V[i][0] = V[i][1];
			}
			break;
		default: break;
	}

	/* LEFT (i minimal) */
	switch (G->bc_left) {
		case STICK:
			for (j = 1; j <= G->jmax; j++) {
				U[0][j] = 0;
				V[0][j] = -V[1][j];
			}
			break;
		case SLIP:
			for (j = 1; j <= G->jmax; j++) {
				U[0][j] = 0;
				V[0][j] = V[1][j];
			}
			break;
		case OUTFLOW:
			for (j = 1; j <= G->jmax; j++) {
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
			break;
		default: break;
	}

	/* RIGHT (i maximal) */
	switch (G->bc_left) {
		case STICK:
			for (j = 1; j <= G->jmax; j++) {
				U[G->imax][j] = 0;
				V[i_ghost][j] = -V[G->imax][j];
			}
			break;
		case SLIP:
			for (j = 1; j <= G->jmax; j++) {
				U[G->imax][j] = 0;
				V[i_ghost][j] = V[G->imax][j];
			}
			break;
		case OUTFLOW:
			for (j = 1; j <= G->jmax; j++) {
				U[G->imax][j] = U[G->imax - 1][j];
				V[i_ghost][j] = V[G->imax][j];
			}
			break;
		default: break;
	}
}
