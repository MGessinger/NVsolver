#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "setup.h"

/* ==================== Creation ==================== */

double ** newDoubleField (int imax, int jmax) {
	double ** X = malloc(imax * sizeof(double*));
	if (X == NULL) {
		fprintf(stderr, "Could not allocate array. Please try again.\n");
		return NULL;
	}

	for (int i = 0; i < imax; i++) {
		X[i] = malloc(jmax * sizeof(double));
		if (X[i] == NULL) {
			fprintf(stderr, "Could not allocate array. Please try again.\n");
			clearDoubleField(X, i);
			return NULL;
		}
		memset(X[i], 0, jmax * sizeof(double));
	}

	return X;
}

lattice * newLattice (double dx, double dy, int imax, int jmax) {
	lattice * G = malloc(sizeof(lattice));
	if (G == NULL) {
		fprintf(stderr, "Could not allocate lattice. Please try again.\n");
		return NULL;
	}

	G->dx = dx;
	G->dy = dy;
	G->imax = imax;
	G->jmax = jmax;

	G->bc_top = STICK;
	G->bc_bottom = STICK;
	G->bc_left = STICK;
	G->bc_right = STICK;

	return G;
}

fluid * newFluid (double Re, int imax, int jmax) {
	fluid * F = malloc(sizeof(fluid));
	if (F == NULL) {
		fprintf(stderr, "Could not allocate fluid. Please try again.\n");
		return NULL;
	}

	F->Reynolds = Re;
	F->U = newDoubleField(imax + 2, jmax + 2);
	F->V = newDoubleField(imax + 2, jmax + 2);
	F->P = newDoubleField(imax + 2, jmax + 2);

	return F;
}

simulation * newSimulationFromFile (char * filename) {
	FILE * parameter_file = fopen(filename, "r");
	if (parameter_file == NULL)
		return NULL;

	char * line = NULL;
	size_t len = 0;

	simulation * S = NULL;
	int imax = 1, jmax = 1;
	double dx = 1, dy = 1;
	double Re = 100;

	for (int k = 0; k < 5; k++) {
		getline(&line, &len, parameter_file);
		if (strncmp("imax", line, 4) == 0)
			imax = atoi(line + 5);
		else if (strncmp("jmax", line, 4) == 0)
			jmax = atoi(line + 5);
		else if (strncmp("Reynolds", line, 8) == 0)
			Re = atof(line + 9);
		else if (strncmp("dx", line, 2) == 0)
			dx = atof(line + 3);
		else if (strncmp("dy", line, 2) == 0)
			dy = atof(line + 3);
		else {
			fprintf(stderr, "This file is formatted incorrectly!\nPlease provide imax, jmax, dx, dy and Re first.\n");
			return NULL;
		}
	}

	S = newSimulation(imax, jmax, dx, dy, Re);
	
	while (getline(&line, &len, parameter_file) > 0) {
		if (strncmp("tau", line, 3) == 0)
			S->tau = atof(line + 4);
		else if (strncmp("GX", line, 2) == 0)
			S->GX = atof(line + 3);
		else if (strncmp("GY", line, 2) == 0)
			S->GY = atof(line + 3);
		else if (strncmp("dt", line, 2) == 0)
			S->dt = atof(line + 3);
		else if (strncmp("alpha", line, 5) == 0)
			S->alpha = atof(line + 6);
		else if (strncmp("omega", line, 5) == 0)
			S->omega = atof(line + 6);
		else if (strncmp("itmax", line, 5) == 0)
			S->max_iterations = atoi(line + 6);
		else if (strncmp("epsilon", line, 7) == 0)
			S->solver_tol = atof(line + 8);
		else
			fprintf(stderr, "Unrecognized option: [%s]\n", line);
	}

	fclose(parameter_file);
	if (line != NULL)
		free(line);
	return S;
}

simulation * newSimulation (int imax, int jmax, double dx, double dy, double Re) {
	simulation * S = malloc(sizeof(simulation));
	if (S == NULL) {
		fprintf(stderr, "Could not allocate simulation. Please try again.\n");
		return NULL;
	}

	S->G = newLattice(dx, dy, imax, jmax);
	S->F = newFluid(Re, imax, jmax);

	S->auxF = newDoubleField(imax + 2, jmax + 2);
	S->auxG = newDoubleField(imax + 2, jmax + 2);
	S->RHS = newDoubleField(imax + 2, jmax + 2);

	S->GX = S->GY = 0;
	S->tau = 0.5;
	S->dt = 0.02;
	S->alpha = 0.5;
	S->omega = 1.7;
	S->max_iterations = 5000;
	S->solver_tol = 0.01;

	return S;
}

/* ==================== Deletion ==================== */

void clearDoubleField (double ** X, int imax) {
	for (int i = 0; i < imax; i++)
		free(X[i]);
	free(X);
}

void clearLattice (lattice * G) {
	free(G);
}

void clearFluid (fluid * F, int imax) {
	clearDoubleField(F->U, imax);
	clearDoubleField(F->V, imax);
	clearDoubleField(F->P, imax);
	free(F);
}

void clearSimulation (simulation * S) {
	clearFluid(S->F, S->G->imax + 2);
	clearDoubleField(S->auxF, S->G->imax + 2);
	clearDoubleField(S->auxG, S->G->imax + 2);
	clearDoubleField(S->RHS, S->G->imax + 2);
	clearLattice(S->G);
	free(S);
}
