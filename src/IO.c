#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "IO.h"
#include "setup.h"

/* ==================== INPUT ==================== */

boundary atobc (char * type) {
	if (strncmp("STICK", type, 5) == 0)
		return STICK;
	else if (strncmp("SLIP", type, 4) == 0)
		return SLIP;
	else if (strncmp("OUTFLOW", type, 7) == 0)
		return OUTFLOW;
	return SPECIAL;
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
		else if (strncmp("bc_top", line, 6) == 0)
			S->G->bc_top = atobc(line + 7);
		else if (strncmp("bc_bottom", line, 9) == 0)
			S->G->bc_bottom = atobc(line + 10);
		else if (strncmp("bc_left", line, 7) == 0)
			S->G->bc_left = atobc(line + 8);
		else if (strncmp("bc_right", line, 8) == 0)
			S->G->bc_right = atobc(line + 9);
		else
			fprintf(stderr, "Unrecognized option: [%s]\n", line);
	}

	fclose(parameter_file);
	if (line != NULL)
		free(line);
	return S;
}

/* ==================== OUTPUT ==================== */

void outputToFile (simulation * S) {
	static unsigned n = 0;

	char p_filename[128] = {'\0'};
	char u_filename[128] = {'\0'};
	char v_filename[128] = {'\0'};

	sprintf(p_filename, "../data/pressure%u.dat", n);
	sprintf(u_filename, "../data/horizontal%u.dat", n);
	sprintf(v_filename, "../data/vertical%u.dat", n);

	FILE * p_file = fopen(p_filename, "w");
	FILE * u_file = fopen(u_filename, "w");
	FILE * v_file = fopen(v_filename, "w");

	fluid * F = S->F;

	fprintf(p_file, "[");
	fprintf(u_file, "[");
	fprintf(v_file, "[");
	for (int i = 1; i <= S->G->imax; i++) {
		fprintf(p_file, "[");
		fprintf(u_file, "[");
		fprintf(v_file, "[");
		for (int j = 1; j <= S->G->jmax; j++) {
			fprintf(p_file, "%g, ", F->P[i][j]);
			fprintf(u_file, "%g, ", (F->U[i - 1][j] + F->U[i][j]) / 2);
			fprintf(v_file, "%g, ", (F->V[i][j - 1] + F->V[i][j]) / 2);
		}
		fprintf(p_file, "],\n");
		fprintf(u_file, "],\n");
		fprintf(v_file, "],\n");
	}
	fprintf(p_file, "]\n");
	fprintf(u_file, "]\n");
	fprintf(v_file, "]\n");

	fclose(p_file);
	fclose(u_file);
	fclose(v_file);

	n++;
}
