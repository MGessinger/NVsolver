#include <stdio.h>
#include "IO.h"

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
