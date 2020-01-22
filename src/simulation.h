#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "types.h"

lattice runSimulation (REAL ***U, REAL ***V, REAL ***P,
			char *scene, char *paramFile, char *imageFile, int output);

int     simulateFluid (REAL **U, REAL **V, REAL **P,
		bndCond* bCond, lattice *grid, fluidSim *sim, MPI_Comm Region,
		REAL t_end, const char *problem, int opt);
/* Put all functions together and simulate the fluid */

#endif /* SIMULATION_H_ */
