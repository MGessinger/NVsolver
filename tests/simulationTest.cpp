#include <iostream>

#include "gtest/gtest.h"
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

TEST(Simulation, Trivial)
{
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	lattice grid = runSimulation(&U, &V, &P, (char*)"Tunnel",
				     (char*)"/home/matthias/Dokumente/Programming/Simulator/data/testfluid", (char*)"None", PRINT);
	REAL err = 0;
	for (int i = 1; i <= grid.deli; i++)
		for (int j = 1; j <= grid.delj; j++)
		{
			err += U[i+1][j]-1;
			err += V[i][j+1];
			err += P[i][j];
		}
	EXPECT_LE(err,1e-10);
	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);
}
