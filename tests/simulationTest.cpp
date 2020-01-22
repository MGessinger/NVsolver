#include <iostream>

#include "gtest/gtest.h"
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

TEST(Setup, Boundary)
{
	lattice grid;
	grid.deli = grid.delj = grid.imax = grid.jmax = 10;
	grid.il = grid.jb = 0;
	grid.edges = LEFT | RIGHT | TOP | BOTTOM;
	grid.delx = grid.dely = 0.1;
	char **FLAG = create2DIntegerField(grid.deli,grid.delj);
	if (FLAG == NULL)
		return;

	initFlags("Step",FLAG,&grid,MPI_COMM_WORLD);
	for (int j = 0; j < grid.jmax+2; j++)
	{
		for (int i = 0; i < grid.imax+2; i++)
			printf("%x,",FLAG[i][j]);
		printf("\n");
	}
	for (int j = 1; j < 4; j++)
		EXPECT_EQ(FLAG[4][j], C_B | B_O);
	for (int i = 1; i < 4; i++)
		EXPECT_EQ(FLAG[i][4], C_B | B_N);
	EXPECT_EQ(FLAG[4][4], C_B | B_N | B_O);
	destroy2Dfield((void**)FLAG,grid.deli+2);
}

TEST(Simulation, Trivial)
{
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	lattice grid = runSimulation(&U, &V, &P, (char*)"Tunnel",
				     (char*)"/home/matthias/Dokumente/Programming/Simulator/data/tunnel", (char*)"None", SILENT);
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
	MPI_Finalize();
}
