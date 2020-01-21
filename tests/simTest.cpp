#include <iostream>

#include "gtest/gtest.h"
#include <mpi/mpi.h>

extern "C" {
#include "types.h"
}

#define DATA_FILE "/home/matthias/Dokumente/Programming/Simulator/data/tunnel"

TEST(Simulation, Trivial)
{
	MPI_Init(nullptr,nullptr);
	boundaryCond bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
	fluidSim sim;
	lattice grid;
	int rank, dims[2];
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	REAL **U1 = nullptr, **V1 = nullptr, **P1 = nullptr;
	REAL init[3];
	REAL delt, t_end;
	ASSERT_GE(readParameters(DATA_FILE,init,&grid,&sim,&bCond,&delt,&t_end),17);
	initUVP(&U,&V,&P,grid.deli,grid.delj,init);
	initUVP(&U1,&V1,&P1,grid.deli,grid.delj,init);
	if (bCond.FLAG == NULL)
		bCond.FLAG = create2DIntegerField(grid.deli+2,grid.delj+2);
	initFlags("Tunnel",bCond.FLAG,&grid,MPI_COMM_WORLD);
	simulateFluid(U,V,P,&bCond,&grid,&sim,MPI_COMM_WORLD,t_end,"Tunnel",0);
	for (int i = 1; i <= grid.deli; i++)
		for (int j = 1; j <= grid.delj; j++)
		{
			EXPECT_LE(U[i+1][j]-U1[i+1][j],1e-10);
			EXPECT_LE(V[i][j+1]-V1[i][j+1],1e-10);
			EXPECT_LE(P[i][j]-P1[i][j],1e-10);
		}
	destroy2Dfield(U,grid.deli+3);
	destroy2Dfield(V,grid.deli+2);
	destroy2Dfield(P,grid.deli+2);

	destroy2Dfield(U1,grid.deli+3);
	destroy2Dfield(V1,grid.deli+2);
	destroy2Dfield(P1,grid.deli+2);
	destroy2DIntegerField(bCond.FLAG,grid.deli+2);
	MPI_Finalize();
}
