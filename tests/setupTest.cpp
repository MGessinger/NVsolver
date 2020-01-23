#include <iostream>

#include "gtest/gtest.h"
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

class Setup : public ::testing::Test {
	protected:
		lattice grid;
		bndCond bCond;
		virtual void SetUp ()
		{
			grid.deli = grid.delj = grid.imax = grid.jmax = 10;
			grid.il = grid.jb = 0;
			grid.edges = LEFT | RIGHT | TOP | BOTTOM;
			grid.delx = grid.dely = 0.1;

			bCond = createBoundCond(NOSLIP,OUTFLOW,FREESLIP,NOSLIP);
			bCond.FLAG = create2DIntegerField(grid.deli,grid.delj);
		}
		virtual	void TearDown ()
		{
			destroy2Dfield((void**)bCond.FLAG,grid.deli+2);
		}
};

int cmpFlags (char **FLAG, int ib, int jb)
{
	for (int j = 0; j < jb; j++)
		if (FLAG[ib][j] != (C_B | B_O))
			return 0;
	for (int i = 0; i < ib; i++)
		if (FLAG[i][jb] != (C_B | B_N))
			return 0;
	if (FLAG[ib][jb] != (C_B | B_N | B_O))
		return 0;
	return 1;
}

TEST_F(Setup, Filling)
{
	REAL **U = create2Dfield(grid.deli,grid.delj);
	fill2Dfield(1.25,U,grid.deli,grid.delj);
	REAL diff = 0;
	for (int i = 0; i < grid.deli; i++)
		for (int j = 0; j < grid.delj; j++)
			diff += sqr(U[i][j]-1.25);
	EXPECT_LE(diff,1e-10);
	destroy2Dfield((void**)U,grid.deli);
}

TEST_F(Setup, Flags)
{
	initFlags("Step",bCond.FLAG,&grid,MPI_COMM_WORLD);
	for (int j = 0; j < grid.jmax+2; j++)
	{
		for (int i = 0; i < grid.imax+2; i++)
			printf("%x,",bCond.FLAG[i][j]);
		printf("\n");
	}
	EXPECT_TRUE(cmpFlags(bCond.FLAG,grid.deli/2,grid.deli/2));
}

TEST_F(Setup, Boundary)
{
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	REAL init[3] = {1.0,1.0,1.0};
	initUVP(&U,&V,&P,grid.imax,grid.jmax,init);
	initFlags("Step",bCond.FLAG,&grid,MPI_COMM_WORLD);

	setBCond(U,V,&grid,&bCond);
	setSpecBCond(U,V,&grid,"Step");
	applyPbndCond(P,&grid,bCond.FLAG);
	print2Dfield(U,grid.deli+3,grid.delj+2);
	print2Dfield(V,grid.deli+2,grid.delj+3);
	print2Dfield(P,grid.deli+2,grid.delj+2);

	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);
}
