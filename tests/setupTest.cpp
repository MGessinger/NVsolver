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
		MPI_Comm Region;
		virtual void SetUp ()
		{
			int rank, dims[2];
			Region = createCommGrid(&rank,dims);

			grid.imax = grid.jmax = 10;
			grid.delx = 0.2;
			grid.dely = 0.1;
			splitRegion(Region,rank,dims,&grid);

			bCond = createBoundCond(NOSLIP,OUTFLOW,FREESLIP,NOSLIP);
			bCond.FLAG = create2DIntegerField(grid.deli,grid.delj);
		}
		virtual	void TearDown ()
		{
			destroy2Dfield((void**)bCond.FLAG,grid.deli+2);
			MPI_Comm_free(&Region);
		}
};

int cmpFlags (char **FLAG, int ib, int jb)
{
	if (ib < 0 || jb < 0)
		return 1;
	for (int j = 1; j <= jb; j++)
		if (FLAG[ib][j] & B_O == 0)
			return 0;
	for (int i = 1; i <= ib; i++)
		if (FLAG[i][jb] & B_N == 0)
			return 0;
	return 1;
}

TEST(MPI, Init)
{
	EXPECT_EQ(MPI_Init(nullptr, nullptr),MPI_SUCCESS);
}

TEST_F(Setup, Filling)
{
	REAL **U = create2Dfield(grid.deli,grid.delj);
	fill2Dfield(1.25,U,grid.deli,grid.delj);
	REAL diff = 0;
	for (int i = 0; i < grid.deli; i++)
		for (int j = 0; j < grid.delj; j++)
			diff += sqr(U[i][j]-1.25);
	destroy2Dfield((void**)U,grid.deli);
	EXPECT_LE(diff,1e-10);
}

TEST_F(Setup, Boundary)
{
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	REAL init[3] = {1.0,1.0,1.0};
	initFlags("Step",bCond.FLAG,&grid,Region);
	int ib = grid.imax/2-grid.il, jb = grid.jmax/2-grid.jb;
	if (ib > grid.deli)
		ib = grid.deli;
	if (jb > grid.delj)
		jb = grid.delj;
	ASSERT_TRUE(cmpFlags(bCond.FLAG,ib,jb));
	initUVP(&U,&V,&P,grid.imax,grid.jmax,init);

	setBCond(U,V,&grid,&bCond);
	setSpecBCond(U,V,&grid,"Step");
	applyPbndCond(P,&grid,bCond.FLAG);
/*
	print2Dfield(U,grid.deli+3,grid.delj+2);
	print2Dfield(V,grid.deli+2,grid.delj+3);
	print2Dfield(P,grid.deli+2,grid.delj+2);
*/
	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);
}

#define OUT_FILE "/home/matthias/Dokumente/Programming/Simulator/data/MomentumField.vtk"
TEST_F(Setup, IO)
{
	int dims[2], periods[2], coords[2];
	MPI_Cart_get(Region,2,dims,periods,coords);
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
	REAL init[3] = {1.7,0.773,2.1};
	initUVP(&U,&V,&P,grid.imax,grid.jmax,init);

	dumpFields(Region,U,V,P,&grid,0);
	MPI_Barrier(Region);
	translateBinary(Region,&grid,1,0,dims);

	int size;
	MPI_Comm_size(Region,&size);
	MPI_Barrier(Region);
	if (size == 1)
		writeVTKfileFor2DvectorField(OUT_FILE,"momentumfield",U,V,&grid);

	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);

	if (size == 0)
		EXPECT_EQ(system("diff -q MomentumField_0.vtk " OUT_FILE),0);
}
