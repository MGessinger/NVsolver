#include <iostream>

#include <gtest/gtest.h>
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

class Computation : public ::testing::Test
{
	protected:
		lattice grid;
		MPI_Comm Region;
		REAL **U = nullptr, **V = nullptr, **P = nullptr;
		virtual void SetUp ()
		{
			int rank, dims[2];
			Region = createCommGrid(&rank,dims);

			grid.imax = grid.jmax = 30;
			grid.delx = grid.dely = 1.0/30;
			splitRegion(Region,rank,dims,&grid);
		}
		virtual	void TearDown ()
		{
			MPI_Comm_free(&Region);
		}
};

static REAL frhs (REAL x, REAL y)
{
	return 2*sqr(2*M_PI)*cos(2*M_PI*x)*cos(2*M_PI*y);
}

static REAL fex (REAL x, REAL y)
{
	return cos(2*M_PI*x)*cos(2*M_PI*y);
}

TEST_F(Computation, Exchange)
{
	int size, rank, coords[2];
	MPI_Comm_size(Region,&size);
	if (size <= 1)
		return;
	MPI_Comm_rank(Region,&rank);
	MPI_Cart_coords(Region,rank,2,coords);
	REAL **P = create2Dfield(grid.deli+2,grid.delj+2);
	fill2Dfield((REAL)coords[0]+coords[1]*100,P,grid.deli+2,grid.delj+2);

	int max = grid.deli;
	if (grid.delj > max)
		max = grid.delj;
	REAL buf1[max+2], buf2[max+2];
	exchangeMat(P,1,1,buf1,buf2,&grid,Region);

	REAL err;
	for (int j = 1; j <= grid.delj; j++)
	{
		if (coords[0] > 0)
			err += sqr(P[0][j] - (100*coords[1])-(coords[0]-1));
	}
	for (int i = 1; i <= grid.deli; i++)
	{
		if (coords[1] > 0)
			err += sqr(P[i][0]-(100*(coords[1]-1))-coords[0]);
	}
	destroy2Dfield((void**)P,grid.deli+2);

	EXPECT_LE(err,0.1);
}

TEST_F(Computation, Timestep)
{
	REAL init[3] = {0.3,0.4,0.5};
	initUVP(&U,&V,&P,grid.deli,grid.delj,init);
	fluidSim sim;
	sim.Re = 1000;
	sim.tau = 0.5;
	sim.dt = 10;
	
	REAL delt = compDelt(&grid,U,V,&sim);
	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)P,grid.deli+2);
	destroy2Dfield((void**)V,grid.deli+2);

	EXPECT_LE(delt-1.0/12,1e-4);
}

TEST_F(Computation, FG)
{
	REAL init[3] = {0.773,14.567,113.4};
	initUVP(&U,&V,&P,grid.deli,grid.delj,init);
	REAL **F = create2Dfield(grid.deli+1, grid.delj+1);
	REAL **G = create2Dfield(grid.deli+1, grid.delj+1);
	REAL **rhs = create2Dfield(grid.deli,grid.delj);
	char **FLAG = create2DIntegerField(grid.deli,grid.delj);

	fluidSim sim;
	sim.Re = 1000;
	sim.GX = 0;
	sim.GY = 0;
	sim.alpha = 0;
	compFG(U,V,F,G,FLAG,0.5,&grid,&sim);
	compRHS(F,G,rhs,FLAG,&grid,0.5);
	REAL errF = 0, errG = 0, errRHS = 0;
	for (int i = 1; i <= grid.deli; i++)
		for (int j = 1; j <= grid.delj; j++)
		{
			errF += sqr(F[i-1][j] - U[i][j]);
			errG += sqr(G[i][j-1] - V[i][j]);
			errRHS += sqr(rhs[i-1][j-1]);
		}

	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);
	destroy2Dfield((void**)F,grid.deli+1);
	destroy2Dfield((void**)G,grid.deli+1);
	destroy2Dfield((void**)rhs,grid.deli);
	destroy2Dfield((void**)FLAG,grid.deli+2);

	EXPECT_LE(errF,1e-4);
	EXPECT_LE(errG,1e-4);
	EXPECT_LE(errRHS,1e-4);
}

TEST_F(Computation, PoissonSolver)
{
	REAL err = 0;
	P = create2Dfield(grid.deli+2,grid.delj+2);
	REAL **rhs = create2Dfield(grid.deli,grid.delj);
	char **FLAG = create2DIntegerField(grid.deli,grid.delj);
	for (int i = 0; i < grid.deli; i++)
		for (int j = 0; j < grid.delj; j++)
			rhs[i][j] = -frhs((i+grid.il+0.5)*grid.delx,(j+grid.jb+0.5)*grid.dely);

	fluidSim sim;
	sim.eps = 1e-5;
	sim.omega = 1.7;
	sim.itmax = 5000;

	REAL buf1[grid.deli+2], buf2[grid.deli+2];
	printf("Performed %i iterations!\n",solveSORforPoisson(P,rhs,FLAG,buf1,buf2,&sim,&grid,Region));
	for (int i = 1; i < grid.deli+1; i++)
		for (int j = 1; j <= grid.delj+1; j++)
			err += sqr(P[i][j] - fex((i+grid.il-0.5)*grid.delx,(j+grid.jb-0.5)*grid.dely));
	printf("The remaining error in L2 was %g.\n",err/sqr(30));
	destroy2Dfield((void**)rhs,grid.deli);
	destroy2Dfield((void**)P,grid.deli+2);
	destroy2Dfield((void**)FLAG,grid.deli+2);
	EXPECT_LE(err/sqr(30),sim.eps);
}
