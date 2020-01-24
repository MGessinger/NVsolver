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

TEST_F(Computation, Timestep)
{
	REAL **U = nullptr, **V = nullptr, **P = nullptr;
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

TEST_F(Computation, PoissonSolver)
{
	REAL err = 0;
	REAL **P = create2Dfield(grid.deli+2,grid.delj+2);
	REAL **rhs = create2Dfield(grid.deli,grid.delj);
	char **FLAG = create2DIntegerField(grid.deli,grid.delj);
	for (int i = 0; i < grid.deli; i++)
		for (int j = 0; j < grid.delj; j++)
			rhs[i][j] = -frhs((i+0.5)*grid.delx,(j+0.5)*grid.dely);

	fluidSim sim;
	sim.eps = 1e-5;
	sim.omega = 1.7;
	sim.itmax = 5000;

	printf("Performed %i iterations!\n",solveSORforPoisson(P,rhs,FLAG,&sim,&grid,Region));
	for (int i = 1; i < grid.deli+1; i++)
		for (int j = 1; j <= grid.delj+1; j++)
			err += sqr(P[i][j] - fex((i-0.5)*grid.delx,(j-0.5)*grid.dely));
	printf("The remaining error in L2 was %g.\n",err/sqr(30));
	destroy2Dfield((void**)rhs,grid.deli);
	destroy2Dfield((void**)P,grid.deli+2);
	destroy2Dfield((void**)FLAG,grid.deli+2);
	EXPECT_LE(err/sqr(30),sim.eps);
}
