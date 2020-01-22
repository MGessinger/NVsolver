#include <iostream>

#include <gtest/gtest.h>
#include <mpi/mpi.h>

extern "C" {
	#include "types.h"
}

static inline REAL sqr(REAL x)
{
	return x*x;
}

static REAL frhs (REAL x, REAL y)
{
	return 2*sqr(2*M_PI)*cos(2*M_PI*x)*cos(2*M_PI*y);
}

static REAL fex (REAL x, REAL y)
{
	return cos(2*M_PI*x)*cos(2*M_PI*y);
}

int sz = 30;

TEST(Computation, PoissonSolver)
{
	MPI_Init(nullptr,nullptr);
	lattice grid;
	grid.deli = grid.imax = sz;
	grid.delj = grid.jmax = sz;
	grid.delx = grid.dely = 1.0/sz;
	grid.il = grid.jb = sz;
	grid.edges = LEFT | RIGHT | TOP | BOTTOM;

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

	printf("Performed %i iterations!\n",solveSORforPoisson(P,rhs,FLAG,&sim,&grid,MPI_COMM_WORLD));
	//print2Dfield(P,grid.deli,grid.delj);
	for (int i = 1; i < grid.deli+1; i++)
		for (int j = 1; j <= grid.delj+1; j++)
			err += sqr(P[i][j] - fex((i-0.5)*grid.delx,(j-0.5)*grid.dely));
	EXPECT_LE(err/sqr(sz),sim.eps);
	printf("The remaining error in L2 was %g.\n",err/sqr(sz));
	destroy2Dfield((void**)rhs,grid.deli);
	destroy2Dfield((void**)P,grid.deli+2);
	destroy2Dfield((void**)FLAG,grid.deli+2);
}
