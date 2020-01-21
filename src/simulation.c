#include "simulation.h"

lattice runSimulation (REAL ***U, REAL ***V, REAL ***P, char *scene, char *paramFile, char *imageFile, int output)
{
	/* Set-Up all memory required to simulate `scene` and then call the simulator */
	lattice grid;
	grid.deli = grid.delj = 0;
	if (!scene || !paramFile || !imageFile)
		return grid;
	if (!U || !V || !P)
		return grid;
	/* Local variables */
	boundaryCond bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
	fluidSim sim;
	REAL init[3];
	REAL delt, t_end;
	int rank, dims[2];

	/* Read parameters from `paramFile` */
	MPI_Comm Region = createCommGrid(&rank,dims);
	if (readParameters(paramFile,init,&grid,&sim,&bCond,&delt,&t_end) < 17)
	{
		printf("The parameter file appears to be incomplete.\n");
		MPI_Abort(Region,0);
	}

	/* Slice Data to process */
	splitRegion(Region, dims, &grid);

	/* Data initialization */
	initUVP(U,V,P,grid.deli,grid.delj,init);
	if (bCond.FLAG == NULL)
		bCond.FLAG = create2DIntegerField(grid.deli+2,grid.delj+2);
	initFlags(scene,bCond.FLAG,&grid,Region);

	/* Run the simulation! */
	int files = simulateFluid(*U,*V,*P,&bCond,&grid,&sim,Region,t_end,scene,output);
	if (output > OUTPUT)
		translateBinary(Region,&grid,files,rank,dims);
	MPI_Comm_free(&Region);
	destroy2Dfield((void**)bCond.FLAG,grid.deli+2);
	return grid;
}

int simulateFluid (REAL **U, REAL **V, REAL **P,
		boundaryCond* bCond, lattice *grid, fluidSim *sim, MPI_Comm Region,
		REAL t_end, const char *problem, int opt)
{
	/* Simulate the fluid characterized by sim, grid and bCond */
	/* Error checking */
	if (!U || !V|| !P)
		return 0;
	if (!bCond || !grid || !sim)
		return 0;
	/* Auxiliary Grids:
	   RHS is used for the Poisson-Solver so no ghost cells are neccessary */
	REAL **F = create2Dfield(grid->deli+1,grid->delj+1);
	REAL **G = create2Dfield(grid->deli+1,grid->delj+1);
	REAL **RHS = create2Dfield(grid->deli,grid->delj);
	if (!F || !G || !RHS)
		return 0;
	int n = 1, rank;
	MPI_Comm_rank(Region,&rank);
	REAL del_vec, dt = 0, delt = sim->dt;
	REAL buf[grid->deli+grid->delj];
	if (opt >= OUTPUT)
		del_vec = t_end/(opt/OUTPUT);
	else
		del_vec = t_end*2;
	/* Begin the simulation */
	if (rank == 0)
		printf("Computing Reynoldsnumber %lg.\n",sim->Re);
	for (REAL time = 0; time <= t_end; time += delt)
	{
		if ((rank == 0) && (opt & PRINT))
			printf("Time is at %lg seconds\n",time);
		/* Update all parameters and fields for the iteration */
		setBCond(U,V,grid,bCond);
		setSpecBCond(U,V,grid,problem);
		compFG(U,V,F,G,bCond->FLAG,delt,grid,sim);
		compRHS(F,G,RHS,bCond->FLAG,grid,delt);

		/* Solve the Poisson Equation */
		solveSORforPoisson(P,RHS,bCond->FLAG,sim,grid, Region);
		/* Update U and V through F,G and P */
		adaptUV(U,V,P,F,G,delt,bCond->FLAG,grid);
		exchangeMat(U,2,1,buf,grid,Region);
		exchangeMat(V,1,2,buf,grid,Region);
		dt = compDelt(grid,U,V,sim);

		if (time > del_vec*n)
		{
			dumpFields(Region,U,V,P,grid,n++);
		}
		MPI_Allreduce(&dt,&delt,1,MPI_DOUBLE,MPI_MIN,Region);
	}
	if (opt >= OUTPUT)
		dumpFields(Region,U,V,P,grid,n++);
	if (rank == 0)
		printf("[Simulation complete!]\n");

	/* Destroy non-simulated grids */
	destroy2Dfield((void**)F,grid->deli+1);
	destroy2Dfield((void**)G,grid->deli+1);
	destroy2Dfield((void**)RHS,grid->deli);
	return n;
}
