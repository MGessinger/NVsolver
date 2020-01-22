#include "poisson.h"

int solveSORforPoisson(REAL **p, REAL **rhs, char **FLAG,
		fluidSim *sim, lattice *grid, MPI_Comm Region)
{
	/* Use a SOR algorithm to solve the poisson equation */
	int  i, j, it = 0;
	REAL temporary, error, eps = sqr(sim->eps);
	REAL invXWidthSqrd = 1/sqr(grid->delx);
	REAL invYWidthSqrd = 1/sqr(grid->dely);
	REAL scale = sim->omega/(2*(invXWidthSqrd+invYWidthSqrd));
	REAL buf[grid->deli + grid->delj + 2];
	/* Count the number of fluid cells */
	int numberOfCells = 0;
	for (i = 1; i <= grid->deli; i++)
		for (j = 1; j <= grid->delj; j++)
		{
			if (FLAG[i][j] == C_F)
				numberOfCells++;
		}
	temporary = eps*numberOfCells;
	MPI_Allreduce(&temporary,&eps,1,MPI_DOUBLE,MPI_SUM,Region);
	do {
		/* Apply the boundary condition */
		applyPboundaryCond(p,grid,FLAG);
		/* Compute the new coefficients iteratively */
		for (i = 1; i <= grid->deli; i++)
			for (j = 1; j <= grid->delj; j++)
			{
				if (FLAG[i][j] != C_F)
					continue;
				temporary = invXWidthSqrd*(p[i+1][j] + p[i-1][j])
					+ invYWidthSqrd*(p[i][j+1] + p[i][j-1])
					- rhs[i-1][j-1];
				p[i][j] = scale*temporary + (1-sim->omega)*p[i][j];
			}
		exchangeMat(p,1,1,buf,grid,Region);
		error = 0;
		/* Calculate the residue with respect to the rhs field */
		for (i = 1; i <= grid->deli; i++)
			for (j = 1; j <= grid->delj; j++)
			{
				if (FLAG[i][j] != C_F)
					continue;
				temporary = invXWidthSqrd*(p[i+1][j] - 2*p[i][j] + p[i-1][j])
					+ invYWidthSqrd*(p[i][j+1] - 2*p[i][j] + p[i][j-1])
					- rhs[i-1][j-1];
				error += sqr(temporary);
			}
		if (++it == sim->itmax)
			break;
		temporary = error;
		MPI_Allreduce(&temporary,&error,1,MPI_DOUBLE,MPI_SUM,Region);
	} while (error > eps);
	return it;
}

REAL compDelt(lattice *grid, REAL **U, REAL **V, fluidSim *sim)
{
	/* Find the optimal step width in time */
	if (sim->tau <= 0)
		return sim->dt;
	REAL dt = sim->Re*(sqr(grid->delx) + sqr(grid->dely))/2;
	REAL utime, vtime;
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			utime = grid->delx/fabs(U[i+1][j]);
			if (utime < dt)
				dt = utime;
			vtime = (grid->dely)/fabs(V[i][j+1]);
			if (vtime < dt)
				dt = vtime;
		}
	dt = sim->tau*dt;
	if (dt > sim->dt)
		return sim->dt;
	return dt;
}

void compRHS(REAL **F, REAL **G, REAL **RHS, char **FLAG, lattice *grid, REAL delt)
{
	REAL facX = delt*(grid->delx);
	REAL facY = delt*(grid->dely);
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			if (FLAG[i][j] != C_F)
				continue;
			RHS[i-1][j-1] = (F[i][j] - F[i-1][j])/facX;
			RHS[i-1][j-1] += (G[i][j] - G[i][j-1])/facY;
		}
	return;
}

void adaptUV(REAL **U, REAL **V, REAL **P, REAL **F, REAL **G,
		REAL delt, char **FLAG, lattice *grid)
{
	REAL facX = delt/(grid->delx);
	REAL facY = delt/(grid->dely);
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			if (FLAG[i][j] != C_F)
				continue;
			U[i+1][j] = F[i][j] - facX*(P[i+1][j] - P[i][j]);
			V[i][j+1] = G[i][j] - facY*(P[i][j+1] - P[i][j]);
		}
	return;
}

REAL delUVbyDelZ(REAL **U, REAL **V, int i, int j, int z, REAL alpha, REAL delz)
{
	REAL duvdz = (U[i+1][j] + U[i+1][j+1])*(V[i][j+1] + V[i+1][j+1]);
	REAL correctionTerm = 0;
	delz *= 4;
	if (z == DERIVE_BY_X)
	{
		duvdz -= (U[i][j] + U[i][j+1])*(V[i-1][j+1] + V[i][j+1]);
		if (alpha == 0)
			return duvdz/delz;
		correctionTerm = fabs(U[i+1][j] + U[i+1][j+1]) * (V[i][j+1] - V[i+1][j+1]);
		correctionTerm -= fabs(U[i][j]+U[i][j+1]) * (V[i-1][j+1] - V[i][j+1]);
	}
	else
	{
		duvdz -= (U[i+1][j-1] + U[i+1][j])*(V[i][j] + V[i+1][j]);
		if (alpha == 0)
			return duvdz/delz;
		correctionTerm = fabs(V[i][j+1] + V[i+1][j+1]) * (U[i+1][j] - U[i+1][j+1]);
		correctionTerm -= fabs(V[i][j] + V[i+1][j]) * (U[i+1][j-1] - U[i+1][j]);
	}
	return (duvdz + alpha*correctionTerm)/delz;
}

REAL delFSqrdByDelZ(REAL **F, int i, int j, int z, REAL alpha, REAL delz)
{
	int dx = (z == DERIVE_BY_X) ? 1 : 0;
	int dy = (z == DERIVE_BY_Y) ? 1 : 0;
	delz *= 4;
	REAL df2dz = sqr(F[i][j] + F[i+dx][j+dy]) - sqr(F[i-dx][j-dy] + F[i][j]);
	if (alpha == 0)
		return df2dz/delz;
	REAL correctionTerm = sqr(F[i][j]) - sqr(F[i+dx][j+dy]);
	if (F[i+dx][j+dy] < -F[i][j])
		correctionTerm *= -1;
	df2dz += alpha*correctionTerm;
	correctionTerm = sqr(F[i-dx][j-dy]) - sqr(F[i][j]);
	if (F[i-dx][j-dy] < - F[i][j])
		correctionTerm *= -1;
	df2dz -= alpha*correctionTerm;
	return df2dz/delz;
}

void    compFG (REAL **U, REAL **V, REAL **F, REAL **G, char **FLAG, REAL delt,
		lattice *grid, fluidSim *simulation)
{
	/* Compute the auxiliary arrays F and G */
	REAL d2ux, d2uy, d2vx, d2vy;
	REAL du2x, dv2y;
	REAL duvx, duvy;
	short flag;
	int i,j;
	for (i = 1; i <= grid->deli; i++)
	{
		for (j = 1; j <= grid->delj; j++)
		{
			flag = FLAG[i][j];
			if (flag == C_B)      /* Boundary cells with no neighboring fluid cells */
				continue;
			else if (flag == C_F) /* Pure fluid cells */
			{
				d2ux = (U[i+2][j] - 2*U[i+1][j] + U[i][j])/sqr(grid->delx);
				d2uy = (U[i+1][j+1] - 2*U[i+1][j] + U[i+1][j-1])/sqr(grid->dely);

				duvy = delUVbyDelZ(U,V,i,j,DERIVE_BY_Y,simulation->alpha,grid->dely);
				du2x = delFSqrdByDelZ(U,i+1,j,DERIVE_BY_X,simulation->alpha,grid->delx);

				F[i][j] = U[i+1][j] + delt*((d2ux+d2uy)/simulation->Re - du2x - duvy + simulation->GX);

				d2vx = (V[i+1][j+1] - 2*V[i][j+1] + V[i-1][j+1])/sqr(grid->delx);
				d2vy = (V[i][j+2] - 2*V[i][j+1] + V[i][j])/sqr(grid->dely);

				duvx = delUVbyDelZ(U,V,i,j,DERIVE_BY_X,simulation->alpha,grid->delx);
				dv2y = delFSqrdByDelZ(V,i,j+1,DERIVE_BY_Y,simulation->alpha,grid->dely);

				G[i][j] = V[i][j+1] + delt*((d2vx+d2vy)/simulation->Re - dv2y - duvx + simulation->GY);
				continue;
			}
			if (flag & B_N)       /* North */
				G[i][j] = V[i][j+1];
			else if (flag & B_S)  /* South */
				G[i-1][j] = V[i-1][j+1];
			if (flag & B_O)       /* East */
				F[i][j] = U[i+1][j];
			else if (flag & B_W)  /* West */
				F[i-1][j] = U[i][j];
		}
	}
	if (grid->edges & LEFT)
	{
		for (j = 1; j<= grid->delj; j++)
			F[0][j] = U[1][j];
	}
	if (grid->edges & RIGHT)
	{
		for (j = 1; j<= grid->delj; j++)
			F[grid->deli][j] = U[grid->deli+1][j];
	}
	if (grid->edges & BOTTOM)
	{
		for (i = 1; i <= grid->deli; i++)
			G[i][0] = V[i][1];
	}
	if (grid->edges & TOP)
	{
		for (i = 1; i <= grid->deli; i++)
			G[i][grid->delj] = V[i][grid->delj+1];
	}
	return;
}
