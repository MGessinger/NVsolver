#include "types.h"

bndCond createBoundCond (int wl, int wr, int wt, int wb)
{
	bndCond bCond;
	bCond.wl = wl;
	bCond.wr = wr;
	bCond.wt = wt;
	bCond.wb = wb;
	bCond.FLAG = NULL;
	return bCond;
}

void applyPbndCond (REAL **P, lattice *grid, char **FLAG)
{
	/* Apply boundary conditions for the pressure field */
	char flag = 0;
	/* First set values on the actual boundary of the region */
	if (grid->edges & LEFT)
		for (int j = 1; j <= grid->delj; j++)
			P[0][j] = P[1][j];
	if (grid->edges & RIGHT)
		for (int j = 1; j <= grid->delj; j++)
			P[grid->deli+1][j] = P[grid->deli][j];
	if (grid->edges & BOTTOM)
		for (int i = 1; i <= grid->deli; i++)
			P[i][0] = P[i][1];
	if (grid->edges & TOP)
		for (int i = 1; i <= grid->deli; i++)
			P[i][grid->delj+1] = P[i][grid->delj];
	REAL dxSqrd = sqr(grid->delx);
	REAL dySqrd = sqr(grid->dely);
	REAL num, denom;
	int delNS, delOW;
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			flag = LO_NIBBLE(FLAG[i][j]);
			if (flag == C_F)
				continue;
			num = 0;
			denom = (flag == C_B) + ((flag & B_NS) != 0)*dxSqrd + ((flag & B_OW) != 0)*dySqrd;
			delNS = ((flag & B_N) != 0) - ((flag & B_S) != 0);
			delOW = ((flag & B_O) != 0) - ((flag & B_W) != 0);

			if (flag & B_NS)
				num += P[i][j+delNS] * dxSqrd;
			if (flag & B_OW)
				num += P[i+delOW][j] * dySqrd;

			P[i][j] = num / denom;
		}
	return;
}

void setBCond (REAL **U, REAL **V, lattice *grid, bndCond *bCond)
{
	/* Boundary conditions on the actual boundary of the region */
	int fac1, fac2;
	if (grid->edges & BOTTOM)
	{
		fac1 = 1-2*(bCond->wb == NOSLIP);
		fac2 = (bCond->wb == OUTFLOW);
		for (int i = 1; i <= grid->deli; i++)
		{
			U[i+1][0] = fac1 * U[i+1][1];
			V[i][1] = fac2 * V[i][2];
		}
	}
	if (grid->edges & TOP)
	{
		fac1 = 1-2*(bCond->wt == NOSLIP);
		fac2 = (bCond->wt == OUTFLOW);
		for (int i = 1; i <= grid->deli; i++)
		{
			U[i+1][grid->delj+1] = fac1 * U[i+1][grid->delj];
			V[i][grid->delj+1] = fac2 * V[i][grid->delj];
		}
	}
	if (grid->edges & LEFT)
	{
		fac1 = (bCond->wl == OUTFLOW);
		fac2 = 1-2*(bCond->wl == NOSLIP);
		for (int j = 1; j <= grid->delj; j++)
		{
			U[1][j] = fac1 * U[2][j];
			V[0][j+1] = fac2 * V[1][j+1];
		}
	}
	if (grid->edges & RIGHT)
	{
		fac1 = (bCond->wr == OUTFLOW);
		fac2 = 1-2*(bCond->wr == NOSLIP);
		for (int j = 1; j <= grid->delj; j++)
		{
			U[grid->deli+1][j] = fac1 * U[grid->deli][j];
			V[grid->deli+1][j+1] = fac2 * V[grid->deli][j+1];
		}
	}
	char flag;
	/* Boundary condtions on obstacles */
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			flag = LO_NIBBLE(bCond->FLAG[i][j]);
			if (flag == C_F)
				continue;
			/* Set boundary conditions */
			switch (flag)
			{
				case B_N:
					V[i][j+1] = 0;
					U[i+1][j] = -U[i+1][j+1];
					U[i][j] = -U[i][j+1];
					break;
				case B_O:
					U[i+1][j] = 0;
					V[i][j+1] = -V[i+1][j+1];
					V[i][j] = -V[i+1][j];
					break;
				case B_S:
					V[i][j] = 0;
					U[i+1][j] = -U[i+1][j-1];
					U[i][j] = -U[i][j-1];
					break;
				case B_W:
					U[i][j] = 0;
					V[i][j+1] = -V[i-1][j+1];
					V[i][j] = -V[i-1][j];
					break;
				case B_N | B_O:
					U[i+1][j] = V[i][j+1] = 0;
					U[i][j] = -U[i][j+1];
					V[i][j] = -V[i+1][j];
					break;
				case B_N | B_W:
					U[i][j] = V[i][j+1] = 0;
					U[i+1][j] = -U[i+1][j+1];
					V[i][j] = -V[i-1][j];
					break;
				case B_S | B_O:
					U[i+1][j] = V[i][j] = 0;
					U[i][j] = -U[i][j-1];
					V[i][j+1] = -V[i+1][j+1];
					break;
				case B_S | B_W:
					U[i][j] = V[i][j] = 0;
					U[i+1][j] = -U[i+1][j-1];
					V[i][j+1] = -V[i-1][j+1];
					break;
				default:
					U[i+1][j] = V[i][j+1] = 0;
					break;
			}
		}
	return;
}

void setSpecBCond (REAL **U, REAL **V, lattice *grid, const char *problem)
{
	/* Set special (e.g. inflow) conditions */
	if (strcmp(problem,"Driven Cavity") == 0)
	{
		if (!(grid->edges & TOP))
			return;
		for (int i = 2; i <= grid->deli+1; i++)
		{
			U[i][grid->delj+1] = 2-U[i][grid->delj];
			V[i-1][grid->delj+2] = 0;
		}
		return;
	}
	if (strcmp(problem,"Tunnel") == 0)
	{
		if (!(grid->edges & LEFT))
			return;
		for (int j = 1; j <= grid->delj; j++)
		{
			U[1][j]  = 1;
			V[0][j+1] = -V[1][j+1];
		}
		return;
	}
	return;
}

void initFlags (const char *problem, char **FLAG, lattice *grid, MPI_Comm Region)
{
	/* Manually sets the flag field for arbitrary generalised geometries.
	 * If the flags are read from a file, set problem to "Image"! */
	if (strcmp(problem,"Step") == 0)
	{
		int max = grid->jmax;
		if (grid->imax < max)
			max = grid->imax;
		max = (max+2)/2;
		for (int i = 0; i < grid->deli+2; i++)
		{
			if (i + grid->il >= max)
				break;
			for (int j = 0; j < grid->delj+2; j++)
			{
				if (j + grid->jb >= max)
					break;
				FLAG[i][j] = C_B;
			}
		}
	}
	else if (strcmp(problem,"Von Karman") == 0)
	{
		for (int i = 0; i < grid->delj+1; i++)
		{
			if (i + grid->il < grid->jmax/3)
				continue;
			if (i + grid->il >= 2*grid->jmax/3)
				break;
			for (int j = -(grid->jmax/4)/2; j <= (grid->jmax/4)/2; j++)
			{
				int x = i+j+grid->il+grid->jb;
				if (x < grid->jmax/3)
					continue;
				if (x >= 2*grid->jmax/3 || x >= grid->deli)
					break;
				FLAG[i+j+1][i+1] = C_B;
			}
		}
	}
	int size = grid->deli+2;
	if (grid->delj+2 > size)
		size = grid->delj+2;
	char buf1[size], buf2[size];
	exchangeIntMat(FLAG,buf1,buf2,grid,Region);
	for (int i = 0; i < grid->deli+1; i++)
	{
		for (int j = 0; j < grid->delj+1; j++)
		{
			if (FLAG[i][j] == C_F)
				continue;
			if ((j != grid->delj) && (FLAG[i][j+1] == C_F))
				FLAG[i][j] |= B_N;
			else if ((j != 0) && (FLAG[i][j-1] == C_F))
				FLAG[i][j] |= B_S;
			if ((i != grid->delj) && (FLAG[i+1][j] == C_F))
				FLAG[i][j] |= B_O;
			else if ((i != 0) && (FLAG[i-1][j] == C_F))
				FLAG[i][j] |= B_W;
		}
	}
	return;
}
