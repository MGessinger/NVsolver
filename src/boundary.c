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
	short flag = 0;
	/* First set values on the actual boundary of the region */
	if (grid->edges & LEFT)
		for (int j = 0; j <= grid->delj; j++)
			P[0][j] = P[1][j];
	if (grid->edges & RIGHT)
		for (int j = 0; j <= grid->delj; j++)
			P[grid->deli+1][j] = P[grid->deli][j];
	if (grid->edges & BOTTOM)
		for (int i = 0; i <= grid->deli; i++)
			P[i][0] = P[i][1];
	if (grid->edges & TOP)
		for (int i = 0; i <= grid->deli; i++)
			P[i][grid->delj+1] = P[i][grid->delj];
	REAL dxSqrd = sqr(grid->delx);
	REAL dySqrd = sqr(grid->dely);
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			flag = FLAG[i][j];
			if (flag == C_F)
				continue;
			switch (flag ^ C_B)
			{
				case B_N:
					P[i][j] = P[i][j+1];
					break;
				case B_O:
					P[i][j] = P[i+1][j];
					break;
				case B_S:
					P[i][j] = P[i][j-1];
					break;
				case B_W:
					P[i][j] = P[i-1][j];
					break;
				case (B_N | B_O):
					P[i][j] = (dxSqrd*P[i][j+1] + dySqrd*P[i+1][j])/(dxSqrd + dySqrd);
					break;
				case (B_N | B_W):
					P[i][j] = (dxSqrd*P[i][j+1] + dySqrd*P[i-1][j])/(dxSqrd + dySqrd);
					break;
				case (B_S | B_O):
					P[i][j] = (dxSqrd*P[i][j-1] + dySqrd*P[i+1][j])/(dxSqrd + dySqrd);
					break;
				case (B_S | B_W):
					P[i][j] = (dxSqrd*P[i][j-1] + dySqrd*P[i-1][j])/(dxSqrd + dySqrd);
					break;
				default:
					P[i][j] = 0;
			}
		}
	return;
}

void setBCond (REAL **U, REAL **V, lattice *grid, bndCond *bCond)
{
	/* Boundary conditions on the actual boundary of the region */
	if (grid->edges & BOTTOM)
		for (int i = 0; i <= grid->deli; i++)
		{
			U[i+1][0] = (bCond->wb == NOSLIP) ? -U[i+1][1] : U[i+1][1];
			V[i][1] = (bCond->wb == OUTFLOW) ? V[i][2] : 0;
		}
	if (grid->edges & TOP)
		for (int i = 0; i <= grid->deli; i++)
		{
			U[i+1][grid->delj+1] = (bCond->wt == NOSLIP) ? -U[i+1][grid->delj] : U[i+1][grid->delj];
			V[i][grid->delj+1] = (bCond->wt == OUTFLOW) ? V[i][grid->delj] : 0;
		}
	if (grid->edges & LEFT)
		for (int j = 0; j <= grid->delj; j++)
		{
			U[1][j] = (bCond->wl == OUTFLOW) ? U[2][j] : 0;
			V[0][j+1] = (bCond->wl == NOSLIP) ? -V[1][j+1] : V[1][j+1];
		}
	if (grid->edges & RIGHT)
		for (int j = 0; j <= grid->delj; j++)
		{
			U[grid->deli+1][j] = (bCond->wr == OUTFLOW) ? U[grid->deli][j] : 0;
			V[grid->deli+1][j+1] = (bCond->wr == NOSLIP) ? -V[grid->deli][j+1] : V[grid->deli][j+1];
		}
	short flag;
	/* Boundary condtions on obstacles */
	for (int i = 1; i <= grid->deli; i++)
		for (int j = 1; j <= grid->delj; j++)
		{
			flag = bCond->FLAG[i][j];
			if (flag == C_F)
				continue;
			/* Set boundary conditions */
			switch (flag ^ C_B)
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
	if (problem[0] == '\0')
		return;
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
	if (strcmp(problem,"Step") == 0)
	{
		if (!(grid->edges & LEFT))
			return;
		for (int j = 1; j <= grid->delj; j++)
		{
			if (j + grid->jb < grid->jmax/2)
				continue;
			U[1][j] = 1;
			V[0][j+1] = -V[1][j+1];
		}
		return;
	}
	if (strcmp(problem,"Tunnel") == 0 || strcmp(problem,"Von Karman") == 0)
	{
		if (!(grid->edges & LEFT))
			return;
		for (int j = 1; j <= grid->delj; j++)
		{
			U[1][j] = U[0][j] = 1;
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
		for (int i = 0; i < grid->deli+1; i++)
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
		for (int i = 1; i < grid->delj+1; i++)
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
				FLAG[i+j][i] = C_B;
			}
		}
	}
	int size = grid->deli+2;
	if (grid->delj+2 > size)
		size = grid->delj+2;
	char buf[size];
	exchangeIntMat(FLAG,buf,grid,Region);
	for (int i = 0; i < grid->deli+1; i++)
	{
		for (int j = 0; j < grid->delj+1; j++)
		{
			if (FLAG[i][j] == C_F)
				continue;
			int acti = i + grid->il;
			int actj = j + grid->jb;
			if ((actj+1) != grid->delj && FLAG[i][j+1] == C_F)
				FLAG[i][j] |= B_N;
			else if (actj != 0 && FLAG[i][j-1] == C_F)
				FLAG[i][j] |= B_S;
			if ((acti+1) != grid->delj && FLAG[i+1][j] == C_F)
				FLAG[i][j] |= B_O;
			else if (acti != 0 && FLAG[i-1][j] == C_F)
				FLAG[i][j] |= B_W;
		}
	}
	return;
}
