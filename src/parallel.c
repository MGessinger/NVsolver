#include "types.h"

MPI_Comm createCommGrid (int *rank, int *dims)
{
	MPI_Comm Region;
	int  nproc, periods[2] = {0,0};
	dims[0] = dims[1] = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Dims_create(nproc,2,dims);
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Region);
	MPI_Comm_rank(Region,rank);
	if (*rank == 0)
		printf("size: %i, dims: [%ix%i]\n",nproc,dims[0],dims[1]);
	return Region;
}

void splitRegion (MPI_Comm Region, int rank, int *dims, lattice *grid)
{
	if (!dims || !grid)
		return;
	grid->edges = 0;
	int coords[2] = {0,0};
	MPI_Cart_coords(Region,rank,2,coords);

	/* Left and right */
	grid->il = coords[0]*grid->imax/dims[0];
	if (grid->il == 0)
		grid->edges |= LEFT;
	int ir = (coords[0]+1)*grid->imax/dims[0];
	if (ir == grid->imax)
		grid->edges |= RIGHT;
	grid->deli = ir - grid->il;

	/* Top and bottom: */
	grid->jb = coords[1]*grid->jmax/dims[1];
	if (grid->jb == 0)
		grid->edges |= BOTTOM;
	int jt = (coords[1]+1)*grid->jmax/dims[1];
	if (jt == grid->jmax)
		grid->edges |= TOP;
	grid->delj = jt - grid->jb;
	return;
}

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf1, REAL *buf2, lattice *grid, MPI_Comm Region)
{
	if (Region == MPI_COMM_WORLD)
		return;
	int prev, next, rank;
	int coords[2];
	MPI_Status st;
	MPI_Comm_rank(Region,&rank);
	MPI_Cart_coords(Region,rank,2,coords);
	int lx = grid->deli+offx+1;
	int ly = grid->delj+offy+1;
	/* Left and right */
	MPI_Cart_shift(Region,0,1,&prev,&next);
	if (coords[0]%2 == 1)
	{
		if (prev >= 0)
			MPI_Sendrecv(mat[offx],ly,MPI_DOUBLE,prev,101,mat[0],ly,MPI_DOUBLE,prev,102,Region,&st);
		if (next >= 0)
			MPI_Sendrecv(mat[lx-offx-1],ly,MPI_DOUBLE,next,104,mat[lx-1],ly,MPI_DOUBLE,next,103,Region,&st);
	}
	else 
	{
		if (next >= 0)
			MPI_Sendrecv(mat[lx-offx-1],ly,MPI_DOUBLE,next,102,mat[lx-1],ly,MPI_DOUBLE,next,101,Region,&st);
		if (prev >= 0)
			MPI_Sendrecv(mat[offx],ly,MPI_DOUBLE,prev,103,mat[0],ly,MPI_DOUBLE,prev,104,Region,&st);
	}
	/* Up and Down */
	MPI_Cart_shift(Region,1,1,&prev,&next);
	if (coords[1]%2 == 1)
	{
		if (prev >= 0)
		{
			for (int i = 0; i < lx; i++)
				buf1[i] = mat[i][offx];
			MPI_Sendrecv(buf1,lx,MPI_DOUBLE,prev,105,buf2,lx,MPI_DOUBLE,prev,106,Region,&st);
			for (int i = 0; i < lx; i++)
				mat[i][0] = buf2[i];
		}
		if (next >= 0)
		{
			for (int i = 0; i < lx; i++)
				buf1[i] = mat[i][ly-offx-1];
			MPI_Sendrecv(buf1,lx,MPI_DOUBLE,next,108,buf2,lx,MPI_DOUBLE,next,107,Region,&st);
			for (int i = 0; i < lx; i++)
				mat[i][ly-1] = buf2[i];
		}
	}
	else
	{
		if (next >= 0)
		{
			for (int i = 0; i < lx; i++)
				buf1[i] = mat[i][ly-offx-1];
			MPI_Sendrecv(buf1,lx,MPI_DOUBLE,next,106,buf2,lx,MPI_DOUBLE,next,105,Region,&st);
			for (int i = 0; i < lx; i++)
				mat[i][ly-1] = buf2[i];
		}
		if (prev >= 0)
		{
			for (int i = 0; i < lx; i++)
				buf1[i] = mat[i][offx];
			MPI_Sendrecv(buf1,lx,MPI_DOUBLE,prev,107,buf2,lx,MPI_DOUBLE,prev,108,Region,&st);
			for (int i = 0; i < lx; i++)
				mat[i][0] = buf2[i];
		}
	}
}

void exchangeIntMat (char **mat, char *bufIn, char *bufOut, lattice *grid, MPI_Comm Region)
{
	if (Region == MPI_COMM_WORLD)
		return;
	int prev, next;
	MPI_Status st;
	int lx = grid->deli+2;
	int ly = grid->delj+2;
	/* Left and right */
	MPI_Cart_shift(Region,0,1,&prev,&next);
	if (prev >= 0)
		MPI_Sendrecv(mat[1],ly,MPI_CHAR,prev,101,mat[0],ly,MPI_CHAR,prev,102,Region,&st);
	if (next >= 0)
		MPI_Sendrecv(mat[grid->deli],ly,MPI_CHAR,next,102,mat[grid->deli+1],ly,MPI_CHAR,next,101,Region,&st);
	/* Up and Down */
	MPI_Cart_shift(Region,1,1,&prev,&next);
	if (prev >= 0)
	{
		for (int i = 0; i < lx; i++)
			bufIn[i] = mat[i][1];
		MPI_Sendrecv(bufIn,lx,MPI_CHAR,prev,103,bufOut,lx,MPI_CHAR,prev,104,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][0] = bufOut[i];
	}
	if (next >= 0)
	{
		for (int i = 0; i < lx; i++)
			bufIn[i] = mat[i][grid->delj];
		MPI_Sendrecv(bufIn,lx,MPI_CHAR,next,104,bufOut,lx,MPI_CHAR,next,103,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][grid->delj+1] = bufOut[i];
	}
	return;
}
