#include "parallel.h"

MPI_Comm createCommGrid(int *rank, int *dims)
{
	MPI_Comm Region;
	int  nproc, periods[2] = {0,0};
	dims[0] = dims[1] = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Dims_create(nproc,2,dims);
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Region);
	MPI_Comm_rank(Region,rank);
	if (*rank == 0)
		printf("size: %i, dims: [%i,%i]\n",nproc,dims[0],dims[1]);
	return Region;
}

void splitRegion(MPI_Comm Region, int *dims, lattice *grid)
{
	if (!dims || !grid)
		return;
	grid->edges = 0;
	int rank;
	int coords[2] = {0,0};
	MPI_Comm_rank(Region,&rank);
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

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf, lattice *grid, MPI_Comm Region)
{
	if (Region == MPI_COMM_WORLD)
		return;
	int prev, next, rank;
	int coords[2];
	MPI_Status st;
	/* Top and bottom */
	MPI_Cart_shift(Region,0,1,&prev,&next);
	MPI_Comm_rank(Region,&rank);
	MPI_Cart_coords(Region,rank,2,coords);
	int lx = grid->deli+2*offx - (offx==2);
	int ly = grid->delj+2*offy - (offy==2);
	if (coords[0]%2 == 1 && prev != -1)
		MPI_Sendrecv(mat[offx],ly,MPI_DOUBLE,prev,101,mat[0],ly,MPI_DOUBLE,prev,101,Region,&st);
	else if (coords[0]%2 == 0 && next != -1)
		MPI_Sendrecv(mat[lx-offx-1],ly,MPI_DOUBLE,next,101,mat[lx-offx],ly,MPI_DOUBLE,next,101,Region,&st);
	if (coords[0]%2 == 0 && prev != -1)
		MPI_Sendrecv(mat[offx],ly,MPI_DOUBLE,prev,102,mat[0],ly,MPI_DOUBLE,prev,102,Region,&st);
	else if (coords[0]%2 == 1 && next != -1)
		MPI_Sendrecv(mat[lx-offx-1],ly,MPI_DOUBLE,next,102,mat[lx-offx],ly,MPI_DOUBLE,next,102,Region,&st);
	/* Left and right */
	return;
	MPI_Cart_shift(Region,1,1,&prev,&next);
	if (prev != -1)
	{
		for (int i = 0; i < lx; i++)
			buf[i] = mat[i][0];
		MPI_Sendrecv(buf,lx,MPI_DOUBLE,prev,103,buf,lx,MPI_DOUBLE,prev,104,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][ly-1] = buf[i];
	}
	if (next != -1)
	{
		for (int i = 0; i < lx; i++)
			buf[i] = mat[i][ly-offx];
		MPI_Sendrecv(buf,lx,MPI_DOUBLE,next,104,buf,lx,MPI_DOUBLE,next,103,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][0] = buf[i];
	}
	return;
}

void exchangeIntMat (char **mat, char *buf, lattice *grid, MPI_Comm Region)
{
	if (Region == MPI_COMM_WORLD)
		return;
	int prev, next, size;
	MPI_Status st;
	/* Top and bottom */
	MPI_Cart_shift(Region,0,1,&next,&prev);
	MPI_Comm_size(Region,&size);
	if (size == 0)
		return;
	int lx = grid->deli+2;
	int ly = grid->delj+2;
	if (prev != -1)
		MPI_Sendrecv(mat[1],ly,MPI_CHAR,prev,101,mat[0],ly,MPI_CHAR,prev,102,Region,&st);
	if (next != -1)
		MPI_Sendrecv(mat[grid->deli],ly,MPI_CHAR,next,102,mat[lx-1],ly,MPI_CHAR,next,101,Region,&st);
	/* Left and right */
	MPI_Cart_shift(Region,1,1,&prev,&next);
	if (prev != -1)
	{
		for (int i = 0; i < lx; i++)
			buf[i] = mat[i][0];
		MPI_Sendrecv(buf,lx,MPI_CHAR,prev,103,buf,lx,MPI_CHAR,prev,104,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][ly-1] = buf[i];
	}
	if (next != -1)
	{
		for (int i = 0; i < lx; i++)
			buf[i] = mat[i][grid->delj];
		MPI_Sendrecv(buf,lx,MPI_CHAR,next,104,buf,lx,MPI_CHAR,next,103,Region,&st);
		for (int i = 0; i < lx; i++)
			mat[i][0] = buf[i];
	}
	return;
}
