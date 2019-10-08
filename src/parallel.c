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

void splitRegion(MPI_Comm Region, int *dims, lattice *grid, char *problem, boundaryCond *bCond)
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
    if (bCond->FLAG == NULL)
        bCond->FLAG = create2DIntegerField(grid->deli,grid->delj);
    initFlags(problem,bCond->FLAG,grid);
}

void exchangeMat (REAL **mat, int offx, int offy, REAL *buf, lattice *grid, MPI_Comm Region)
{
    int prev, next, size;
    MPI_Status st;
    /* Top and bottom */
    MPI_Cart_shift(Region,0,1,&prev,&next);
    MPI_Comm_size(Region,&size);
    int lx = grid->deli+2*offx - (offx==2);
    int ly = grid->delj+2*offy - (offy==2);
    if (prev != -1)
        MPI_Send(mat[offx],ly,MPI_DOUBLE,prev,101,Region);
    if (next != -1)
    {
        MPI_Recv(mat[lx-1],ly,MPI_DOUBLE,next,101,Region,&st);
        MPI_Send(mat[lx-offx-1],ly,MPI_DOUBLE,next,102,Region);
    }
    if (prev != -1)
        MPI_Recv(mat[0],ly,MPI_DOUBLE,prev,102,Region,&st);
    /* Left and right */
    MPI_Cart_shift(Region,1,1,&prev,&next);
    if (prev != -1)
    {
        for (int i = 0; i < lx; i++)
            buf[i] = mat[i][0];
        MPI_Send(buf,lx,MPI_DOUBLE,prev,103,Region);
    }
    if (next != -1)
    {
        MPI_Recv(buf,lx,MPI_DOUBLE,next,103,Region,&st);
        for (int i = 0; i < lx; i++)
        {
            mat[i][0] = buf[i];
            buf[i] = mat[i][ly-offx-1];
        }
        MPI_Send(buf,lx,MPI_DOUBLE,next,104,Region);
    }
    if (prev != -1)
    {
        MPI_Recv(buf,lx,MPI_DOUBLE,prev,104,Region,&st);
        for (int i = 0; i < lx; i++)
            mat[i][ly-1] = buf[i];
    }
    return;
}
