#include "parallel.h"

MPI_Comm createCommGrid(int *rank, int *dims)
{
    MPI_Comm Region;
    int  nproc, periods[2] = {0,0};
    dims[0] = dims[1] = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Dims_create(nproc,2,dims);
    printf("size: %i, dims: [%i,%i]\n",nproc,dims[0],dims[1]);
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&Region);
    MPI_Comm_rank(Region,rank);
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
    if (!buf || !grid)
        return;
    int prev, next, size;
    MPI_Status st;
    /* Top and bottom */
    MPI_Cart_shift(Region,1,1,&prev,&next);
    MPI_Comm_size(Region,&size);
    if (prev != -1)
        MPI_Send(&(mat[offx][1]),grid->delj,MPI_DOUBLE,prev,101,Region);
    if (next != -1)
    {
        MPI_Recv(&(mat[grid->deli+2*offx-1][1]),grid->delj,MPI_DOUBLE,next,101,Region,&st);
        MPI_Send(&(mat[grid->deli+offx][1]),grid->delj,MPI_DOUBLE,next,102,Region);
    }
    if (prev != -1)
        MPI_Recv(&(mat[0][1]),grid->delj,MPI_DOUBLE,prev,102,Region,&st);
    /* Left and right */
    MPI_Cart_shift(Region,0,1,&prev,&next);
    if (prev != -1)
    {
        for (int i = 1; i <= grid->deli; i++)
            buf[i-1] = mat[i][offy];
        MPI_Send(buf,grid->deli,MPI_DOUBLE,prev,103,Region);
    }
    if (next != -1)
    {
        MPI_Recv(buf,grid->deli,MPI_DOUBLE,next,103,Region,&st);
        for (int i = 1; i <= grid->deli; i++)
        {
            mat[i][0] = buf[i-1];
            buf[i-1] = mat[i][grid->delj+offy];
        }
        MPI_Send(buf,grid->deli,MPI_DOUBLE,next,104,Region);
    }
    if (prev != -1)
    {
        MPI_Recv(buf,grid->deli,MPI_DOUBLE,prev,104,Region,&st);
        int y = grid->delj+2*offy-1;
        for (int i = 1; i <= grid->deli; i++)
            mat[i][y] = buf[i-1];
    }
    return;
}








