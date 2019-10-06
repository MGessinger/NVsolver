#include "types.h"

lattice *createLattice()
{
    lattice *grid = malloc(sizeof(lattice));
    if (grid == NULL)
        return NULL;
    grid->imax = grid->jmax = 0;
    grid->edges = 0;
    return grid;
}

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
    int rank, size = dims[0]*dims[1];
    int coords[2];
    MPI_Comm_rank(Region,&rank);
    MPI_Cart_coords(Region,rank,size,coords);
    grid->il = coords[0]*grid->imax/dims[0];
    if (grid->il == 0)
        grid->edges |= LEFT;
    int ir = (coords[0]+1)*grid->imax/dims[0];
    if (ir == grid->imax)
        grid->edges |= RIGHT;
    grid->deli = ir - grid->il;
    grid->jb = coords[1]*grid->jmax/dims[1];
    if (grid->jb == 0)
        grid->edges |= BOTTOM;
    int jt = (coords[1]+1)*grid->jmax/dims[1];
    if (jt == grid->jmax)
        grid->edges |= TOP;
    grid->delj = jt - grid->jb;
    printf("%i ]\t%i,%i|%i,%i\n",rank,grid->il,ir,grid->jb,jt);
    if (bCond->FLAG == NULL)
        bCond->FLAG = create2DIntegerField(grid->deli,grid->delj);
    initFlags(problem,bCond->FLAG,grid);
}

int main (int argc, char **argv)
{
    if (argc < 3 || argv[1][0] == '-')
    {
        printf("Usage: simulator <scene>\n"
               "                 [-p\"parameter_file\"]\n"
             /*"                 [-i\"image_file]\"\n"*/
               "                 [number_of_frames]\n");
        printf("When specifying both an image and a parameter file, "
               "the sizes specified from the image take precedence!\n");
        printf("Type simulator -h for extended help.\n");
        return 0;
    }
    int rank, dims[2];
    MPI_Init(&argc,&argv);
    MPI_Comm Region = createCommGrid(&rank,dims);
    REAL **U = NULL, **V = NULL, **P = NULL;
    REAL init[3];
    REAL delt, t_end;
    int out = SILENT;

    boundaryCond *bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
    lattice *grid = createLattice();
    fluidSim sim;
    char problem[128];
    strcpy(problem,argv[1]);
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] != '-')
        {
            out = PRINT | atoi(argv[i])*OUTPUT;
            continue;
        }
        if (argv[i][1] == 'p')
        {
            /* Read prameters from a file */
            if (readParameters(argv[i]+2,init,grid,&sim,bCond,&delt,&t_end) < 17)
            {
                printf("The parameter file appears to be incomplete.\n");
                destroyBoundCond(bCond,grid->imax);
                free(grid);
                MPI_Abort(Region,0);
            }
        }
        /*else if (argv[i][1] == 'i')
        {
            strcpy(problem,"Image");
            int width, height;
            bCond->FLAG = readGeometry(argv[i]+2,&width,&height);
            findOptimalFlags(image,height,width,&(grid->imax),&(grid->jmax));
            bCond->FLAG = adjustFlags(image,height,width,grid->imax,grid->jmax);
        }*/
    }
    /* Slice Data to process */
    splitRegion(Region, dims, grid, problem, bCond);
    initUVP(&U,&V,&P,grid->deli,grid->delj,init);
    int files = simulateFluid(U,V,P,bCond,grid,&sim,t_end,problem,out);
    translateBinary(Region,grid,files,rank,dims);
    /* Destroy simulated grids */
    if (grid != NULL)
    {
        destroy2Dfield(U,grid->deli+2);
        destroy2Dfield(V,grid->deli+2);
        destroy2Dfield(P,grid->deli+2);
        destroyBoundCond(bCond,grid->deli);
        free(grid);
    }
    MPI_Finalize();
    return 0;
}
