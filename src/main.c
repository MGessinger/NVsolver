#include "types.h"

int main (int argc, char **argv)
{
    if (argc < 3 || argv[1][0] == '-')
    {
        printf("Usage: simulator <scene>\n"
               "                 [-p\"parameter_file\"]\n"
               "                 [-i\"image_file]\"\n"
               "                 [number_of_frames]\n");
        printf("When specifying both an image and a parameter file, "
               "the sizes specified from the image take precedence!\n");
        printf("Type simulator -h for extended help.\n");
        return 0;
    }
    boundaryCond bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
    lattice grid;
    fluidSim sim;
    int rank, dims[2];
    REAL **U = NULL, **V = NULL, **P = NULL;
    REAL init[3];
    REAL delt, t_end;
    int out = SILENT;

    MPI_Init(&argc,&argv);
    MPI_Comm Region = createCommGrid(&rank,dims);
    for (int i = 2; i < argc; i++)
    {
        if (argv[i][0] != '-')
            out |= atoi(argv[i])*OUTPUT;
        else if (argv[i][1] == 'v')
            out |= PRINT;
        else if (argv[i][1] == 'p')
        {
            /* Read parameters from a file */
            if (readParameters(argv[i]+2,init,&grid,&sim,&bCond,&delt,&t_end) < 17)
            {
                printf("The parameter file appears to be incomplete.\n");
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
    splitRegion(Region, dims, &grid);
    initUVP(&U,&V,&P,grid.deli,grid.delj,init);
    if (bCond.FLAG == NULL)
        bCond.FLAG = create2DIntegerField(grid.deli+2,grid.delj+2);
    initFlags(argv[1],bCond.FLAG,&grid,Region);
    int files = simulateFluid(U,V,P,&bCond,&grid,&sim,Region,t_end,argv[1],out);
    if (out > OUTPUT)
        translateBinary(Region,&grid,files,rank,dims);
    /* Destroy simulated grids */
    destroy2Dfield(U,grid.deli+3);
    destroy2Dfield(V,grid.deli+2);
    destroy2Dfield(P,grid.deli+2);
    destroy2DIntegerField(bCond.FLAG,grid.deli+2);
    MPI_Finalize();
    return 0;
}

















