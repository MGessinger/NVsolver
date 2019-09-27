#include "types.h"

lattice *createLattice()
{
    lattice *grid = malloc(sizeof(lattice));
    if (grid == NULL)
        return NULL;
    grid->imax = grid->il = 0;
    grid->jmax = grid->jb = 0;
    return grid;
}

int main (int argc, char **argv)
{
    if (argc < 3)
    {
        printf("Usage: simulator [-p \"parameter_file\"]\n"
               "                 [-i \"image_file]\"\n"
               "                 [number_of_frames]\n");
        printf("When specifying both an image and a parameter file, "
               "the sizes specified from the image take precedence!\n");
        return 0;
    }
    REAL **U = NULL, **V = NULL, **P = NULL;
    REAL init[3];
    REAL delt, t_end;
    int out = SILENT;

    boundaryCond *bCond = createBoundCond(NOSLIP,NOSLIP,NOSLIP,NOSLIP);
    lattice *grid = createLattice();
    fluidSim sim;
    char problem[128];
    for (int i = 0; i < argc; i++)
    {
        if (argv[i][0] != '-')
        {
            out = PRINT | atoi(argv[i])*OUTPUT;
            continue;
        }
        if (argv[i][1] == 'p' && i+1 < argc)
        {
            /* Read prameters from a file */
            if (readParameters(argv[i+1],init,grid,&sim,bCond,&delt,&t_end,problem) < 17)
            {
                destroyBoundCond(bCond,grid->imax);
                free(grid);
                return 0;
            }
        }
        else if (argv[i][1] == 'i' && i+1 < argc)
        {
            strcpy(problem,"Image");
            int width, height;
            bCond->FLAG = readGeometry(argv[i+1],&width,&height);
            /*findOptimalFlags(image,height,width,&(grid->imax),&(grid->jmax));
            bCond->FLAG = adjustFlags(image,height,width,grid->imax,grid->jmax);*/
        }
    }
    initUVP(&U,&V,&P,grid->imax,grid->jmax,init);
    if (bCond->FLAG == NULL)
    {
        bCond->FLAG = create2DIntegerField(grid->imax,grid->jmax);
        initFlags(problem,bCond->FLAG,grid->imax,grid->jmax);
    }
    simulateFluid(U,V,P,bCond,grid,&sim,delt,t_end,problem,out);
    outputVec(U,V,P,grid,0);
    /* Destroy simulated grids */
    if (grid != NULL)
    {
        destroy2Dfield(U,grid->imax+2);
        destroy2Dfield(V,grid->imax+2);
        destroy2Dfield(P,grid->imax+2);
        destroyBoundCond(bCond,grid->imax);
        free(grid);
    }
    return 0;
}
