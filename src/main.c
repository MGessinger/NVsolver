#include "types.h"

void writeHorizontal (REAL **U, int sizeX, int sizeY, char *fileName)
{
	FILE *out = fopen(fileName,"w");
	fprintf(out,"[");
	for (int j = 1; j <= sizeY; j++)
		fprintf(out,"%g,",U[2+sizeX/2][j]);
	fprintf(out,"]");
	fclose(out);
	return;
}

int main (int argc, char **argv)
{
	MPI_Init(&argc,&argv);
	if (argc < 3)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		if (rank == 0)
		{
			printf("Usage: simulator <scene>\n"
					"                 [-p\"parameter_file\"]\n"
					"                 [-i\"image_file]\"\n"
					"                 [number_of_frames]\n");
			printf("When specifying both an image and a parameter file, "
					"the sizes specified from the image take precedence!\n");
		}
		MPI_Finalize();
		return 0;
	}
	REAL **U = NULL, **V = NULL, **P = NULL;
	char *paramFile = NULL;
	char *imageFile = NULL;
	char *scene = NULL;
	int output = SILENT;

	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] != '-')
		{
			char *endptr = "\0";
			int out = strtol(argv[i],&endptr,0);
			if (endptr[0] == '\0')
				output |= out*OUTPUT;
			else
				scene = argv[i];
			continue;
		}
		switch (argv[i][1])
		{
			case 'v':
				output |= PRINT;
				break;
			case 'p':
				paramFile = argv[i]+2;
				break;
			case 'i':
				imageFile = argv[i]+2;
				break;
			default:
				printf("Unrecognized option %s\n",argv[i]);
				break;
		}
	}
	lattice grid = runSimulation(&U,&V,&P,scene,paramFile,imageFile,output);
	//writeHorizontal(U,grid.deli,grid.delj,"hor.dat");
	/* Destroy simulated grids */
	destroy2Dfield((void**)U,grid.deli+3);
	destroy2Dfield((void**)V,grid.deli+2);
	destroy2Dfield((void**)P,grid.deli+2);
	MPI_Finalize();
	return 0;
}
