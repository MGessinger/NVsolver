#include "types.h"


int main (int argc, char **argv)
{
	MPI_Init(&argc,&argv);
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
	REAL **U = NULL, **V = NULL, **P = NULL;
	char *paramFile = NULL;
	char *imageFile = NULL;
	int output = SILENT;

	for (int i = 2; i < argc; i++)
	{
		if (argv[i][0] != '-')
		{
			output |= atoi(argv[i])*OUTPUT;
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
	lattice grid = runSimulation(&U,&V,&P,argv[1],paramFile,output);
	/* Destroy simulated grids */
	destroy2Dfield(U,grid.deli+3);
	destroy2Dfield(V,grid.deli+2);
	destroy2Dfield(P,grid.deli+2);
	MPI_Finalize();
	return 0;
}
