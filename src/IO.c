#include "types.h"

void print2Dfield (REAL** field, int sizeX, int sizeY)
{
	if (field == NULL)
	{
		printf("The field is invalid. Please check your input!\n");
		return;
	}
	for (int j = 0; j < sizeY; j++)
	{
		for (int i = 0; i < sizeX; i++)
		{
			printf("%5g\t",field[i][j]);
		}
		printf("\n");
	}
	printf("\n"); /* Add another empty line for better readability in the output */
	return;
}

void write2Dfield (const char* fileName, REAL** field, size_t sizeX, size_t sizeY, size_t offX, size_t offY, const char *mode)
{
	if (field == NULL)
	{
		printf("The field is invalid.\n");
		return;
	}
	FILE *out = fopen(fileName,mode);
	if (out == NULL)
	{
		printf("The file %s could not be opened.\n",fileName);
		return;
	}
	fwrite(&sizeX,sizeof(size_t),1,out);
	fwrite(&sizeY,sizeof(size_t),1,out);
	for (size_t i = 0; i < sizeX; i++)
	{
		fwrite(field[i+offX]+offY,sizeof(REAL),sizeY,out);
	}
	fclose(out);
	return;
}

int read2Dfield (FILE *in, REAL **field, size_t offX, size_t offY)
{
	if (!in || !field)
		return -1;
	size_t size[2];
	if (fread(size,sizeof(size_t),2,in) < 2)
		return -1;
	for (size_t k = 0; k < size[0]; k++) /* Loop over lines */
	{
		if (fread(field[k+offX]+offY,sizeof(REAL),size[1],in) != size[1])
			return -1;
	}
	return 1;
}

void writeVTKfileFor2DintegerField (const char* fileName, const char* description, char **field, lattice *grid)
{
	FILE* vtkFile = fopen(fileName, "w");
	if (vtkFile == NULL)
		return;

	fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
	fprintf(vtkFile, "Scalar Field\n");
	fprintf(vtkFile, "ASCII\n");

	fprintf(vtkFile, "DATASET RECTILINEAR_GRID\n");
	fprintf(vtkFile, "DIMENSIONS %d %d 1\n", grid->imax, grid->jmax);
	fprintf(vtkFile, "X_COORDINATES %d double\n", grid->imax);
	for (REAL i = 0; i < grid->imax; i++)
		fprintf(vtkFile, "%lf ", grid->delx*i);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Y_COORDINATES %d double\n", grid->jmax);
	for (REAL j = 0; j < grid->jmax; j++)
		fprintf(vtkFile, "%lf ", grid->dely*j);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Z_COORDINATES 1 double\n");
	fprintf(vtkFile, "0.0\n");

	fprintf(vtkFile, "POINT_DATA %d\n", grid->imax * grid->jmax);
	fprintf(vtkFile, "SCALARS %s double 1\n", description);
	fprintf(vtkFile, "LOOKUP_TABLE default\n");
	for(int j = 0; j < grid->delj; j++)
	{
		for(int i = 0; i < grid->deli; i++)
		{
			fprintf(vtkFile, "%i\n", (field[i][j]!=C_F));
		}
	}
	fprintf(vtkFile,"\n");
	fclose(vtkFile);
	return;
}

void writeVTKfileFor2DscalarField (const char* fileName, const char* description, REAL** field, lattice *grid)
{
	FILE* vtkFile = fopen(fileName, "w");
	if (vtkFile == NULL)
		return;

	fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
	fprintf(vtkFile, "Scalar Field\n");
	fprintf(vtkFile, "ASCII\n");

	fprintf(vtkFile, "DATASET RECTILINEAR_GRID\n");
	fprintf(vtkFile, "DIMENSIONS %d %d 1\n", grid->imax, grid->jmax);
	fprintf(vtkFile, "X_COORDINATES %d double\n", grid->imax);
	for (REAL i = 0; i < grid->imax; i++)
		fprintf(vtkFile, "%lf ", grid->delx*i);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Y_COORDINATES %d double\n", grid->jmax);
	for (REAL j = 0; j < grid->jmax; j++)
		fprintf(vtkFile, "%lf ", grid->dely*j);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Z_COORDINATES 1 double\n");
	fprintf(vtkFile, "0.0\n");

	fprintf(vtkFile, "POINT_DATA %d\n", grid->imax * grid->jmax);
	fprintf(vtkFile, "SCALARS %s double 1\n", description);
	fprintf(vtkFile, "LOOKUP_TABLE default\n");
	for(int j = 0; j < grid->jmax; j++)
	{
		for(int i = 0; i < grid->imax; i++)
		{
			fprintf(vtkFile, "%e\n", field[i][j]);
		}
	}
	fprintf(vtkFile,"\n");
	fclose(vtkFile);
	return;
}

void writeVTKfileFor2DvectorField (const char* fileName, const char* description,
		REAL** fieldU, REAL** fieldV, lattice *grid)
{
	FILE* vtkFile = fopen(fileName, "w");
	if (vtkFile == NULL)
		return;

	fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
	fprintf(vtkFile, "Vector Field\n");
	fprintf(vtkFile, "ASCII\n");

	fprintf(vtkFile, "DATASET RECTILINEAR_GRID\n");
	fprintf(vtkFile, "DIMENSIONS %d %d 1\n", grid->imax, grid->jmax);
	fprintf(vtkFile, "X_COORDINATES %d double\n", grid->imax);
	for (REAL i = 0; i < grid->imax; i++)
		fprintf(vtkFile, "%lf ", grid->delx*i);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Y_COORDINATES %d double\n", grid->jmax);
	for (REAL j = 0; j < grid->jmax; j++)
		fprintf(vtkFile, "%lf ", grid->dely*j);
	fprintf(vtkFile, "\n");
	fprintf(vtkFile, "Z_COORDINATES 1 double\n");
	fprintf(vtkFile, "0.0\n");

	fprintf(vtkFile, "POINT_DATA %d\n", grid->imax * grid->jmax);
	fprintf(vtkFile, "VECTORS %s double\n", description);
	for(int j = 0; j < grid->jmax; j++)
	{
		for(int i = 0; i < grid->imax; i++)
		{
			fprintf(vtkFile, "%e %e 0.0\n", fieldU[i][j], fieldV[i][j]);
		}
	}
	fprintf(vtkFile,"\n");
	fclose(vtkFile);
	return;
}

void writeParticle (particle *parts, int partcount, int n)
{
	if (!parts|| !partcount)
		return;
	char fileName[128];
	sprintf(fileName, "ParticleField_%i.vtk",n);
	FILE *out = fopen(fileName,"w");
	if (out == NULL)
		return;
	fprintf(out,"# vtk DataFile Version 3.0\n");
	fprintf(out,"Partikel\n");
	fprintf(out,"ASCII\n");

	fprintf(out,"DATASET POLYDATA\n");
	fprintf(out,"POINTS %i double\n",partcount);
	for (int p = 0; p < partcount; p++)
	{
		if (!parts[p].onScreen)
			fprintf(out,"0.1 0.1 -0.1\n");
		else
			fprintf(out,"%e %e 0.0\n",parts[p].x,parts[p].y);
	}
	fclose(out);
}

int check_if_png (const char *fileName, FILE **file)
{
	unsigned char buf[8];
	/* Open the alleged PNG file. */
	if ((*file = fopen(fileName, "rb")) == NULL)
		return NOT_PNG;

	/* Read the first eight bytes and compare them to the signature */
	if (fread(buf, 1, 8, *file) != 8)
		return NOT_PNG;
	return !png_sig_cmp(buf, 0, 8);
}

void readImageData (FILE *flagData, png_structpp png_ptr, png_infopp info_ptr)
{
	if (flagData == NULL || png_ptr == NULL || info_ptr == NULL)
		return;
	/* Initialization */
	*png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
	if (*png_ptr == NULL)
		return;
	*info_ptr = png_create_info_struct(*png_ptr);
	if (*info_ptr == NULL)
	{
		png_destroy_read_struct(png_ptr, NULL, NULL);
		return;
	}
	if (setjmp(png_jmpbuf(*png_ptr)))
	{
		/* Free all of the memory associated with the png_ptr and info_ptr */
		png_destroy_read_struct(png_ptr, info_ptr, NULL);
		/* If we get here, we had a problem reading the file. */
		printf("The has been an issue when reading the file. If only I knew...\n");
		return;
	}
	png_init_io(*png_ptr, flagData);

	/* Start reading data from the file (except for the signature, which has been read before) */
	png_set_sig_bytes(*png_ptr, 8);
	png_read_png(*png_ptr, *info_ptr, PNG_TRANSFORM_PACKING | PNG_TRANSFORM_GRAY_TO_RGB, NULL);
	fclose(flagData);
	return;
}

char** readGeometry (const char *flagFile, int *height, int *width)
{
	/* Variables */
	FILE *flagData;
	png_structp png_ptr;
	png_infop info_ptr;

	if (!check_if_png(flagFile,&flagData))
	{
		printf("The given file does not seem to be a png file. Please check your input.\n");
		return NULL;
	}
	printf("Reading flags from \"%s\"",flagFile);
	readImageData(flagData,&png_ptr,&info_ptr);

	*height = png_get_image_height(png_ptr,info_ptr);
	*width = png_get_image_width(png_ptr,info_ptr);

	png_bytepp rows;
	rows = png_get_rows(png_ptr,info_ptr);

	char **FLAG = create2DIntegerField(*width,*height);
	for (int i = 0; i < *height; i++)
	{
		for (int j = 0; j+2 < 3*(*width); j+=3)
		{
			if ( rows[i][j] < 0x90 || rows[i][j+1] < 0x90 || rows[i][j+2] < 0x90)
				FLAG[j/3][*height-1-i] = C_B;
			else
				FLAG[j/3][*height-1-i] = C_F;
		}
	}

	png_destroy_read_struct(&png_ptr,&info_ptr,NULL);
	return FLAG;
}

char** adjustFlags (char **FLAG, int height, int width, int imax, int jmax)
{
	/* Adjust the number of cells to a predefined number */
	if (FLAG == NULL)
		return FLAG;
	char **newFLAG = create2DIntegerField(imax,jmax);
	if (newFLAG == NULL)
		return FLAG;
	short average;
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			average = 0;
			for (int k = (i*width)/imax; k <= ((i+1)*width)/imax; k++)
				for (int l = (j*height)/jmax; l <= ((j+1)*height)/jmax; l++)
				{
					if (k == width || l == height)
						continue;
					average += FLAG[k][l];
				}
			if (average > (width*height)/(imax*jmax*2))
				newFLAG[i][j] = C_B;
			else
				newFLAG[i][j] = C_F;
		}
	}
	destroy2Dfield((void**)FLAG,width);
	return newFLAG;
}

void findOptimalFlags (char **FLAG, int height, int width, int *imax, int *jmax)
{
	if (!FLAG || !imax || !jmax)
		return;
	if (!height || !width)
		return;
	/* Find the structure size in y-direction */
	int blockType, blockSize = height;
	int newHeight, newWidth;
	int i, j;
	short *structureVec = malloc(((short)height)*sizeof(short));
	if (structureVec == NULL)
		return;
	for (j = 0; j < height; j++)
		structureVec[j] = 0;
	for (i = 0; i < width; i++)
	{
		blockType = FLAG[i][0];
		for (j = 1; j < height; j++)
		{
			if (FLAG[i][j] != blockType)
			{
				structureVec[j] = 1;
				blockType = FLAG[i][j];
			}
		}
	}
	for (j = 0; j < height; j+=i)
	{
		i = 1;
		while (j+i < height && structureVec[j+i] == 0)
		{
			i++;
		}
		if (i < blockSize)
			blockSize = i;
		if (blockSize == 1)
			break;
	}
	free(structureVec);
	/* Extract the height necessary to represent the image, but no less than *jmax */
	newHeight = (height * 2)/blockSize;
	if (newHeight < *jmax)
		newHeight = *jmax;
	/* Now do exactly the same thing for the x-direction */
	blockSize = width;
	structureVec = malloc((short)width*sizeof(short));
	if (structureVec == NULL)
		return;
	for (i = 0; i < width; i++)
		structureVec[i] = 0;
	for (j = 0; j < height; j++)
	{
		blockType = FLAG[0][j];
		for (i = 1; i < width; i++)
		{
			if (FLAG[i][j] != blockType)
				structureVec[i] = 1;
			blockType = FLAG[i][j];
		}
	}
	for (i = 0; i < width; i+=j)
	{
		j = 1;
		while (j+i < width && structureVec[j+i] == 0)
		{
			j++;
		}
		if (j < blockSize)
			blockSize = j;
		if (blockSize == 1)
			break;
	}
	free(structureVec);
	/* Extract the width necessary to represent the image, but no less than *imax */
	newWidth = (width * 2)/blockSize;
	if (newWidth < *imax)
		newWidth = *imax;
	*jmax = newHeight;
	*imax = newWidth;
	return;
}

int readParameters (const char *inputFile, REAL *init,
		lattice *grid, fluidSim *sim, bndCond *bCond,
		REAL *delt, REAL *t_end)
{
	FILE *input = fopen(inputFile,"r");
	if (input == NULL)
	{
		printf("The file \"%s\" probably doesn't exist. Please check your input!\n",inputFile);
		return 0;
	}
	REAL value;
	REAL xlength = 0, ylength = 0;
	int readVars = 0;

	char *variable = NULL;
	size_t len = 0;

	while (getline(&variable, &len, input) > 0)
	{
		if (variable[0] == '#')
			continue;
		if (sscanf(variable, "%*s%*[ \t]%lg", &value) != 1)
		{
			printf("Malformed input: %s\n",variable);
			continue;
		}
		switch(variable[0])
		{
			case 'x':
				xlength = value;
				break;
			case 'y':
				ylength = value;
				break;
			case 'i':
				if (variable[1] == 'm') grid->imax = (int)value;
				else sim->itmax = (int)value;
				break;
			case 'j':
				grid->jmax = (int)value;
				break;
			case 'e':
				sim->eps = value;
				break;
			case 'o':
				sim->omega = value;
				break;
			case 'a':
				sim->alpha = value;
				break;
			case 'd':
				*delt = sim->dt = value;
				break;
			case 't':
				if (variable[1] == 'a') sim->tau = value;
				else *t_end = value;
				break;
			case 'R':
				sim->Re = value;
				break;
			case 'U':
				init[0] = value;
				break;
			case 'V':
				init[1] = value;
				break;
			case 'P':
				init[2] = value;
				break;
			case 'G':
				if (variable[1] == 'X') sim->GX = value;
				else sim->GY = value;
				break;
			case 'w':
				if (variable[1] == 't') bCond->wt = (unsigned)value;
				else if (variable[1] == 'b') bCond->wb = (unsigned)value;
				else if (variable[1] == 'r') bCond->wr = (unsigned)value;
				else if (variable[1] == 'l') bCond->wl = (unsigned)value;
				break;
			default:
				printf("Found unexpected Variable %s of value %lg.\n",variable,value);
				continue;
		}
		readVars++;
	}
	fclose(input);
	free(variable);
	grid->delx = xlength/grid->imax;
	grid->dely = ylength/grid->jmax;
	return readVars;
}

void outputVec (REAL **U, REAL **V, REAL **P, lattice *grid, int n)
{
	int imax = grid->imax;
	int jmax = grid->jmax;
	REAL **S = create2Dfield(imax,jmax);
	REAL **T = create2Dfield(imax,jmax);
	char fileName[32];
	sprintf(fileName,"PressureField_%i.vtk",n);
	for (int i = 0;  i < imax; i++)
		for (int j = 0; j < jmax; j++)
			T[i][j] = P[i+1][j+1];
	writeVTKfileFor2DscalarField(fileName,"pressurefield",T,grid);
	sprintf(fileName,"MomentumField_%i.vtk",n);
	for (int i = 1;  i <= imax; i++)
		for (int j = 1; j <= jmax; j++)
		{
			T[i-1][j-1] = (U[i+1][j] + U[i][j])/2;
			S[i-1][j-1] = (V[i][j+1] + V[i][j])/2;
		}
	writeVTKfileFor2DvectorField(fileName,"momentumfield",T,S,grid);
	destroy2Dfield((void**)T,imax);
	destroy2Dfield((void**)S,imax);
	return;
}

int dumpFields (MPI_Comm Region, REAL **U, REAL **V, REAL **P, lattice *grid, int n)
{
	int rank, size, send = 1;
	MPI_Comm_rank(Region,&rank);
	MPI_Comm_size(Region,&size);
	MPI_Status st;
	char pFile[32], uFile[32];
	char *mode = (rank == 0 ? "wb" : "ab");
	sprintf(pFile,"PressureField_%i.vtk",n);
	sprintf(uFile,"MomentumField_%i.vtk",n);
	if (rank > 0)
		MPI_Recv(&send,1,MPI_INT,rank-1,WRITE + 0, Region, &st);
	write2Dfield(pFile,P,grid->deli,grid->delj,1,1,mode);
	if (rank+1 != size)
		MPI_Send(&send,0,MPI_INT,rank+1,WRITE + 0, Region);
	if (rank > 0)
		MPI_Recv(&send,0,MPI_INT,rank-1,WRITE + 1,Region,&st);
	write2Dfield(uFile,U,grid->deli+1,grid->delj+1,1,0,mode);
	write2Dfield(uFile,V,grid->deli+1,grid->delj+1,0,1,"ab");
	if (rank+1 != size)
		MPI_Send(&send,0,MPI_INT,rank+1,WRITE + 1,Region);
	return send;
}

void translateBinary (MPI_Comm Region, lattice *grid, int files, int rank, int *dims)
{
	if (files <= rank)
		return;
	char pFile[32], uFile[32];
	FILE *PF, *UF;
	REAL **U, **V, **P;
	initUVP(&U,&V,&P,grid->imax,grid->jmax,NULL);
	int coords[2], nproc = dims[0]*dims[1];
	int il[nproc], jb[nproc];
	if (Region == MPI_COMM_WORLD)
		il[0] = jb[0] = 0;
	else
		for (int k = 0; k < nproc; k++)
		{
			MPI_Cart_coords(Region,k,2,coords);
			il[k] = coords[0]*grid->imax/dims[0];
			jb[k] = coords[1]*grid->jmax/dims[1];
		}
	for (int i = rank; i < files; i += nproc) /* Loop over files */
	{
		sprintf(pFile,"PressureField_%i.vtk",i);
		sprintf(uFile,"MomentumField_%i.vtk",i);
		PF = fopen(pFile,"rb");
		UF = fopen(uFile,"rb");
		for (int j = 0; j < nproc; j++) /* Loop over sub-matrices */
		{
			if (read2Dfield(PF,P,il[j]+1,jb[j]+1) < 0)
				printf("Skipping file %s...\n",pFile);
			if (read2Dfield(UF,U,il[j]+1,jb[j]) < 0)
				printf("Skipping file %s...\n",uFile);
			if (read2Dfield(UF,V,il[j],jb[j]+1) < 0)
				printf("Skipping file %s...\n",uFile);
		}
		fclose(PF);
		fclose(UF);
		outputVec(U,V,P,grid,i);
	}
	destroy2Dfield((void**)U,grid->imax+3);
	destroy2Dfield((void**)V,grid->imax+2);
	destroy2Dfield((void**)P,grid->imax+2);
	return;
}
