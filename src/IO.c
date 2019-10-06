#include "IO.h"

void print1Dfield(REAL* field, int size)
{
    if (field == NULL)
    {
        printf("The field is invalid. Please check your input!\n");
        return;
    }
    for (int i = 0; i < size; i++)
    {
        printf("%g\n",field[i]);
    }
    return;
}

void print2Dfield(REAL** field, int sizeX, int sizeY)
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

void write1Dfield(const char *fileName, REAL* field, int size)
{
    if (field == NULL)
    {
        printf("The field is invalid.\n");
        return;
    }
    FILE *out = fopen(fileName,"wb");
    if (out == NULL)
    {
        printf("The file %s could not be opened.\n",fileName);
        return;
    }
    fwrite(&size,sizeof(int),1,out);
    fwrite(field,sizeof(REAL),size,out);
    fclose(out);
    return;
}

void write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY, const char *mode)
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
    fwrite(&sizeX,sizeof(int),1,out);
    fwrite(&sizeY,sizeof(int),1,out);
    for (int i = 0; i < sizeX; i++)
    {
        fwrite(field[i],sizeof(REAL),sizeY,out);
    }
    fclose(out);
    return;
}

REAL* read1Dfield(const char* fileName, int* size)
{
    if (size == NULL)
    {
        printf("There is no way to save the size. Please provide a pointer.\n");
        return NULL;
    }
    FILE *in = fopen(fileName,"rb");
    if (in == NULL)
    {
        printf("The file %s could not be read.\n",fileName);
        return NULL;
    }
    /* The file has been opened successfully so data can be stored */
    if (fread(size,sizeof(int),1,in) != 1)
    {
        printf("The file has incorrect format.\n");
        fclose(in);
        return NULL;
    }
    REAL *field = create1Dfield(*size);
    if (field == NULL)
    {
        printf("Could not allocate memory. Please try again.\n");
        fclose(in);
        return NULL;
    }
    int read = fread(field,sizeof(REAL),*size,in);
    if (read != *size)
    {
        printf("[WARNING] Only %i out of %i values were read correctly.\n",read,*size);
    }
    fclose(in);
    return field;
}

REAL** read2Dfield(const char *fileName, int *sizeX, int *sizeY)
{
    if (sizeX == 0 || sizeY == NULL)
    {
        printf("Cannot save the size of the array. Please check the pointers.\n");
        return NULL;
    }
    FILE *in = fopen(fileName,"rb");
    if (in == NULL)
    {
        printf("The file %s could not be read.\n",fileName);
        return NULL;
    }
    /* The file has been opened successfully so data can be stored */
    if (fread(sizeX,sizeof(int),1,in) == 0 || fread(sizeY,sizeof(int),1,in) == 0)
    {
        printf("The file has incorrect format.\n");
        fclose(in);
        return NULL;
    }
    REAL **field = create2Dfield(*sizeX,*sizeY);
    if (field == NULL)
    {
        printf("Could not allocate memory. Please try again.\n");
        fclose(in);
        return NULL;
    }
    int read = 0;
    for (int i = 0; i < *sizeY; i++)
    {
        read += fread(field[i],sizeof(REAL),*sizeY,in);
    }
    if (read != (*sizeX)*(*sizeY))
    {
        printf("Only %i out of %i variables were read correctly.\n",read,(*sizeX)*(*sizeY));
    }
    fclose(in);
    return field;
}

void writeVTKfileFor2DintegerField(const char* fileName, const char* description, short** field, lattice *grid)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Scalar Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", grid->imax, grid->jmax);
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
    fprintf(vtkFile, "LOOKUP_TABLE default \n");
    for(int j=0;j<grid->delj;j++)
    {
        for(int i=0;i<grid->deli;i++)
        {
            fprintf(vtkFile, "%i\n", (field[i][j]!=C_F));
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, lattice *grid)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Scalar Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", grid->imax, grid->jmax);
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
    fprintf(vtkFile, "LOOKUP_TABLE default \n");
    for(int j = 0; j < grid->jmax; j++)
    {
        for(int i = 0; i < grid->imax; i++)
        {
            fprintf(vtkFile, "%e\n", field[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DvectorField(const char* fileName, const char* description,
                                  REAL** fieldU, REAL** fieldV, lattice *grid)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Vector Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", grid->imax, grid->jmax);
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
    fprintf(vtkFile, "VECTORS %s double \n", description);
    for(int j = 0; j < grid->jmax; j++)
    {
        for(int i = 0; i < grid->imax; i++)
        {
            fprintf(vtkFile, "%e %e 0.0\n", fieldU[i][j], fieldV[i][j]);
        }
    }
    fprintf(vtkFile,"\n");
    fclose(vtkFile);
}

void WriteParticle (particle *parts, int partcount, int n)
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

int check_if_png(const char *fileName, FILE **file)
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

short** readGeometry (const char *flagFile, int *height, int *width)
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

    short **FLAG = create2DIntegerField(*width,*height);
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

short** adjustFlags(short **FLAG, int height, int width, int imax, int jmax)
{
    /* Adjust the number of cells to a predefined number */
    if (FLAG == NULL)
        return FLAG;
    short **newFLAG = create2DIntegerField(imax,jmax);
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
    destroy2DIntegerField(FLAG,width);
    return newFLAG;
}

void findOptimalFlags(short **FLAG, int height, int width, int *imax, int *jmax)
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

int readParameters(const char *inputFile, REAL *init,
                   lattice *grid, fluidSim *sim, boundaryCond *bCond,
                   REAL *delt, REAL *t_end)
{
    FILE *input = fopen(inputFile,"r");
    if (input == NULL)
    {
        printf("The file \"%s\" probably doesn't exist. Please check your input!\n",inputFile);
        return 0;
    }
    char variable[32];
    REAL value;
    REAL xlength = 0, ylength = 0;
    int readVars = 0;

    while (fscanf(input,"%s%*[A-Za-z ]%lg\n",variable,&value) == 2)
    {
        if (variable[0] == '#')
            continue;
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
            if (variable[1] == 't') bCond->wt = value;
            else if (variable[1] == 'b') bCond->wb = value;
            else if (variable[1] == 'r') bCond->wr = value;
            else if (variable[1] == 'l') bCond->wl = value;
            break;
        default:
            printf("Found unexpected Variable %s of value %lg.\n",variable,value);
            continue;
        }
        readVars++;
    }
    fclose(input);
    grid->delx = xlength/grid->imax;
    grid->dely = ylength/grid->jmax;
    return readVars;
}

void outputVec(REAL **U, REAL **V, REAL **P, lattice *grid, int n)
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
            T[i-1][j-1] = (U[i][j] + U[i-1][j])/2;
            S[i-1][j-1] = (V[i][j] + V[i][j-1])/2;
        }
    writeVTKfileFor2DvectorField(fileName,"momentumfield",T,S,grid);
    destroy2Dfield(T,imax);
    destroy2Dfield(S,imax);
    return;
}

int dumpFields(REAL **U, REAL **V, REAL **P, lattice *grid, int n)
{
    int rank, size, send = 1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status st;
    char pFile[32], uFile[32], vFile[32];
    char *mode = (rank == 0 ? "wb" : "ab");
    sprintf(pFile,"P%i",n);
    sprintf(uFile,"U%i",n);
    sprintf(vFile,"V%i",n);
    if (rank > 0)
        MPI_Recv(&send,1,MPI_INT,rank-1,WRITE + 0, MPI_COMM_WORLD, &st);
    else
        printf("%s, %s, %s\n",pFile,uFile,vFile);
    write2Dfield(pFile,P,grid->deli+2,grid->delj+2,mode);
    if (rank+1 != size)
        MPI_Send(&send,1,MPI_INT,rank+1,WRITE + 0, MPI_COMM_WORLD);
    if (rank > 0)
        MPI_Recv(&send,1,MPI_INT,rank-1,WRITE + 1,MPI_COMM_WORLD,&st);
    write2Dfield(uFile,U,grid->deli+2,grid->delj+2,mode);
    if (rank+1 != size)
        MPI_Send(&send,1,MPI_INT,rank+1,WRITE + 1, MPI_COMM_WORLD);
    if (rank > 0)
        MPI_Recv(&send,1,MPI_INT,rank-1,WRITE + 2,MPI_COMM_WORLD,&st);
    write2Dfield(vFile,V,grid->deli+2,grid->delj+2,mode);
    if (rank+1 != size)
        MPI_Send(&send,1,MPI_INT,rank+1,WRITE + 2,MPI_COMM_WORLD);
    return send;
}

void translateBinary (MPI_Comm Region, lattice *grid, int files, int rank, int *dims)
{
    char pFile[32], uFile[32], vFile[32];
    REAL **U, **V, **P;
    REAL size[2], init[3] = {0,0,0};
    initUVP(&U,&V,&P,grid->imax,grid->jmax,init);
    FILE *PF, *UF, *VF;
    int coords[2], nproc = dims[0]*dims[1];
    int il[nproc], jb[nproc];
    for (int k = 0; k < nproc; k++)
    {
        MPI_Cart_coords(Region,k,2,coords);
        il[k] = coords[0]*grid->imax/dims[0];
        jb[k] = coords[1]*grid->jmax/dims[1];
    }
    for (int i = rank; i < files; i += nproc) /* Loop over files */
    {
        sprintf(pFile,"P%i",i);
        sprintf(uFile,"U%i",i);
        sprintf(vFile,"V%i",i);
        PF = fopen(pFile,"rb");
        UF = fopen(uFile,"rb");
        VF = fopen(vFile,"rb");
        if (!PF || !UF || !VF)
            continue;
        for (int k = 0; k < nproc; k++) /* Loop over matrices */
        {
            if (fread(size,sizeof(REAL),2,PF) == 0)
                break;
            for (int j = 0; j < size[0]; j++) /* Loop over lines */
            {
                if (fread(&(P[j+il[k]][jb[k]]),sizeof(REAL),size[1],PF) == 0)
                    break;
                if (fread(&(U[j+il[k]][jb[k]]),sizeof(REAL),size[1],UF) == 0)
                    break;
                if (fread(&(V[j+il[k]][jb[k]]),sizeof(REAL),size[1],VF) == 0)
                    break;
            }
        }
        fclose(PF);
        fclose(UF);
        fclose(VF);
        outputVec(U,V,P,grid,i);
    }
    return;
}
