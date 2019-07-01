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
    printf("\n"); /*Add another empty line for better readability in the output */
    return;
}

void printVector(REAL* vector, int len)
{
    print1Dfield(vector,len);
}

void printMatrix(REAL** matrix, int rows, int cols)
{
    print2Dfield(matrix,rows,cols);
}

void write1Dfield(const char *fileName, REAL* field, int size)
{
    if (field == NULL)
    {
        printf("The field in invalid.\n");
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

void write2Dfield(const char* fileName, REAL** field, int sizeX, int sizeY)
{
    if (field == NULL)
    {
        printf("The field in invalid.\n");
        return;
    }
    FILE *out = fopen(fileName,"wb");
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
        return NULL;
    }
    REAL *field = create1Dfield(*size);
    if (field == NULL)
    {
        printf("Could not allocate memory. Please try again.\n");
        return NULL;
    }
    int read = fread(field,sizeof(REAL),*size,in);
    if (read != *size)
    {
        printf("Only %i out of %i values were read correctly.\n",read,*size);
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
        return NULL;
    }
    REAL **field = create2Dfield(*sizeX,*sizeY);
    if (field == NULL)
    {
        printf("Could not allocate memory. Please try again.\n");
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

void writeVTKfileFor2DintegerField(const char* fileName, const char* description, short** field,
                                   int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Scalar Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", sizeX, sizeY);
    fprintf(vtkFile, "X_COORDINATES %d double\n", sizeX);
    for (int i=0;i<sizeX;i++)
        fprintf(vtkFile, "%lf ", dx*(double)i);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Y_COORDINATES %d double\n", sizeY);
    for (int j=0;j<sizeY;j++)
        fprintf(vtkFile, "%lf ", dy*(double)j);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Z_COORDINATES 1 double\n");
    fprintf(vtkFile, "0.0\n");

    fprintf(vtkFile, "POINT_DATA %d\n", 1 * sizeX * sizeY);
    fprintf(vtkFile, "SCALARS %s double 1\n", description);
    fprintf(vtkFile, "LOOKUP_TABLE default \n");
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%i\n", (field[i][j]!=0));
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field,
                                  int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Scalar Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", sizeX, sizeY);
    fprintf(vtkFile, "X_COORDINATES %d double\n", sizeX);
    for (int i=0;i<sizeX;i++)
        fprintf(vtkFile, "%lf ", dx*(double)i);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Y_COORDINATES %d double\n", sizeY);
    for (int j=0;j<sizeY;j++)
        fprintf(vtkFile, "%lf ", dy*(double)j);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Z_COORDINATES 1 double\n");
    fprintf(vtkFile, "0.0\n");

    fprintf(vtkFile, "POINT_DATA %d\n", 1 * sizeX * sizeY);
    fprintf(vtkFile, "SCALARS %s double 1\n", description);
    fprintf(vtkFile, "LOOKUP_TABLE default \n");
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%e\n", field[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV,
                                  int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = fopen(fileName, "w");
    if (vtkFile == NULL)
        return;

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Vector Field\n");
    fprintf(vtkFile, "ASCII\n");

    fprintf(vtkFile, "DATASET RECTILINEAR_GRID \n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1 \n", sizeX, sizeY);
    fprintf(vtkFile, "X_COORDINATES %d double\n", sizeX);
    for (int i=0;i<sizeX;i++)
        fprintf(vtkFile, "%lf ", dx*(double)i);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Y_COORDINATES %d double\n", sizeY);
    for (int j=0;j<sizeY;j++)
        fprintf(vtkFile, "%lf ", dy*(double)j);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "Z_COORDINATES 1 double\n");
    fprintf(vtkFile, "0.0\n");


    fprintf(vtkFile, "POINT_DATA %d\n", 1 * sizeX * sizeY);
    fprintf(vtkFile, "VECTORS %s double \n", description);
    for(int j=0;j<sizeY;j++)
    {
        for(int i=0;i<sizeX;i++)
        {
            fprintf(vtkFile, "%e %e 0.0\n", fieldU[i][j], fieldV[i][j]);
        }
    }
    fprintf(vtkFile,"\n");

    fclose(vtkFile);
}

void WriteParticle (particle *parts, int partcount, int n)
{
    if (parts == NULL || partcount == 0)
        return;
    char fileName[128];
    if (n != 0)
        sprintf(fileName, "ParticleField_%i.vtk",n);
    else
        sprintf(fileName, "ParticleField.vtk");
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
        if (parts[p].onScreen == 0)
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
   return(!png_sig_cmp(buf, 0, 8));
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

short** readGeometry (const char *flagFile, int *width, int *height)
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
    if (height == 0 || width == 0)
        return;
    /* Find the structure size in y-direction */
    short blockType, blockSize = height;
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
                structureVec[j] = 1;
            blockType = FLAG[i][j];
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
    structureVec = malloc(width*sizeof(short));
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
    *jmax = newWidth;
    *imax = newHeight;
    return;
}

int readParameters(const char *inputFile, REAL ***U, REAL ***V, REAL ***P,
                   lattice *grid, fluidSim *fluid, boundaryCond **bCond,
                   REAL *delt, REAL *t_end, char *problem)
{
    FILE *input = fopen(inputFile,"r");
    if (input == NULL)
    {
        printf("The requested file probably doesn't exist. Please check your input!\n");
        return 0;
    }
    char variableType[128];
    REAL value;
    REAL xlength = 0, ylength = 0;
    REAL UI = 0, VI = 0, PI = 0;
    int height = 0, width = 0;
    int readVars = 0;
    if (fscanf(input,"%[^\n\r]\n",problem) == 0)
        printf("The problem could not be detected. Assuming trivial fluid.\n");
    else
        printf("Initialising problem \"%s\"\n",problem);
    while (fscanf(input,"%s%*[^0-9-]%lg\n",variableType,&value) == 2)
    {
        switch(variableType[0])
        {
        case 'x':
            xlength = value;
            break;
        case 'y':
            ylength = value;
            break;
        case 'i':
            if (variableType[1] == 'm') grid->imax = (int)value;
            else fluid->itmax = (int)value;
            break;
        case 'j':
            grid->jmax = (int)value;
            break;
        case 'e':
            fluid->eps = value;
            break;
        case 'o':
            fluid->omega = value;
            break;
        case 'a':
            fluid->alpha = value;
            break;
        case 'd':
            *delt = value;
            break;
        case 't':
            if (variableType[1] == '_') *t_end = value;
            else fluid->tau = value;
            break;
        case 'R':
            fluid->Re = value;
            break;
        case 'U':
            UI = value;
            break;
        case 'V':
            VI = value;
            break;
        case 'P':
            PI = value;
            break;
        case 'G':
            if (variableType[1] == 'X') fluid->GX = value;
            else fluid->GY = value;
            break;
        default:
            printf("Found unexpected Variable %s of value %lg.\n",variableType,value);
            readVars--;
        }
        readVars++;
    }
    fclose(input);

    *bCond = createBoundCond(NOSLIP,OUTFLOW,NOSLIP,NOSLIP);
    short **image = readGeometry(variableType,&width,&height);
    findOptimalFlags(image,height,width,&(grid->imax),&(grid->jmax));
    (*bCond)->FLAG = adjustFlags(image,height,width,grid->imax,grid->jmax);
    initFlags("Image",(*bCond)->FLAG,grid->imax,grid->jmax);
    initUVP(U,V,P,grid->imax,grid->jmax,UI,VI,PI);
    grid->delx = xlength/grid->imax;
    grid->dely = ylength/grid->jmax;
    return readVars;
}

void outputVec(REAL **U, REAL **V, REAL **P, particle *parts, lattice *grid, int partcount, int n)
{
    if (grid == NULL)
        return;
    int imax = grid->imax;
    int jmax = grid->jmax;
    REAL **S = createMatrix(imax,jmax);
    REAL **T = createMatrix(imax,jmax);
    char fileName[64];
    if (P != NULL)
    {
        if (n != 0) sprintf(fileName,"PressureField_%i.vtk",n);
        else sprintf(fileName,"PressureField.vtk");
        for (int i = 0;  i < imax; i++)
            for (int j = 0; j < jmax; j++)
                T[i][j] = P[i+1][j+1];
        writeVTKfileFor2DscalarField(fileName,"pressurefield",T,imax,jmax,grid->delx,grid->dely);
    }
    if (U == NULL || V == NULL)
    {
        destroyMatrix(T,imax);
        destroyMatrix(S,imax);
        return;
    }
    if (n != 0) sprintf(fileName,"MomentumField_%i.vtk",n);
    else sprintf(fileName,"MomentumField.vtk");
    for (int i = 1;  i <= imax; i++)
        for (int j = 1; j <= jmax; j++)
        {
            T[i-1][j-1] = (U[i][j] + U[i-1][j])/2;
            S[i-1][j-1] = (V[i][j] + V[i][j-1])/2;
        }
    writeVTKfileFor2DvectorField(fileName,"momentumfield",T,S,imax,jmax,grid->delx,grid->dely);
    WriteParticle(parts,partcount,n);
    destroyMatrix(T,imax);
    destroyMatrix(S,imax);
    return;
}
