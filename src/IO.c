#include "IO.h"

FILE *open_file(const char *fileName, const char *mode)
{
    char filePath[256];
    #ifdef linux
    strcpy(filePath,"/media/");
    strcat(filePath,USER);
    strcat(filePath,"/STICK/");
    #else
    strcpy(filePath,"E:/");
    #endif
    strcat(filePath,"Programming/Praktikum/Blatt1/data/");
    strcat(filePath,fileName);
    return fopen(filePath,mode);
}

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
            printf("%g\t",field[i][j]);
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
    FILE *out = open_file(fileName,"wb");
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
    FILE *out = open_file(fileName,"wb");
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
    FILE *in = open_file(fileName,"rb");
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

REAL** read2Dfield(const char* fileName, int* sizeX, int* sizeY)
{
    if (sizeX == 0 || sizeY == NULL)
    {
        printf("Cannot save the size of the array. Please check the pointers.\n");
        return NULL;
    }
    FILE *in = open_file(fileName,"rb");
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


void writeVTKfileFor2DscalarField(const char* fileName, const char* description, REAL** field, int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = open_file(fileName, "w");
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

void writeVTKfileFor2DvectorField(const char* fileName, const char* description, REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy)
{
    FILE* vtkFile = open_file(fileName, "w");
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

int readParameters(const char *inputFile, lattice *grid, fluidSim *fluid,
                   REAL *delx, REAL *dely, REAL *delt, REAL *tau,
                   REAL *UI, REAL *VI, REAL *PI)
{
    FILE *input = open_file(inputFile,"r");
    if (input == NULL)
        return 0;
    char variableType[128];
    REAL value;
    int readVars = 0;
    while (fscanf(input,"%s%*[^0-9]%lg\n",variableType,&value) != EOF)
    {
        switch(variableType[0])
        {
        case 'x':
            grid->xlength = value;
            break;
        case 'y':
            grid->ylength = value;
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
            if (variableType[1] == '_') fluid->t_end = value;
            else *tau = value;
            break;
        case 'R':
            fluid->Re = value;
            break;
        case 'U':
            *UI = value;
            break;
        case 'V':
            *VI = value;
            break;
        case 'P':
            *PI = value;
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
    *delx = grid->xlength/grid->imax;
    *dely = grid->ylength/grid->jmax;
    fclose(input);
    return readVars;
}

void outputVec(REAL **U, REAL **V, REAL **P, lattice *grid, int n)
{
    if (grid == NULL)
        return;
    int imax = grid->imax;
    int jmax = grid->jmax;
    REAL **S = createMatrix(imax,jmax);
    REAL **T = createMatrix(imax,jmax);
    REAL dx = grid->xlength/imax;
    REAL dy = grid->ylength/jmax;
    char fileName[64];
    if (P != NULL)
    {
        if (n != 0) sprintf(fileName,"PressureField_%i.vtk",n);
        else sprintf(fileName,"PressureField.vtk");
        for (int i = 0;  i < imax; i++)
            for (int j = 0; j < jmax; j++)
                T[i][j] = P[i+1][j+1];
        writeVTKfileFor2DscalarField(fileName,"pressurefield",T,imax,jmax,dx,dy);
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
    writeVTKfileFor2DvectorField(fileName,"momentumfield",T,S,imax,jmax,dx,dy);
    destroyMatrix(T,imax);
    destroyMatrix(S,imax);
    return;
}
