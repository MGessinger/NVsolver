#include "fields.h"

REAL* create1Dfield(int size)
{
    if (size <= 0) return NULL;
    REAL *array = malloc(size*sizeof(REAL));
    fill1Dfield(0,array,size);
    return array; /* Checking for NULL is done elsewhere */
}

REAL** create2Dfield(int sizeX, int sizeY)
{
    if (sizeX <= 0 || sizeY <= 0) return NULL;
    REAL **matrix = malloc(sizeX*sizeof(REAL*));
    for (int i = 0; i < sizeX; i++)
    {
        matrix[i] = create1Dfield(sizeY);
        if (matrix[i] == NULL)
        {
            destroy2Dfield(matrix,i);
            return NULL;
        }
    }
    fill2Dfield(0,matrix,sizeX,sizeY);
    return matrix;
}

short** create2DIntegerField(int imax, int jmax)
{
    if (imax <= 0 || jmax <= 0)
        return NULL;
    short **matrix = malloc(imax*sizeof(short*));
    for (int i = 0; i < sizeX; i++)
    {
        matrix[i] = malloc(jmax*sizeof(short));
        if (matrix[i] == NULL)
        {
            destroy2DIntegerField(matrix,i);
            return NULL;
        }
    }
    return matrix;
}

REAL* createVector(int len)
{
    return create1Dfield(len);
}

REAL** createMatrix(int rows, int cols)
{
    return create2Dfield(rows, cols);
}

void destroy1Dfield(REAL* field)
{
    if (field == NULL)
        return;
    free(field);
    return;
}

void destroy2Dfield(REAL** field, int sizeX)
{
    if (field == NULL)
        return;
    for (int i = 0; i < sizeX; i++)
    {
        destroy1Dfield(field[i]);
    }
    free(field);
    return;
}

void destroy2DIntegerField(short **field, int imax)
{
    if (field == NULL)
        return;
    for (int i = 0; i < imax; i++)
    {
        free(field[i]);
    }
    free(field);
    return;
}

void destroyVector(REAL* vector)
{
    destroy1Dfield(vector);
}

void destroyMatrix(REAL** matrix, int rows)
{
    destroy2Dfield(matrix,rows);
}

void fill1Dfield(REAL value, REAL* field, int size)
{
    if (field == NULL)
        return;
    for (int i = 0; i < size; i++)
    {
        field[i] = value;
    }
    return;
}

void fill2Dfield(REAL value, REAL** field, int sizeX, int sizeY)
{
    if (field == NULL)
        return;
    for (int i = 0; i < sizeX; i++)
    {
        fill1Dfield(value,field[i],sizeY);
    }
    return;
}

int isEqualScalar(REAL x, REAL y, REAL eps)
{
    /* In any case, one of these inequalities is trivial, but both have to hold */
    if (x-y <= eps && y-x <= eps)
        return 1;
    return 0;
}

int isEqual1Dfield(REAL* field1, REAL* field2, int size, REAL eps)
{
    if (field1 == NULL || field2 == NULL)
        return 0;
    for (int i = 0; i < size; i++)
    {
        if (isEqualScalar(field1[i],field2[i],eps) == 0)
            return 0;
    }
    return 1;
}

int isEqual2Dfield(REAL** field1, REAL** field2, int sizeX, int sizeY, REAL eps)
{
    if (field1 == NULL || field2 == NULL)
        return 0;
    for (int i = 0; i < sizeX; i++)
    {
        if (isEqual1Dfield(field1[i],field2[i],sizeY,eps) == 0)
            return 0;
    }
    return 1;
}

void applyFunctionTo1Dfield(REAL (*func)(REAL), REAL* field, int size)
{
    if (field == NULL)
        return;
    if (func == NULL)
        return;
    for (int i = 0; i < size; i++)
        field[i] = func(field[i]);
    return;
}

void	applyFunctionTo2Dfield(REAL (*func)(REAL), REAL** field, int sizeX, int sizeY)
{
    if (field == NULL)
        return;
    if (func == NULL)
        return;
    for (int i = 0; i < sizeX; i++)
        applyFunctionTo1Dfield(func,field[i],sizeY);
    return;
}

REAL** sampleFDgridOnCellCorners(REAL (*func)(REAL,REAL), REAL ilength, REAL jlength,int imax, int jmax)
{
    REAL dx = ilength/(imax+1);
    REAL dy = jlength/(jmax+1);
    REAL **A = createMatrix(imax,jmax);
    if (A == NULL)
        return NULL;
    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
            A[i][j] = func(dx*(i+1),dy*(j+1));
    return A;
}

REAL** sampleFDgridOnCellCenters(REAL (*func)(REAL,REAL), REAL ilength, REAL jlength,int imax, int jmax)
{
    REAL dx = ilength/imax;
    REAL dy = jlength/jmax;
    REAL **A = createMatrix(imax,jmax);
    if (A == NULL)
        return NULL;
    for (int i = 0; i < imax; i++)
        for (int j = 0; j < jmax; j++)
            A[i][j] = func(dx*(i+1./2),dy*(j+1./2));
    return A;
}

void initUVP(REAL **U, REAL **V, REAL **P, int imax, int jmax, REAL UI, REAL VI, REAL PI)
{
    if (U == NULL || V == NULL || P == NULL)
        return;
    fill2Dfield(UI,U,imax+2,jmax+2);
    fill2Dfield(VI,V,imax+2,jmax+2);
    fill2Dfield(PI,P,imax+2,jmax+2);
    return;
}
