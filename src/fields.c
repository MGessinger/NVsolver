#include "types.h"

REAL* create1Dfield (int size)
{
	if (size <= 0)
		return NULL;
	REAL *array = malloc(size*sizeof(REAL));
	fill1Dfield(0,array,size);
	return array; /* Checking for NULL is done elsewhere */
}

REAL** create2Dfield (int sizeX, int sizeY)
{
	if (sizeX <= 0 || sizeY <= 0)
		return NULL;
	REAL **matrix = malloc(sizeX*sizeof(REAL*));
	for (int i = 0; i < sizeX; i++)
	{
		matrix[i] = create1Dfield(sizeY);
		if (matrix[i] == NULL)
		{
			destroy2Dfield((void**)matrix,i);
			return NULL;
		}
	}
	fill2Dfield(0,matrix,sizeX,sizeY);
	return matrix;
}

char **create2DIntegerField (int imax, int jmax)
{
	/* Automatically add ghost cells (For exchange with other processes) */
	if (imax <= 0 || jmax <= 0)
		return NULL;
	imax += 2;
	jmax += 2;
	char **matrix = malloc(imax*sizeof(short*));
	for (int i = 0; i < imax; i++)
	{
		matrix[i] = malloc(jmax*sizeof(short));
		if (matrix[i] == NULL)
		{
			destroy2Dfield((void**)matrix,i);
			return NULL;
		}
		for (int j = 0; j < jmax; j++)
			matrix[i][j] = C_F;
	}
	return matrix;
}

void destroy1Dfield (void* field)
{
	if (field == NULL)
		return;
	free(field);
	return;
}

void destroy2Dfield (void** field, int sizeX)
{
	if (field == NULL)
		return;
	for (int i = 0; i < sizeX; i++)
	{
		destroy1Dfield((void*)field[i]);
	}
	free(field);
	return;
}

void fill1Dfield (REAL value, REAL* field, int size)
{
	if (field == NULL)
		return;
	for (int i = 0; i < size; i++)
	{
		field[i] = value;
	}
	return;
}

void fill2Dfield (REAL value, REAL** field, int sizeX, int sizeY)
{
	if (field == NULL)
		return;
	for (int i = 0; i < sizeX; i++)
		fill1Dfield(value,field[i],sizeY);
	return;
}

void initUVP (REAL ***U, REAL ***V, REAL ***P, int imax, int jmax, REAL *init)
{
	if (!U || !V || !P)
		return;
	*P = create2Dfield(imax+2,jmax+2);
	*U = create2Dfield(imax+3,jmax+2);
	*V = create2Dfield(imax+2,jmax+3);
	if (!init)
		return;
	fill2Dfield(init[2],*P,imax+2,jmax+2);
	fill2Dfield(init[0],*U,imax+3,jmax+2);
	fill2Dfield(init[1],*V,imax+2,jmax+3);
	return;
}
