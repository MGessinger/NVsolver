#include "blas.h"

void axpy(REAL alpha, REAL* x, REAL* y, int len)
{
    if (x == NULL || y == NULL)
    {
        printf("A field is invalid. Please check your input.\n");
        return;
    }
    for (int i = 0; i < len; i++)
    {
        y[i] = alpha*x[i]+y[i];
    }
    return;
}

REAL dot(REAL* x, REAL* y, int len)
{
    if (x == NULL || y == NULL)
    {
        printf("A field is invalid. Please check your input.\n");
        return 0;
    }
    REAL alpha = 0;
    for (int i = 0; i < len; i++)
    {
        alpha += x[i]*y[i];
    }
    return alpha;
}

REAL nrm2(REAL* x, int len)
{
    return sqrt(dot(x,x,len));
}

void copy(REAL* x, REAL* y, int len)
{
    if (x == NULL || y == NULL)
    {
        printf("A field is invalid. Please check your input.\n");
        return;
    }
    for (int i = 0; i < len; i++)
    {
        y[i] = x[i];
    }
    return;
}

void scal(REAL alpha, REAL* x, int len)
{
    if (x == NULL)
    {
        printf("The field is invalid. Please check your input.\n");
        return;
    }
    for (int i = 0; i < len; i++)
    {
        x[i] *= alpha;
    }
    return;
}

void gemv(REAL alpha, REAL** A, REAL* x, REAL beta, REAL* y, int rows, int cols)
{
    if (x == NULL || y == NULL || A == NULL)
    {
        printf("A field is invalid. Please chech your input.\n");
        return;
    }
    for (int i = 0; i < rows; i++)
    {
        y[i] *= beta;
        if (alpha != 0) y[i] += alpha*dot(A[i],x,cols);
    }
    return;
}

void scal2Dfield(REAL alpha, REAL** X, int sizeX, int sizeY)
{
    if (X == NULL)
    {
        printf("The input field ins invalid.\n");
        return;
    }
    for (int i = 0; i < sizeX; i++)
    {
        scal(alpha,X[i],sizeY);
    }
    return;
}

void axpy2Dfield(REAL alpha, REAL** X, REAL** Y, int sizeX, int sizeY)
{
    if (X == NULL || Y == NULL)
    {
        printf("A field is invalid. Please check your input.\n");
        return;
    }
    for (int i = 0; i < sizeX; i++)
    {
        axpy(alpha,X[i],Y[i],sizeY);
    }
    return;
}
