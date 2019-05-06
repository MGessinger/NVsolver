#ifndef BLAS_H_
#define BLAS_H_

#include <math.h>
#include <stdio.h>
#include "real.h"

void	axpy (REAL alpha, REAL *x, REAL *y, int len);
/* Computes y = alpha*x + y */

REAL	dot (REAL *x, REAL *y, int len);
/* Return <x,y> */

REAL	nrm2 (REAL *x, int len);
/* Returns the 2-norm of x via ||x|| = sqrt(<x,x>) */

void	copy (REAL *src, REAL *dest, int len);
/* Copies data from x to y */

void	scal (REAL alpha, REAL *x, int len);
/* Scales x by alpha */

void	gemv (REAL alpha, REAL **A, REAL *x, REAL beta, REAL *y, int rows, int cols);
/* Computes y = alpha*A*x + beta*y */

void	scal2Dfield (REAL alpha, REAL **X, int sizeX, int sizeY);
/* Scales the 2D-field X by alpha */

void	axpy2Dfield (REAL alpha, REAL **X, REAL **Y, int sizeX, int sizeY);
/* Computes Y = alpha*X + Y */

#endif /* BLAS_H_ */
