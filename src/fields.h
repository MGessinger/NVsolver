#ifndef FIELDS_H_
#define FIELDS_H_

#include "types.h"

REAL*	create1Dfield (int size);
REAL**	create2Dfield (int sizeX, int sizeY);
char** create2DIntegerField (int imax, int jmax);
/* Creates arrays of the given size */

void 	destroy1Dfield (void *field);
void 	destroy2Dfield (void **field, int sizeX);
/* Clears the memory allocated for field */

void	fill1Dfield (REAL value, REAL *field, int size);
void	fill2Dfield (REAL value, REAL **field, int sizeX, int sizeY);
/* Sets every entry in field to value */

void    initUVP (REAL ***U, REAL ***V, REAL ***P, int imax, int jmax, REAL *init);
/* Initialise the fields U, V and P to UI, VI and PI respectively */

#endif /* FIELDS_H_ */
