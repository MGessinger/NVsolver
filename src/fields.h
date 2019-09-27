#ifndef FIELDS_H_
#define FIELDS_H_

#include "types.h"

REAL*	create1Dfield (int size);
REAL**	create2Dfield (int sizeX, int sizeY);
short** create2DIntegerField(int imax, int jmax);
/* Creates arrays of the given size */

void 	destroy1Dfield (REAL *field);
void 	destroy2Dfield (REAL **field, int sizeX);
void    destroy2DIntegerField(short **field, int imax);
/* Clears the memory allocated for field */

void	fill1Dfield (REAL value, REAL *field, int size);
void	fill2Dfield (REAL value, REAL **field, int sizeX, int sizeY);
/* Sets every entry in field to value */

int     isEqualScalar (REAL x, REAL y, REAL eps);
int 	isEqual1Dfield (REAL *field1, REAL *field2, int size, REAL eps);
int 	isEqual2Dfield (REAL **field1, REAL **field2, int sizeX, int sizeY,  REAL eps);
/* Returns 1 iff the two inputs agree up to an error of eps */

void	applyFunctionTo1Dfield (REAL (*func)(REAL), REAL *field, int size);
void	applyFunctionTo2Dfield (REAL (*func)(REAL), REAL **field, int sizeX, int sizeY);
/* Applies a given funtion to every entry of field */

REAL**  sampleFDgridOnCellCorners (REAL (*func)(REAL,REAL), lattice *grid);
REAL**  sampleFDgridOnCellCenters (REAL (*func)(REAL,REAL), lattice *grid);
/* Returns a matrix containing the discrete values of func on a lattice */

void    initUVP(REAL ***U, REAL ***V, REAL ***P, int imax, int jmax, double *init);
/* Initialise the fields U, V and P to UI, VI and PI respectively */

#endif /* FIELDS_H_ */
