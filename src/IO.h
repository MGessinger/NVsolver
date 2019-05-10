#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <string.h>
#include "real.h"
#include "fields.h"

#define USER "gessinge"

FILE    *open_file (const char *fileName, const char *mode);
/* Redirects fileName to the /data/ directory on the USB stick */

void    print1Dfield (REAL* field, int size);
void    print2Dfield (REAL** field, int sizeX, int sizeY);
void    printVector (REAL* vector, int len);
void    printMatrix (REAL** matrix, int rows, int cols);
/* Prints a field to stdio */

void    write1Dfield (const char* fileName, REAL* field, int size);
void    write2Dfield (const char* fileName, REAL** field, int sizeX, int sizeY);
/* Writes a field to fileName as a binary file */

REAL*   read1Dfield (const char* fileName, int* size);
REAL**  read2Dfield (const char* fileName, int* sizeX, int* sizeY);
/* Reads a field from fileName */

void    writeVTKfileFor2DscalarField (const char* fileName, const char* description,
                                      REAL** field, int sizeX, int sizeY, REAL dx, REAL dy);
void    writeVTKfileFor2DvectorField (const char* fileName, const char* description,
                                      REAL** fieldU, REAL** fieldV, int sizeX, int sizeY, REAL dx, REAL dy);
void    outputVec (REAL **U, REAL **V, REAL **P, lattice *grid, int n);
/* Outputs the given fields as VTK file for visualisation with Paraview */

int     readParameters (const char *inputFile, lattice *grid, fluidSim *fluid,
                        REAL *delx, REAL *dely, REAL *delt, REAL *tau, REAL *UI, REAL *VI, REAL *PI, char *problem);
/* Read parameters for a simulation from inputFile */

#endif /* IO_H_ */
