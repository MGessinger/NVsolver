#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <png.h>
#include "types.h"

#define NOT_PNG (0)

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

void    writeVTKfileFor2DscalarField (const char* fileName, const char* description, REAL** field, lattice *grid);
void    writeVTKfileFor2DintegerField(const char* fileName, const char* description, short** field, lattice *grid);
void    writeVTKfileFor2DvectorField (const char* fileName, const char* description,
                                      REAL** fieldU, REAL** fieldV, lattice *grid);
void    WriteParticle (particle *parts, int partcount, int n);
void    outputVec (REAL **U, REAL **V, REAL **P, lattice *grid, int n);
/* Outputs the given fields as VTK file for visualisation with Paraview */

int     check_if_png(const char *fileName, FILE **file);
void    readImageData (FILE *flagData, png_structpp png_ptr, png_infopp info_ptr);
short** readGeometry(const char *flagFile, int *width, int *height);
/* Opens, confirms and reads a png-file into a flag array */

short** adjustFlags(short **FLAG, int height, int width, int imax, int jmax);
void    findOptimalFlags(short **FLAG, int height, int width, int *imax, int *jmax);
/* Changes the number of geomatry cells */

int     readParameters(const char *inputFile, REAL *init,
                       lattice *grid, fluidSim *sim, boundaryCond *bCond,
                       REAL *delt, REAL *t_end);
/* Read parameters for a simulation from inputFile */

#endif /* IO_H_ */
