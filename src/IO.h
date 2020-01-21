#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <png.h>
#include "types.h"

#define NOT_PNG (0)
#define WRITE   (0x100)
#define WRITTEN (0x200)

void    print2Dfield (REAL** field, int sizeX, int sizeY);
void    write2Dfield (const char* fileName, REAL** field, int sizeX, int sizeY, const char *mode);
/* Output a matrix, either to stdout or to a file */

REAL**  read2Dfield (const char* fileName, int* sizeX, int* sizeY);
/* Reads a field from fileName */

void    writeVTKfileFor2DscalarField (const char* fileName, const char* description, REAL** field, lattice *grid);
void    writeVTKfileFor2DintegerField(const char* fileName, const char* description, char** field, lattice *grid);
void    writeVTKfileFor2DvectorField (const char* fileName, const char* description,
		REAL** fieldU, REAL** fieldV, lattice *grid);
void    WriteParticle (particle *parts, int partcount, int n);
void    outputVec (REAL **U, REAL **V, REAL **P, lattice *grid, int n);
/* Outputs the given fields as VTK file for visualisation with Paraview */

int     dumpFields(MPI_Comm Region, REAL **U, REAL **V, REAL **P, lattice *grid, int n);
void    translateBinary (MPI_Comm Region, lattice *grid, int files, int rank, int *dims);
/* Dump all matrices as binaries and translate them back to VTK */

int     check_if_png(const char *fileName, FILE **file);
void    readImageData (FILE *flagData, png_structpp png_ptr, png_infopp info_ptr);
char**  readGeometry(const char *flagFile, int *width, int *height);
/* Opens, confirms and reads a png-file into a flag array */

char**  adjustFlags(char **FLAG, int height, int width, int imax, int jmax);
void    findOptimalFlags(char **FLAG, int height, int width, int *imax, int *jmax);
/* Changes the number of geomatry cells */

int     readParameters(const char *inputFile, REAL *init,
		lattice *grid, fluidSim *sim, boundaryCond *bCond,
		REAL *delt, REAL *t_end);
/* Read parameters for a simulation from inputFile */

#endif /* IO_H_ */
