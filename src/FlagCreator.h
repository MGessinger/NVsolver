#ifndef FLAGCREATOR_H_
#define FLAGCREATOR_H_

#include <stdlib.h>
#include <stdio.h>
#include "IO.h"
#define C_F (0)
#define C_B (1)

void WriteFlag(const char* fileName, short** flag, int sizeX, int sizeY); //writes Flag in file (binary)
int** ReadFlag(const char* fileName, int* sizeX, int* sizeY); //reads Flag from file (binary)
int** mallocFlag(int sizeX, int sizeY); //allocates memory for flag field
void CreateFlag(); // allows user to set blocks freely
short CorrectnessCheck(short** flag, int sizeX, int sizeY);
void EditFlag(short** flag, int sizeX, int sizeY, double dx, double dy);

#endif
