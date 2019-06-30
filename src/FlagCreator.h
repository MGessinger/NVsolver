#ifndef FLAGCREATOR_H_
#define FLAGCREATOR_H_

#include <stdlib.h>
#include <stdio.h>

void WriteFlag(const char* fileName, int** flag, int sizeX, int sizeY); //writes Flag in file (binary)
int** ReadFlag(const char* fileName, int* sizeX, int* sizeY); //reads Flag from file (binary)
int** mallocFlag(int sizeX, int sizeY); //allocates memory for flag field
int** CreateFlag(int sizeX, int sizeY); // allows user to set blocks freely

#endif
