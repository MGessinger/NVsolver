#ifndef FLAGCREATOR_H_
#define FLAGCREATOR_H_

#include <stdlib.h>
#include <stdio.h>
#include "IO.h"
#include "real.h"
#include <string.h> //fuer strlength und strtok
#define C_F (0)
#define C_B (1)

void WriteFlag(const char* fileName, short** flag, int sizeX, int sizeY); //writes Flag in file (binary)
short** ReadFlag(const char* fileName, int* sizeX, int* sizeY); //reads Flag from file (binary)
short **mallocFlag(int sizeX, int sizeY); //allocates memory for flag field
int CreateFlag(); // allows user to set blocks freely
short CorrectnessCheck(short** flag, int sizeX, int sizeY); //checks if a given flag field is correct
void EditFlag(short** flag, int sizeX, int sizeY, double dx, double dy); //edits an existing field
void EditFlagFromFile(); //opens a file from console to edit it
void AddRectangle(short** flag, int isize,int jsize,int ipos,int jpos,short BorF); //add a rectangle with given size to given position
void CreativeMode(short** flag, int sizeX,int sizeY, int ipos, int jpos); //CREATIVE MODE

#endif
