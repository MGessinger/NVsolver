#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <png.h>
#include <zlib.h>
#include "boundary.h"
#include "IO.h"
#include "real.h"

#define NOT_PNG (0)

int check_if_png(const char *fileName, FILE **file);

void readImageData (FILE *flagData, png_structpp png_ptr, png_infopp info_ptr);

void readGeometry (const char *flagFile, int minimumHeight, int minimumWidth);

#endif /* GEOMETRY_H_ */
