#include "geometry.h"

int check_if_png(const char *fileName, FILE **file)
{
   unsigned char buf[8];
   /* Open the alleged PNG file. */
   if ((*file = open_file(fileName, "rb")) == NULL)
      return NOT_PNG;

   /* Read the first eight bytes and compare them to the signature */
   if (fread(buf, 1, 8, *file) != 8)
      return NOT_PNG;
   return(!png_sig_cmp(buf, 0, 8));
}

void readImageData (FILE *flagData, png_structpp png_ptr, png_infopp info_ptr)
{
    if (flagData == NULL || png_ptr == NULL || info_ptr == NULL)
        return;
    /* Initialization */
    *png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    if (*png_ptr == NULL)
       return;
    *info_ptr = png_create_info_struct(*png_ptr);
    if (*info_ptr == NULL)
    {
        png_destroy_read_struct(png_ptr, NULL, NULL);
        return;
    }
    if (setjmp(png_jmpbuf(*png_ptr)))
    {
        /* Free all of the memory associated with the png_ptr and info_ptr */
        png_destroy_read_struct(png_ptr, info_ptr, NULL);
        /* If we get here, we had a problem reading the file. */
        printf("The has been an issue when reading the file. If only I knew...\n");
        return;
    }
    png_init_io(*png_ptr, flagData);

    /* Start reading data from the file (except for the signature, which has been read before) */
    png_set_sig_bytes(*png_ptr, 8);
    png_read_png(*png_ptr, *info_ptr, PNG_TRANSFORM_PACKING, NULL);
    fclose(flagData);
    return;
}

short** readGeometry (const char *flagFile, int *minimumHeight, int *minimumWidth)
{
    /* Variables */
    FILE *flagData;
    png_structp png_ptr;
    png_infop info_ptr;
    int height, width;

    if (!check_if_png(flagFile,&flagData))
    {
        printf("The given file does not seem to be a png file. Please check your input.\n");
        return NULL;
    }
    readImageData(flagData,&png_ptr,&info_ptr);

    height = png_get_image_height(png_ptr,info_ptr);
    width = png_get_image_width(png_ptr,info_ptr);

    if (height < *minimumHeight || width < *minimumWidth)
        printf("The image has height %i and width %i.\n",height,width);
    *minimumHeight = height;
    *minimumWidth = width;

    png_bytepp rows;
    rows = png_get_rows(png_ptr,info_ptr);

    short **FLAG = create2DIntegerField(height,width);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < 3*width; j+=3)
        {
            printf("(%2x,%2x,%2x), ",rows[i][j],rows[i][j+1],rows[i][j+2]);
            if ( rows[i][j] < 0x90 || rows[i][j+1] < 0x90 || rows[i][j+2] < 0x90)
                FLAG[i][j/3] = C_B;
            else
                FLAG[i][j/3] = C_F;
        }
        printf("\n");
    }

    png_destroy_read_struct(&png_ptr,&info_ptr,NULL);
    return FLAG;
}
