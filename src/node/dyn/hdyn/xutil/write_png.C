
#ifndef TOOLBOX

// Write_png - Write an image using libpng -- stripped-down version.
// This should eventually find its way to src/gfx/util.
//
// *** Requires the png and zlib libraries. ***
// *** Does not need any Starlab libraries. ***

#include "png.h"

#define ERROR	1
#define OK	0

void set_palette_colors(png_color palette[], int n)
{
    int i;
    for (i = 0; i < n; i++)
	palette[i].red = palette[i].green = palette[i].blue = i;
    // printf("Using %d-color palette.\n", n);
}

// Use separate red, green, and blue on input because existing software
// (see e.g. snap_to_image.C) uses this format.

// The image is expected to be stored starting at the bottom left.
// Conversion to "image" format (start at top left) is automatic.

// Image is unsigned char, dimensions width by height.  Currently,
// no option exists for output to stdout -- must specify a file name.

int write_png(png_bytep image, png_uint_32 width, png_uint_32 height,
	      unsigned char red[],
	      unsigned char green[],
	      unsigned char blue[],
	      char *file_name)
{
    FILE *fp;
    png_structp png_ptr;
    png_infop info_ptr;
    int bit_depth = 8;
    int bytes_per_pixel = 1;
    png_colorp palette;
    png_uint_32 k;
    png_bytep row_pointers[height];

    if (!image || width <= 0 || height <= 0) return (ERROR);
    if (!file_name) return (ERROR);
   
    /* Open the output file. */

    fp = fopen(file_name, "wb");
    if (!fp) return (ERROR);

    /* Create and initialize the png_struct. */

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
				      NULL, NULL, NULL);
    if (!png_ptr) {
	fclose(fp);
	return (ERROR);
    }

    /* Allocate and initialize the image information data. */

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
	fclose(fp);
	png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
	return (ERROR);
    }

    /* Set error handling. */

    if (setjmp(png_jmpbuf(png_ptr))) {

	/* If we get here, we had a problem reading the file. */

	fclose(fp);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	return (ERROR);
    }

    /* Set up output control using standard C streams. */

    png_init_io(png_ptr, fp);

    /* Set image information. */

   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth,
		PNG_COLOR_TYPE_PALETTE,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_BASE,
		PNG_FILTER_TYPE_BASE);

   /* Set up the palette. */

   palette = (png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH
				    		* sizeof (png_color));
   set_palette_colors(palette, PNG_MAX_PALETTE_LENGTH);
   png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);

   /* Write the file header information. */

   png_write_info(png_ptr, info_ptr);

   for (k = 0; k < height; k++)
       row_pointers[k] = image + k*width*bytes_per_pixel;

   /* Note: Image is written from the top down, so the image
    * 	    array starts at the bottom left. */

   for (k = 0; k < height/2; k++) {
       png_bytep tmp = row_pointers[k];
       row_pointers[k] = row_pointers[height-1-k];
       row_pointers[height-1-k] = tmp;
   }

   /* Write the image. */

   png_write_image(png_ptr, row_pointers);

   /* Finish writing the rest of the file. */

   png_write_end(png_ptr, info_ptr);

   /* Clean up and exit. */

   png_free(png_ptr, palette);
   palette = NULL;
   png_destroy_write_struct(&png_ptr, &info_ptr);
   fclose(fp);

   return (OK);
}

#else

#include "write_png.h"
#define N 512

main()
{
    unsigned char image[N*N];
    unsigned char red[256], green[256], blue[256];
    int i, j;

    for (j = N-1; j >= 0; j--)
	for (i = 0; i < N; i++) {
#if 0
	    if (i == j)
//	    if (i + j == N - 1)
		image[i+j*N] = i%256;
	    else
		image[i+j*N] = 0;
#else
	    image[i+j*N] = i%256;
#endif
	}

    write_png(image, N, N, red, green, blue, "xxx.png");
}

#endif
