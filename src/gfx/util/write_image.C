
//  WRITE_IMAGE:  Convert an array of real values a[m, n] into a color image.
//		  Assume that the numbers are already in the range [0,1].
//
//		  This version is called from Fortran or C.
//		  For a UNIX pipe version, use make_image.

#include "stdinc.h"

// From make_header.C:

void make_header(int m, int n, FILE* out_file,
		 char* colormap_file);
void make_header(int m, int n, FILE* out_file,
		 unsigned char *red,
		 unsigned char *green,
		 unsigned char *blue);

// The description and construction of the header for SUN raster
// images is now in make_header.c

#define BUFSIZE 4096

static void write_line(float* a, int m, FILE* file,
		       int scale, float amin, float amax)
{
    int i, modd = 0;
    unsigned char c[BUFSIZE];
    float afac, as[BUFSIZ];

    int cmin, cmax;

    if (m%2 != 0) modd = 1;

    // Options: no scaling, expect input data in range [0, 1)
    // 		with scaling, will scale min to 0, max to 255

    if (!scale || amax <= amin)
	for (i = 0; i < m; i++) as[i] = 255.0*a[i];	  // will vectorize
    else {
	afac = 255.0/(amax-amin);
	for (i = 0; i < m; i++) as[i] = afac*(a[i]-amin); // will vectorize
    }

    for (i = 0; i < m; i++)
	c[i] = (unsigned char) as[i]; // Won't vectorize

    fwrite(&c, 1, m, file);			// may vectorize?
    if (modd) fwrite(&c, 1, 1, file);		// pad to 16-bit boundary.
}

static void get_limits(float* a, int n, float* amin, float* amax)
{
    int i;

    *amin = 1.0e30;
    *amax = -1.e30;
    for (i = 0; i < n; i++) {
	if (a[i] < *amin) *amin = a[i];
	if (a[i] > *amax) *amax = a[i];
    }
}

// These "write_image" functions are identical except for the arguments
// and the calls to make_header.  (Sorry...)

void write_image(float* a, int m, int n, char* filename, int scale,
		 char *colormap)

// Write the contents of the 2-D array a to a file as a SUN raster image.
// If no file is specified, use stdout.
// Scale the data before writing if scale is set.

{
    FILE* file;
    int file_open = 0;
    int j, mm, nn;
    float amin, amax;

    if (m > BUFSIZE) exit(1);

    if (filename != NULL) {
	if ((file = fopen(filename, "w")) == NULL) exit(1);
    } else
	file = stdout;

    // Make dimensions even (pad if necessary):

    mm = m;
    if (mm%2 != 0) mm++;

    nn = n;
    if (nn%2 != 0) nn++;

    make_header(mm, nn, file, colormap);

    // Scan 2-D array top to bottom, left to right.

//    if (scale) {
	get_limits(a, m*n, &amin, &amax);
//	fprintf(stderr, "%s: min = %f max = %f\n", filename, amin, amax);
//    }

    // Note: j ordering goes from top to bottom...

    for (j = n - 1; j >= 0; j--) write_line(a+m*j, m, file,
					    scale, amin, amax);

    if (nn != n) write_line(a, m, file,		// Pad to 16-bit boundary.
			    scale, amin, amax);
	
    if (filename != NULL) fclose(file);
}

void write_image(float* a, int m, int n, char* filename, int scale,
		 unsigned char *red,
		 unsigned char *green,
		 unsigned char *blue)

// Write the contents of the 2-D array a to a file as a SUN raster image.
// If no file is specified, use stdout.
// Scale the data before writing if scale is set.

{
    FILE* file;
    int file_open = 0;
    int j, mm, nn;
    float amin, amax;

    if (m > BUFSIZE) exit(1);

    if (filename != NULL) {
	if ((file = fopen(filename, "w")) == NULL) exit(1);
    } else
	file = stdout;

    // Make dimensions even (pad if necessary):

    mm = m;
    if (mm%2 != 0) mm++;

    nn = n;
    if (nn%2 != 0) nn++;

    make_header(mm, nn, file, red, green, blue);

    // Scan 2-D array top to bottom, left to right.

//    if (scale) {
	get_limits(a, m*n, &amin, &amax);
//	fprintf(stderr, "%s: min = %f max = %f\n", filename, amin, amax);
//    }

    // Note: j ordering goes from top to bottom...

    for (j = n - 1; j >= 0; j--) write_line(a+m*j, m, file,
					    scale, amin, amax);

    if (nn != n) write_line(a, m, file,		// Pad to 16-bit boundary.
			    scale, amin, amax);
	
    if (filename != NULL) fclose(file);
}


//**********************************************************************
//
//		Fortran interfaces to the C++ routines:
//
//**********************************************************************

// Need to take care of name mangling...  (Not done yet.)

void write_image_f77(float* a, int* m, int* n, char* filename,
		     int* scale, long int nl)
{
    int j, file_open = 0;
    char* f = NULL;

    for (j = 0; j < nl; j++)
	if (*(filename+j) > ' ') file_open = 1;

    if (file_open) {

	// Make a C string for the filename.

	if ((f = (char*)malloc(1+nl)) == NULL) exit(1);

	// Make a new string, since filename isn't null terminated...

	for (j = 0; j < nl; j++) *(f+j) = *(filename+j);
	*(f+nl) = '\0';
    }

    // Use the C routine to do the work.

    write_image(a, *m, *n, f, *scale, NULL);
}

// Some brain-dead Fortrans (e.g. on the Sun and Cray) want global
// names to be uppercase or terminated with "_"...

void WRITE_IMAGE_F77(float* a, int* m, int* n, char* filename,
		 int* scale, long int nl)
{
    write_image_f77(a, m, n, filename, scale, nl);
}

void write_image_f77_(float* a, int* m, int* n, char* filename,
		 int* scale, long int nl)
{
    write_image_f77(a, m, n, filename, scale, nl);
}
