/*
 *  MAKE_IMAGE:  Convert the input stream into a Sun rasterfile color image.
 *		 Assume that the numbers are *already* in the range [0,255].
 *
 *		 Use the ImageMagick "convert" program to convert the
 *		 image to other formats (e.g. gif, etc.)
 *
 *		 For a compiled version for use with C or Fortran
 *		 programs, see write_image.
 */

#include "stdinc.h"

// From make_header.C:

void make_header(int m, int n, FILE* output_file, char* colormap_file);

/* Convenient shorthand: */

#define PAD_ROW if (pad_row && (++nchar)%m_orig == 0) putchar((unsigned char)0)
#define USAGE \
"Usage:  make_image  [-c]  [-f]  [-i]  [-m color-map]  [-s sizex [sizey]]\n"

main(int argc, char** argv)
{
    /* Variables that may be set from the command line: */

    int m = 0, n = 0;		/* Image dimensions (no default) */
    int which = 2;		/* 0 ==> char, 1 ==> int, 2 ==> float [def] */
    char* colormap = NULL;	/* Use standard built-in colormap by default */

    int pad_row = 0, pad_col = 0;
    int nchar = 0;

    int m_orig;
    int i;
    float f;

    /* Parse the argument list. */

    for (i = 0; i < argc; i++)
	if (argv[i][0] == '-') {
	    switch (argv[i][1]) {
		case 'c':	which = 0;
				break;
		case 'f':	which = 2;
				break;
		case 'i':	which = 1;
				break;
		case 'm':	colormap = argv[++i];
				break;
		case 's':	m = atoi(argv[++i]);
				if (i+1 < argc && argv[i+1][0] != '-')
				    n = atoi(argv[++i]);
				else
				    n = m;
				break;

		default:	fprintf(stderr, USAGE);
				exit(1);

	    }
	}

    if (m <= 0 || n <= 0) {
	fprintf(stderr, "Must specify image dimensions.\n");
	exit(1);
    }

    /* NOTE that a peculiarity of the Sun rasterfile format is
       that the numbers of rows and columns both must be even... */

    if (m%2 != 0 || n%2 != 0) {
	m_orig = m;
	if (m%2 != 0) {
	    pad_row = 1;
	    m++;
	}
	if (n%2 != 0) {
	    pad_col = 1;
	    n++;
	}
	fprintf(stderr, 
		"Image dimensions must be even. Padding data to %dx%d\n",
		m, n);
    }

    make_header(m, n, stdout, colormap);

    if (which == 0)
	while ((i = getchar()) != EOF) {
	    putchar((unsigned char)i);
	    PAD_ROW;
	}
    else if (which == 1)
	while (scanf("%d", &i) != EOF) {
	    putchar((unsigned char)i);
	    PAD_ROW;
	}
    else
	while (scanf("%f", &f) != EOF) {
	    putchar((unsigned char)f);
	    PAD_ROW;
	}

    if (pad_col)
	for (i = 0; i < m; i++) putchar((unsigned char)0);
}
