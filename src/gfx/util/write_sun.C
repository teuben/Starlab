
#include "stdinc.h"
#include "image_fmt.h"

/* Convenient shorthand: */

#define PAD_ROW if (pad_row && (i+1)%m_orig == 0) putchar((unsigned char)0)

int write_sun(FILE *dst, int m, int n,
	      unsigned char *pixels, char *colormap_file)
{
    // Note that a peculiarity of the Sun rasterfile format is
    // that the numbers of rows and columns both must be even...

    int m_orig = 1;
    bool pad_row = false, pad_col = false;

    if (m%2 != 0 || n%2 != 0) {
	m_orig = m;
	if (m%2 != 0) {
	    pad_row = true;
	    m++;
	}
	if (n%2 != 0) {
	    pad_col = true;
	    n++;
	}

	cerr << "Image dimensions must be even. Padding data to "
	     << m << " x " << n << endl;
    }

    write_sun_header(m, n, dst, colormap_file);

    for (int i = 0; i < m*n; i++) {
	fputc(pixels[i], dst);
	PAD_ROW;
    }

    if (pad_col)
	for (int i = 0; i < m; i++) fputc((unsigned char)0, dst);

    return OK;
}

// This differs from the above only in the arguments and sun_header call...

int write_sun(FILE *dst, int m, int n,
	      unsigned char *pixels,
	      unsigned char *red,
	      unsigned char *green,
	      unsigned char *blue)
{
    // Note that a peculiarity of the Sun rasterfile format is
    // that the numbers of rows and columns both must be even...

    int m_orig = 1;
    bool pad_row = false, pad_col = false;

    if (m%2 != 0 || n%2 != 0) {
	m_orig = m;
	if (m%2 != 0) {
	    pad_row = true;
	    m++;
	}
	if (n%2 != 0) {
	    pad_col = true;
	    n++;
	}

	cerr << "Image dimensions must be even. Padding data to "
	     << m << " x " << n << endl;
    }

    write_sun_header(m, n, dst, red, green, blue);

    for (int i = 0; i < m*n; i++) {
	fputc(pixels[i], dst);
	PAD_ROW;
    }

    if (pad_col)
	for (int i = 0; i < m; i++) fputc((unsigned char)0, dst);

    return OK;
}
