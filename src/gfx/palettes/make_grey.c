
#include <stdio.h>

/*
 *	Make_grey:  Make a greyscale map (write to stdout).
 *
 *	NOTE:  0 is black, 255 is white, by (my) convention!
 */

make_colormap (unsigned char* red, unsigned char* green, unsigned char* blue)
{
    int i;

    for (i = 0; i < 256; i++) {
	red[i] = green[i] = blue[i] = i;
    }
}

main()
{
    int i;
    unsigned char red[256], green[256], blue[256];

    make_colormap(red, green, blue);

    /* Format:  256 red bytes, 256 green, 256 blue. */

    for (i = 0; i < 256; i++) putchar(red[i]);
    for (i = 0; i < 256; i++) putchar(green[i]);
    for (i = 0; i < 256; i++) putchar(blue[i]);
}
