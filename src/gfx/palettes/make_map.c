
#include <stdio.h>

/*
 *	Make_map:  Make a color map (write to stdout).
 *
 *	NOTE:  0 is black, 255 is white, by (my) convention!
 */

make_colormap (unsigned char* red, unsigned char* green, unsigned char* blue)
{
    int i;

    for (i = 0; i < 32; i++) {
	red[i] = 0;
	green[i] = 0;
	blue[i] = 8*i;
    }
    for (i = 32; i < 96; i++) {
	red[i] = 0;
	green[i] = 4*(i-32);
	blue[i] = 255;
    }
    for (i = 96; i < 160; i++) {
	red[i] = 4*(i-96);
	green[i] = 255;
	blue[i] = 255 - red[i];
    }
    for (i = 160; i < 224; i++) {
	red[i] = 255;
	green[i] = 255 - 4*(i-160);
	blue[i] = 0;
    }
    for (i = 224; i < 256; i++) {
	red[i] = 255;
	blue[i] = green[i] = 8*(i-224);
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
