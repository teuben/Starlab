
#include <stdio.h>

/*
 *	random_map:  Make a random color map (write to stdout).
 *
 *	NOTE:  0 is black, 255 is white, by (my) convention!
 */

make_colormap (unsigned char* red, unsigned char* green, unsigned char* blue)
{
    int i;

    red[0] = green[0] = blue[0] = 0;
    for (i = 1; i < 255; i++) {
      red[i] = random() % 256;
      green[i] = random() % 256;
      blue[i] = random() % 256;
      while (red[i] + green[i] + blue[i] < 128) {
	red[i] *= 2;
	green[i] *= 2;
	blue[i] *= 2;
      }
    }
    red[255] = green[255] = blue[255] = 255;
}

main(int argc, char *argv[])
{
    int i;
    unsigned char red[256], green[256], blue[256];

    if (argc <= 1)
      srand(42);
    else
      srand(atoi(argv[1]));

    make_colormap(red, green, blue);

    /* Format:  256 red bytes, 256 green, 256 blue. */

    for (i = 0; i < 256; i++) putchar(red[i]);
    for (i = 0; i < 256; i++) putchar(green[i]);
    for (i = 0; i < 256; i++) putchar(blue[i]);
}
