
#include "stdinc.h"

void make_standard_colormap(unsigned char* red,
			    unsigned char* green,
			    unsigned char* blue)
{
    // Make a simple standard colormap.

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

inline void swap(unsigned char *a, unsigned char *b)
{
    unsigned char tmp = *a;
    *a = *b;
    *b = tmp;
}

void make_alternate_colormap(unsigned char* red,
			     unsigned char* green,
			     unsigned char* blue)
{
    // Make a "reverse standard colormap" (red --> blue).

    int i;

    for (i = 0; i < 48; i++) {				// violet to blue
	red[i] = 5*(48-i)-1;				// 239 119 255
	green[i] = red[i]/2;				//   4   2 255
	blue[i] = 255;
    }
    for (i = 48; i < 86; i++) {				// blue to blue-green
        red[i] = 0;					//   0   0 255
        green[i] = (unsigned char)(6.84*(i-48));	//   0 251 137
        blue[i] = 255 - (unsigned char)(2.5*(i-48));
    }
    for (i = 86; i < 100; i++) {			// blue-green to green
        red[i] = 0;					//   0 251 180
        green[i] = 255 - (100-i)/3;			//   0 255   5
        blue[i] = 135 - 10*(i-86);
    }
    for (i = 100; i < 130; i++) {			// green to yellow
        red[i] = (unsigned char)(8.5*(i-100));		//   0 255   0
        green[i] = 255;					// 246 p255   0
        blue[i] = 0;
    }
    for (i = 130; i < 215; i++) {			// yellow to red
        red[i] = 255;					// 255 255   0
        green[i] = 255 - 3*(i-130);			// 255   3   0
        blue[i] = 0;
    }
    for (i = 215; i < 256; i++) {			// red to dark red
        red[i] = 255 - 2*(i-215);			// 255 255   0
        green[i] = 0;					// 255 175   0
        blue[i] = 0;
    }

    for (i = 0; i < 128; i++) {
	swap(red+i, red+255-i);
	swap(green+i, green+255-i);
	swap(blue+i, blue+255-i);
    }
}

void make_greymap(unsigned char* red,
		  unsigned char* green,
		  unsigned char* blue)
{
    // Make a simple greyscale colormap.

    for (int i = 0; i < 256; i++)
        red[i] = green[i] = blue[i] = i;
}

void get_colormap(unsigned char *red,
		  unsigned char *green,
		  unsigned char *blue,
		  char *colormap_file)		// default = NULL
{
    // Construct a (SUN-style) color map.

    FILE* map = NULL;

    // Attempt to read a color map.

    if (colormap_file) {

	if ((map = fopen(colormap_file, "r")) == NULL)
	    cerr << "Can't open color map file " << colormap_file
		 << ": using default" << endl;

	else {
	    if (fread(red, 1, 256, map) != 256) {
		cerr << "Error reading color map: using default" << endl;
		colormap_file = NULL;
	    } else {
		if (fread(green, 1, 256, map) != 256) {
		    cerr << "Error reading color map: using default" << endl;;
		    colormap_file = NULL;
		} else {
		    if (fread(blue, 1, 256, map) != 256) {
			cerr << "Error reading color map: using default"
			     << endl;
			colormap_file = NULL;
		    }
		}
	    }
	}
    }

    if (map) fclose(map);

    if (!colormap_file) {

	// Default color map:

	make_standard_colormap(red, green, blue);
        // make_greymap(red, green, blue);
	// make_alternate_colormap(red, green, blue);
    }

}
