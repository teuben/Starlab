
/* ximage.c: Simple program to plot simple image data...
 *
 *	arguments:	nx, ny		- specify image dimensions
 *			-s nx [ny]	- same (ny = nx if unspecified)
 *			-c colormap	- specify colormap file
 */

#define USAGE "[nx] [ny] [-c colormapfile] [-o xo [yo]] [-s nx [ny]]"

#include "win.h"
extern lux_wins *get_currentwin();

/* For convenience... */

#include "Xnames.h"

#define DEBUG		0

#define NX_DEFAULT	128
#define NY_DEFAULT	128
#define COLORMAP_DIR	"/src/gfx/palettes/"
#define COLORMAP_FILE	"Standard"

#define NDIR		3
static char *colormap_dir[] = {
  "./",
  "../",
  COLORMAP_DIR
};

void get_cmapfilename(char *name)
{
  /* Try various directories prepended to the input file name if
     it doesn't exist in the current directory.  Change name to
     the correct name, if found. */

  int i;
  char tempname[256];
  FILE *cmap;

  for (i = 0; i < NDIR; i++) {

    if (i == 2) {
      char *local_directory;
      if ((local_directory = getenv("STARLAB_PATH")) == NULL) {
	name[0] = '\0';
	return;
      }
      strcpy(tempname, local_directory);
      strcat(tempname, colormap_dir[i]);
    } else
      strcpy(tempname, colormap_dir[i]);

    strcat(tempname, name);
    if ((cmap = fopen(tempname, "r"))) {
      fclose(cmap);
      strcpy(name, tempname);
      return;
    }
  }

  name[0] = '\0';

}

static unsigned char red[256], green[256], blue[256];

void read_colormap(char *filename)
{
  FILE *cmap;

  if ((cmap = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Unable to open colormap file.\n");
    return;
  }

  fread(red,   1, 256, cmap);
  fread(green, 1, 256, cmap);
  fread(blue,  1, 256, cmap);

  fclose(cmap);
}

int get_data(unsigned char *data, int n, int nr)
{
  /* Read the next chunk of data (the next image). 
     Only nr = 1 and nr = 4 are actually supported... */

  int i = 0, c;

  for (i = 0; i < n; i++) {

    if ((c = getchar()) == EOF) break;

    if (nr == 1)

      *(data++) = c;

    else {

      /* Turn the input value into a 3-byte (B, G, R) sequence,
	 followed by a null byte. */

      *(data++) = blue[c];
      *(data++) = green[c];
      *(data++) = red[c];
      *(data++) = 0;

    }

  }
  return i*nr;
}

main(int argc, char* argv[])
{
  int i, j;
  char colormapfile[256];
  char *dir;
  int xo, yo, nx, ny;
  Window window;
  XEvent report;

  /* Establish defaults. */

  xo = 50;
  yo = 50;
  nx = NX_DEFAULT;
  ny = NY_DEFAULT;
  strcpy(colormapfile, COLORMAP_FILE);
  j = 0;

  /* Parse the argument list. */

  for (i = 1; i < argc; i++) {

    if (argv[i][0] == '-') {

      switch (argv[i][1]) {

      case 'c':	strcpy(colormapfile, argv[++i]);
		break;

      case 'o':	xo = atoi(argv[++i]);
		if (i < argc-1 && argv[i+1][0] != '-')
		  yo = atoi(argv[++i]);
		else
		  yo = xo;
		break;

      case 's':	nx = atoi(argv[++i]);
		if (i < argc-1 && argv[i+1][0] != '-')
		  ny = atoi(argv[++i]);
		else
		  ny = nx;

		break;

      case 'h':
      case '-':
      default:	dir = rindex(argv[0], '/');
	        if (dir == NULL) dir = argv[0] - 1;
		fprintf(stderr, "%s: %s\n", dir+1, USAGE);
		exit(0);

      }

    } else {

      /* Interpret argument as nx or ny. */

      if (j == 0)
	nx = atoi(argv[i]);
      else
	ny = atoi(argv[i]);
      j++;

    }
  }

  /* Open an X-window. */

  if ( (window = lux_openwin(xo, yo, nx, ny)) <= 0) {

      fprintf(stderr, "Error opening X window.\n");
      exit(1);

  } else {

    int ns = 0, nbytes, nrep;
    lux_wins *current, *repwin;
    Pixmap pixmap;
    XImage *image;
    unsigned char *data;

    current = get_currentwin(window);

    pixmap = XCreatePixmap(current->win.display, current->win.window, 
			   nx, ny, current->win.window_depth);
    image = XGetImage(current->win.display, pixmap, 
		      0, 0, nx, ny, AllPlanes, ZPixmap);

    fprintf(stderr, "Image width = %d\n", image->width);
    fprintf(stderr, "Image height = %d\n", image->height);
    fprintf(stderr, "Image depth = %d\n", image->depth);
    fprintf(stderr, "Bits per pixel = %d\n", image->bits_per_pixel);

    /* Getting the colors right is a kludge.  The lux package apparently
       only knows how to handle 8-bit color, so use it only in that case.
       For 24-bit color, for now at least, use the default colormap
       and modify the data appropriately as it is read in... */

    get_cmapfilename(colormapfile);	/* Note: changes colormapfile */

    if (colormapfile[0] == '\0')

      fprintf(stderr, "Can't load color map file.\n");

    else {

      fprintf(stderr, "Color map file is %s\n", colormapfile);

      if (image->depth == 8) {

	/* Load in a standard 8-bit colormap. */

	lux_set_window_colormap(window, colormapfile);

      } else {

	/* Read the colormap file for use in manipulating the data. */

	read_colormap(colormapfile);
      }

    }

    /* Create the data array.  Note that 24-bit applications use
       4 bytes per pixel... */

    nbytes = nx*ny;
    nrep = 1;
    if (image->bits_per_pixel > 8) {
      nrep = image->bits_per_pixel / 8;
      nbytes *= nrep;
    }

    if ((data = (unsigned char*)malloc(nbytes)) == NULL) {
      fprintf(stderr, "Error allocating memory for image.\n");
      exit(1);
    }

    while (get_data(data, nx*ny, nrep) == nbytes) {

      image->data = (char*)data;

      XPutImage(current->win.display, window, current->win.gc,
		image, 0, 0, 0, 0, nx, ny);

      /* Sleep for ns microseconds between frames. */

      lux_pause(ns);

    }

    /* Enter idle mode before quitting. */
    
    fprintf(stderr,
	  "Press any key or button in display window to quit current plot\n");

    while (1) {

      /* Get the next event from the queue. */

      XNextEvent(current->win.display, &report);

      repwin = get_currentwin(report.xany.window);

      if (repwin->win.window == window) {

	if (DEBUG)
	  fprintf(stderr, "event %d (%s)\n",
		  report.type, event_name[report.type]);

	if (report.type == Expose || report.type == ConfigureNotify)

	  XPutImage(current->win.display, window, current->win.gc,
		    image, 0, 0, 0, 0, nx, ny);

	else if (report.type == KeyPress || report.type == ButtonPress)

	  break;

      }
    }
  }
}
