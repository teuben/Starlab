
#include "win.h"
extern lux_wins *get_currentwin();

/* For convenience... */

#include "Xnames.h"

#define DEBUG		0

#define COLORMAP_DIR	"/home/beowulf2/starlab/2.1/src/gfx/palettes/"
#define COLORMAP_FILE	"Standard"

#define NDIR		3
static char *colormap_dir[] = {
  "./",
  "../",
  COLORMAP_DIR
};

int get_cmapfilename(char *name)
{
  /* Try various directories prepended to the input file name if
     it doesn't exist in the current directory. */

  int i;
  char tempname[256];
  FILE *cmap;

  for (i = 0; i < NDIR; i++) {
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

int read_colormap(char *filename)
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

static char colormapfile[256];
static int xo, yo, nx, ny;
static Window window;
static XEvent report;

static lux_wins *current;
static Pixmap pixmap;
static XImage *image;

static int nrep = 1;

int ximage_init(int sizex, int sizey)
{
  int i, j;

  /* Establish defaults. */

  xo = 50;
  yo = 50;
  nx = sizex;
  ny = sizey;
  strcpy(colormapfile, COLORMAP_FILE);

  /* Open an X-window. */

  if ( (window = lux_openwin(xo, yo, nx, ny)) <= 0) {

      fprintf(stderr, "Error opening X window.\n");
      exit(1);

  }

  current = get_currentwin(window);

  pixmap = XCreatePixmap(current->win.display, current->win.window, 
			 nx, ny, current->win.window_depth);
  image = XGetImage(current->win.display, pixmap, 
		    0, 0, nx, ny, AllPlanes, ZPixmap);

  /* Getting the colors right is a kludge.  The lux package apparently
     only knows how to handle 8-bit color, so use it only in that case.
     For 24-bit color, for now at least, use the default colormap
     and modify the data appropriately as it is read in... */

  get_cmapfilename(colormapfile);

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

  if (image->bits_per_pixel > 8) nrep = image->bits_per_pixel / 8;

}

int ximage(char *data)
{
  image->data = (char*)data;

  XPutImage(current->win.display, window, current->win.gc,
	    image, 0, 0, 0, 0, nx, ny);
}

int ximage_quit()
{
    lux_wins *repwin;

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

/*------------------------------------------------------------------------*/

/* Accessor functions. */

int get_nrep()
{
  return nrep;
}

int get_cmap(unsigned char *r, unsigned char *g, unsigned char *b)
{
  int i;
  for (i = 0; i < 255; i++) {
    r[i] = red[i];
    g[i] = green[i];
    b[i] = blue[i];
  }
  return 0;
}
