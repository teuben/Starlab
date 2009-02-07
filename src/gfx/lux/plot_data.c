/*
 *  plot_data:	create a simple plot of the input data stream.
 *		Type plot_data --help to see current options.
 *
 *		This version plots multiple columns and allows
 *		inline commands to be embedded in the data.
 *
 *  Author:	Copyright (c) Steve McMillan, Drexel University, 1994-2000.
 */


/* Wish list:

   1. Tolerance of errors in data
   2. log plots
   3. dialogs, etc.

*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define NMAX	 50000
#define NYMAX	 10
#define NBUFFER  100
#define INFINITY 1.e15	/* Close enough! */

#define USAGE 	 "\n\
Usage:  plot_data  [-c xcol ycol1 ycol2 ...]  [-c[c][xy]]  [-C color1...] \n\
		   [-e]  [-h header]  [-i]  [-l xmin xmax ymin ymax] \n\
		   [-L scale]  [-N nmax]  [-o]  [-O xo yo]  [-p[p]] \n\
 		   [-P point_size]  [-q]  [-s xs ys]  [-S skip]  [-t ntrail] \n\
		   [-w[w][xy]]  [-W]  [-x x-label]  [-X]  [-y y-label] \n\
		   [-z zcol] [[-]-help] \n"

typedef unsigned long Window;

typedef struct {float xmin;
		float xmax;
		float ymin;
		float ymax;
		int xmin_from_data;
		int xmax_from_data;
		int ymin_from_data;
		int ymax_from_data;
	    } limits;

typedef struct {Window xwin;		/* X window ID		      */

		int xcol;		/* X-column to graph	      */
		int nmax;		/* Length of storage arrays   */
		int ycol[NYMAX];	/* Y-column(s) to graph.      */
		int ny;			/* Number of y columns	      */
		int zcol;		/* Colors from this column    */

		int lines;		/* Plot lines	 	      */
		float point_size;	/* Point size (0 ==> 1 pixel) */
		int ntrail;		/* >0 ==> trail length	      */
		float* xsave;		/* For saving trail data      */
		float* ysave;
		int output;		/* Don't pass data to stdout  */

		char xlabel[256];	/* X-axis label		      */
		int  xlabel_set;
		char ylabel[256];	/* Y-axis label 	      */
		int  ylabel_set;
		char header[256];	/* Overall label 	      */

		char color[NYMAX][16];	/* Point color(s) (by name)   */
		int icolor[NYMAX];	/* X indices of colors	      */
		int ibox;		/* Color index of box	      */
		int iback;		/* Color index of background  */
		int curr_color;		/* Current color index	      */

		int limits_from_data;	/* Get limits from input data */
		limits lims;		/* X and Y plotting limits    */

		int crop_data;		/* Crop x or y to plot limits */
		int crop_x;		/* Crop x to plot limits      */
		int crop_y;		/* Crop y to plot limits      */
		int wrap_data;		/* Wrap x or y to plot limits */
		int wrap_x;		/* Wrap x to plot limits      */
		int wrap_y;		/* Wrap y to plot limits      */
		float yscale;		/* "Aspect ratio"	      */

		int xorigin;		/* Top left of plot window    */
		int yorigin;		/* Top left of plot window    */
		int xsize;		/* Size of plot window	      */
		int ysize;		/* Size of plot window	      */

		int quiet;		/* Suppress most output	      */
		int ignore_inline;	/* Ignore inline commands     */
		int skip;		/* Skip leading lines	      */

		int modify;		/* Internal flag	      */
		int pause;		/* Internal flag	      */

		int echo;

	    } plot_params;

static int init_arr_size, arr_size;	/* Ugly but simple! */
static int in_line = 0;

void err_exit(char* s)
{
    fprintf(stderr, "Error: %s\n", s);
    if (!in_line) exit(1);
}

static int counter = 0;
static char last_err[1024];

void print_help(char* s)
{
    if (!s) counter = 0;
    if (counter == 1 && strcmp(s, last_err)) counter = 0;

    /* Print error message only if s is not the same as last time around. */

    if (counter == 0) {

	if (s) fprintf(stderr, "Unknown option %s\n\n", s);
	fprintf(stderr, USAGE);

	fprintf(stderr, "\nOptions:\n");
	fprintf(stderr,
"\t-c xcol ycol1...  plot data in xcol horizontally, ycol vertically [1 2]\n");
	fprintf(stderr,
		"\t-c[c][xy]\tcrop [x or y] data to plot limits\n");
	fprintf(stderr,
		"\t-C color1...\tspecify line/point colors [all black]\n");
	fprintf(stderr,
		"\t-e\t\techo current settings\n");
	fprintf(stderr,
		"\t-h header\tspecify overall label for plot [none]\n");
	fprintf(stderr,
		"\t-i\t\tignore inline commands [don't ignore]\n");
	fprintf(stderr,
		"\t-l xmin xmax ymin ymax\tspecify limits for plot [get from data]\n");
	fprintf(stderr,
		"\t-ll\t\tforce plot lines only [plot lines]\n");
	fprintf(stderr,
		"\t-L scale\tspecify limits to be +/- scale for both axes\n");
	fprintf(stderr,
		"\t-N nmax\t\tspecify maximum number of points to store [%d]\n",
		NMAX);
	fprintf(stderr,
		"\t-o\t\techo stdin to stdout [do not echo]\n");
	fprintf(stderr,
		"\t-O xo yo\tspecify top left corner of plotting box [150, 50]\n");
	fprintf(stderr,
		"\t-p\t\ttoggle plot points only [plot lines]\n");
	fprintf(stderr,
		"\t-pp\t\tforce plot points only [plot lines]\n");
	fprintf(stderr,
		"\t-P size\t\tspecify point size, in x-axis units [0 ==> pixel]\n");
	fprintf(stderr,
		"\t-q\t\tsuppress most output [don't suppress]\n");
	fprintf(stderr,
		"\t-Q\t\tquit\n");
	fprintf(stderr,
		"\t-s xs ys\tspecify box size [500, 500]\n");
	fprintf(stderr,
		"\t-S skip\t\tskip leading lines [0]\n");
	fprintf(stderr,
		"\t-t ntrail\tspecify number of trailing points [infinite]\n");
	fprintf(stderr,
		"\t-w[w][xy]\twrap [x or y] data to plot limits\n");
	fprintf(stderr,
		"\t-W\t\twait for keyboard input [inline only]\n");
	fprintf(stderr,
		"\t-x xlabel\tspecify label for x-axis [\"column 'xcol'\"]\n");
	fprintf(stderr,
		"\t-X\t\tclear the display and redraw axes [inline only]\n");
	fprintf(stderr,
		"\t-y ylabel\tspecify label for y-axis [\"column 'ycol1...'\"]\n");
	fprintf(stderr,
		"\t-z zcol\t\tspecify column for color data [none]\n");
	fprintf(stderr,
		"\t[-]-help\t\tprint this help message\n");
    }

    counter++;
    if (s) strcpy(last_err, s);
    
}

void swap (float* x, float* y)
{
    float z;

    z = *x;
    *x = *y;
    *y = z;
}

void set_or_flag(char* arg, int* flag, float* limit)
{
    if (!strcmp(arg, "."))
	*flag = 1;
    else {
	*flag = 0;
	*limit = atof(arg);
    }
}

char *get_version(char *cvs_id)
{
    // Extract a version string from the CVS id.
    // CVS id format is "$Revision$".

    if (!cvs_id) return NULL;

    char *start = strstr(cvs_id, "Revision:");
    if (!start) return NULL;

    start += 9;
    while (*start > 0 && *start <= ' ') start++;

    if (start - cvs_id > strlen(cvs_id)) return NULL;

    char *end = start;
    while (*end > ' ') end++;

    if (*end == ' ') {
	int n = end-start+1;
	char *version = (char*)malloc(n*sizeof(char));
	strncpy(version, start, n-1);
	version[n-1] = 0;
	return version;
    } else {

	// Didn't find another space.  Assume no version was found.

	char *version = (char*)malloc(4*sizeof(char));
	strcpy(version, "0.0");
	return version;
    }
}

void parse_command_line(int argc, char** argv, plot_params* params)
{
    int i, j;
    double temp;

    /* Parse the argument list. */

    for (i = 1; i < argc; i++)

	if (argv[i][0] == '-') {

	    switch (argv[i][1]) {

	        case '-':	if (!strcmp(argv[i], "--help")) {
		    		    print_help(NULL);
				    exit(0);

				} else if (!strcmp(argv[i], "--version")) {
				    printf("plot_data (CVS ID %s)\n\n",
					   get_version("$Revision$"));

				    // GNU boilerplate:

				    printf(
"Copyright (C) 1994-2004, the Starlab development group.\n\
This is free software; see the source for copying conditions\n\
There is NO warranty; not even for MERCHANTABILITY or FITNESS\n\
A PARTICULAR PURPOSE.\n"
					);
					exit(0);
				}

		case 'c':	if (argv[i][2] == 'x')		       /*  cx */
				    params->crop_x = 1;
				else if (argv[i][2] == 'y')	       /*  cy */
		    		    params->crop_y = 1;
				else if (argv[i][2] == 'c') {
				    if (argv[i][3] == 'x')	       /* ccx */
					params->crop_x = 0;
				    else if (argv[i][3] == 'y')	       /* ccy */
					params->crop_y = 0;
				    else {
					params->crop_x = params->crop_y = 1;
				    }
				} else if (i < argc-1) {
		    		    params->xcol = atoi(argv[++i]);
				    if (params->xcol <= 0) {
					fprintf(stderr, USAGE);
					if (!in_line) exit(1);
				    }
				    params->ny = 0;
				    while (i < argc-1
					   && sscanf(argv[i+1], "%d", &j) > 0) {
					if (j == params->xcol) {
					    fprintf(stderr, USAGE);
					    if (!in_line) exit(1);
					} else if (j > 0) {
					    params->ycol[(params->ny)++] = j;
					    if (params->ny > 10) err_exit(
					   "Maximum of 10 y-columns allowed.");
					    if (params->ny > 1)
					    strcpy(params->color[params->ny-1],
						  params->color[params->ny-2]);
					    params->icolor[params->ny-1] = -1;
					}
					i++;
				    }
				} else
				    fprintf(stderr,
					    "-c: No parameters specified.\n");

				params->crop_data = (params->crop_x
						      || params->crop_y);
				if (params->crop_data) params->wrap_data = 0;

				break;

		case 'C':	if (i < argc-1) {
				    for (j = 0; j < params->ny; j++) {
					if (i < argc-1 && argv[i+1][0] != '-')
					    strcpy(params->color[j], argv[++i]);
					else if (j > 0)
					    strcpy(params->color[j],
						   params->color[j-1]);
					params->icolor[j] = -1;
				    }
				} else
				    fprintf(stderr,
					    "-C: No parameter specified.\n");

				break;

		case 'e':	if (argv[i][2] <= ' ') params->echo = 1;
				break;

		case 'h':       if (!strcmp(argv[i], "-help")) {
		    		    print_help(NULL);
				    exit(0);
				}
		                if (i < argc-1) {
				    strcpy(params->header, argv[++i]);
				} else
				    fprintf(stderr,
					    "-h: No parameter specified.\n");

				break;

		case 'i':	params->ignore_inline =
					1 - params->ignore_inline;
				break;

		case 'l':	if (argv[i][2] == 'l')
				    params->lines = 1;
				else {
				    if (i < argc-1) {

					set_or_flag(argv[++i],
						   &params->lims.xmin_from_data,
						   &params->lims.xmin);

					if (i < argc-1)
					    set_or_flag(argv[++i],
						   &params->lims.xmax_from_data,
						   &params->lims.xmax);

					if (i < argc-1)
					    set_or_flag(argv[++i],
						   &params->lims.ymin_from_data,
						   &params->lims.ymin);

					if (i < argc-1)
					    set_or_flag(argv[++i],
						   &params->lims.ymax_from_data,
						   &params->lims.ymax);

					if (!params->lims.xmin_from_data
					    && !params->lims.xmax_from_data
					    && params->lims.xmax
					        <= params->lims.xmin)
					    swap(&params->lims.xmax,
						 &params->lims.xmin);

					if (!params->lims.ymin_from_data
					    && !params->lims.ymax_from_data
					    && params->lims.ymax
					        <= params->lims.ymin)
					    swap(&params->lims.ymax,
						 &params->lims.ymin);

					params->limits_from_data = 
					    params->lims.xmin_from_data
					    || params->lims.xmax_from_data
					    || params->lims.ymin_from_data
					    || params->lims.ymax_from_data;

					if (!params->limits_from_data)
					    if (params->lims.xmax
					        <= params->lims.xmin
						|| params->lims.ymax
					        <= params->lims.ymin)
						err_exit("Illegal limits");

					params->yscale = -1.0;
				    }
				}

				break;

		case 'L':	if (i < argc-1) {
				    temp = atof(argv[++i]);
				    if (temp <= 0) temp = 1.0;	   /* "-L" OK */
				} else
				     temp = 1.0;

				params->lims.xmin = -temp;
				params->lims.xmax = temp;
				params->lims.ymin = -temp;
				params->lims.ymax = temp;

		    		params->limits_from_data = 0;
				params->lims.xmin_from_data = 0;
				params->lims.xmax_from_data = 0;
				params->lims.ymin_from_data = 0;
				params->lims.ymax_from_data = 0;

				params->yscale = -1.0;
				break;

		case 'N':	if (i < argc-1)
				    params->nmax = atoi(argv[++i]);
				else
				    fprintf(stderr,
					    "-N: No parameter specified.\n");

				break;

		case 'o':	params->output = 1 - params->output;
				break;

		case 'O':	if (i < argc-1)
				    params->xorigin = atoi(argv[++i]);
		                if (i < argc-1)
				    params->yorigin = atoi(argv[++i]);
				break;

		case 'p':	if (argv[i][2] == 'p')
		    		    params->lines = 0;
		    		else
				    params->lines = 1 - params->lines;

				break;

		case 'P':	if (i < argc-1)
				    params->point_size = atof(argv[++i]);
				if (params->point_size < 0.0)
				    params->point_size = 0.0;
				break;

		case 'q':	params->quiet = 1 - params->quiet;
				break;

		case 'Q':	exit(0);

		case 's':	if (i < argc-1)
				    params->xsize = atoi(argv[++i]);
				if (i+1 < argc && argv[i+1][0] != '-')
				    params->ysize = atoi(argv[++i]);
				else
				    params->ysize = params->xsize;
				break;

		case 'S':	if (i < argc-1)
				    params->skip = atoi(argv[++i]);
				break;

		case 't':	if (i < argc-1) {
				    j = atoi(argv[++i]);
				    if (j < 0) j = 0;		   /* "-t" OK */
				} else
				     j = 0;
				params->ntrail = j;
				break;

		case 'W':	params->pause = 1;
				break;

		case 'w':	if (argv[i][2] == 'x')		       /*  wx */
		    		    params->wrap_x = 1;
				else if (argv[i][2] == 'y')	       /*  wy */
		    		    params->wrap_y = 1;
				else if (argv[i][2] == 'w') {
				    if (argv[i][3] == 'x')	       /* wwx */
					params->wrap_x = 0;
				    else if (argv[i][3] == 'y')	       /* wwy */
					params->wrap_y = 0;
				}

				params->wrap_data = (params->wrap_x
						      || params->wrap_y);

				break;

		case 'x':	if (i < argc-1) {
				    strcpy(params->xlabel, argv[++i]);
				    params->xlabel_set = 1;
				    params->yscale = -1.0;
				} else
				    fprintf(stderr,
					    "-x: No parameter specified.\n");

				break;

		case 'X':	params->yscale = -1.0;
				break;

		case 'y':	if (i < argc-1) {
				    strcpy(params->ylabel, argv[++i]);
				    params->ylabel_set = 1;
				    params->yscale = -1.0;
				} else
				    fprintf(stderr,
					    "-y: No parameter specified.\n");

				break;

		case 'z':	if (i < argc-1)
				    params->zcol = atoi(argv[++i]);
				break;

		case 'H':
		default:	if (!params->quiet) print_help(argv[i]);
				if (strcmp(argv[0], "INLINE")) exit(1);

	    }

	}
}

void echo_parameters(plot_params* params)
{
    int i;

    if (!params->echo) return;

    fprintf(stderr, "\nCurrent parameter settings:\n");
    fprintf(stderr, "    xcol = %d, ycol =", params->xcol);
    if (params->ny > 0)
	for (i = 0; i < params->ny; i++)
	    fprintf(stderr, " %d (%s)", params->ycol[i], params->color[i]);
    else
	fprintf(stderr, " (not set)");
    fprintf(stderr, ", zcol = %d\n", params->zcol);

    fprintf(stderr, "    xmin = %f, xmax = %f, ymin = %f, ymax = %f\n",
	   params->lims.xmin, params->lims.xmax,
	   params->lims.ymin, params->lims.ymax);
    fprintf(stderr, "    crop_x = %d, crop_y = %d\n",
	    params->crop_x, params->crop_y);
    fprintf(stderr, "    wrap_x = %d, wrap_y = %d\n",
	    params->wrap_x, params->wrap_y);

    fprintf(stderr, "    x-label = \"%s\", y-label = \"%s\", header = \"%s\"\n",
	   params->xlabel, params->ylabel,
	   (params->header ? params->header : "(none)"));

    fprintf(stderr, "    lines = %d, ntrail = %d, point_size = %f\n",
	   params->lines, params->ntrail, params->point_size);

    fprintf(stderr, "    quiet = %d, ignore_inline = %d, output = %d\n",
	   params->quiet, params->ignore_inline, params->output);

    fprintf(stderr, "\n");
    params->echo = 0;
}

void create_labels(plot_params* params)
{
    int j;

    if (!params->xlabel_set) sprintf(params->xlabel, "column %d", params->xcol);
    if (!params->ylabel_set) {
	sprintf(params->ylabel, "column");
	if (params->ny > 1) strcat(params->ylabel, "s");
	for (j = 0; j < params->ny; j++) {
	    char temp[4];
	    sprintf(temp, " %d", params->ycol[j]);
	    strcat(params->ylabel, temp);
	}
    }
}

void allocate_save(plot_params* params)
{
    if ((params->xsave = (float*)malloc(params->ntrail*sizeof(float)))
	== NULL)
	err_exit("Can't allocate x-storage space.");
    if ((params->ysave = (float*)malloc(NYMAX*params->ntrail*sizeof(float)))
	== NULL)
	err_exit("Can't allocate y-storage space.");
}

void skip_lines(int skip, int output)
{
    /* Place pointer after the skip-th newline. */

    int i;

    while (skip-- > 0) {
	while ((i = getchar()) != '\n' && i != EOF)
	    if (output) putchar(i);
	if (output) putchar('\n');
    }
}

void split_inline_command(char* temp, int* argc, char** argv)
{
    *argc = 0;
    argv[(*argc)++] = "INLINE";

    while (*temp != '\0') {
	while (*temp == ' ') temp++;
	argv[(*argc)++] = temp;
	while (*temp != ' ' && *temp != '\0') temp++;
	if (*temp != '\0') {
	    *temp = '\0';
	    temp++;
	}
    }    
}

void parse_inline_command(char* temp, plot_params* params)
{
    int argc;
    char* argv[256];	/* Seems like a conservative limit! */

    /* Parse an inline command.  Assume syntax identical to the
       UNIX command-line switches (use parse_command_line to do
       the work).
     */

    split_inline_command(temp, &argc, argv);
    parse_command_line(argc, argv, params);
}

void set_window_parameters(plot_params* params)
{
    /* Draw a box, with tick marks. */

    lux_setup_region(params->xwin, 1.5, 1.5, 7.5, 7.5);
    lux_setup_axis(params->xwin,
		   params->lims.xmin, params->lims.xmax,
		   params->lims.ymin, params->lims.ymax);
    lux_draw_axis(params->xwin);

    /* Add some (optional) labels. */

    if (params->xlabel)
	lux_draw_string(params->xwin,
			0.5*(params->lims.xmin + params->lims.xmax),
			params->lims.ymin - 0.075*(params->lims.ymax
						   - params->lims.ymin),
			-1., params->xlabel, 0);
    if (params->ylabel)
	lux_draw_vstring(params->xwin,
			 params->lims.xmin - 0.077*(params->lims.xmax
						    - params->lims.xmin),
			 0.5*(params->lims.ymin + params->lims.ymax),
			 0., params->ylabel, 0);

    /* Add an optional header. */

    if (params->header)
	lux_draw_string(params->xwin,
			0.5*(params->lims.xmin + params->lims.xmax),
			params->lims.ymax + 0.05*(params->lims.ymax
						  - params->lims.ymin),
			0., params->header, 0);
}

void process_inline_command(char *input_line,
			    plot_params* params, int* reinit)
{
    /* Decode an embedded command line. */

    char xlabel[256], ylabel[256];

    if (!params->limits_from_data) {

	int xcol, ycol[NYMAX], zcol, iy, ny, ntrail;	/* Local copies */

	xcol = params->xcol;
	ny = params->ny;
	for (iy = 0; iy < params->ny; iy++) ycol[iy] = params->ycol[iy];
	zcol = params->zcol;
	ntrail = params->ntrail;

	strcpy(xlabel, params->xlabel);
	strcpy(ylabel, params->ylabel);

 	if (!params->quiet)
	    fprintf(stderr, "Read inline command \"%s\"\n", input_line);

	parse_inline_command(input_line, params);

	if (params->ntrail != ntrail) {
	    if (ntrail == 0) {
		arr_size = 1;
		allocate_save(params);
	    }
		
	    if (params->ntrail == 0) {
		arr_size = init_arr_size;
		free(params->xsave);
		free(params->ysave);
	    }
	}

	/* Force recomputation of read_list and index (saved as static
	   variables for efficiency) if xcol, ycol, or zcol are modified. */

	if (params->xcol != xcol || params->ny != ny || params->zcol != zcol)
	    *reinit = -1;
	else
	    for (iy = 0; iy < params->ny; iy++)
		if (params->ycol[iy] != ycol[iy]) *reinit = -1;

	if (*reinit == -1) create_labels(params);

	if (strcmp(xlabel, params->xlabel) || strcmp(xlabel, params->xlabel))
	    params->yscale = -1;

	/* (Setting yscale < 0 forces the box to be redrawn.) */

	/* Check to see if any color was changed. */

	for (iy = 0; iy < params->ny; iy++)
	    if (params->icolor[iy] < 0)
		params->icolor[iy] = lux_lookup_color(params->xwin,
						      params->color[iy]);

	/* Check to see if limits were changed. */

	if (params->yscale < 0) {

	    lux_reset_window(params->xwin);

	    params->curr_color = params->ibox;
	    lux_set_color(params->xwin, params->curr_color);
	    set_window_parameters(params);

	    params->yscale = (params->lims.ymax - params->lims.ymin)
		/ (params->lims.xmax - params->lims.xmin);
	}

	if (params->pause) {

	    fprintf(stderr,
		    "Press any key in display window to continue\n");
	    while(!win_getkey(params->xwin));

	    params->pause = 0;
	}

	if (params->skip) {
	    skip_lines(params->skip, params->output);
	    params->skip = 0;
	}

	if (params->echo) echo_parameters(params);

    } else

	if (!params->quiet)
	    fprintf(stderr, "Read and skipped inline command \"%s\"\n",
		    input_line);
}

void force_limits(float* xmin, int set_xmin, float* xmax, int set_xmax)
{
    double dx, dy, z;
    int idy, i;

    if (*xmax <= *xmin) return;

    /* Simple criterion: take the next "tick-mark" up or down, in units
       of the power of 10 closest to one-twentieth the given interval. */

    dx = 0.05 * (*xmax - *xmin);

    /* Find the power of 10 closest to dx. */

    dy = log10(dx);
    idy = dy + 0.5;
    if (dy < -0.5) idy = idy - 1;

    /* Don't trust floating-point work with floats...

    dy = pow(10., (float)idy);
    fprintf(stderr, "idy = %d, dy = %f\n", idy, dy);

    */

    dy = 1.0;
    if (idy > 0)
	for (i = 0; i < idy; i++) dy *= 10.0;
    else if (idy < 0)
	for (i = 0; i < -idy; i++) dy /= 10.0;

    /* fprintf(stderr, "idy = %d, dy = %f\n", idy, dy); */

    /* Start at zero, go up above xmax, then down below xmin. */

    /* fprintf(stderr, "*xmin, *xmax = %f %f\n", *xmin, *xmax); */

    for (z = 0.0; z < *xmax + 0.1*dx; z += dy);
    for (; z > *xmin - 0.1*dx; z -= dy);
    if (set_xmin) *xmin = z;

    for (; z < *xmax + 0.1*dx; z += dy);
    if (set_xmax) *xmax = z;

    /* fprintf(stderr, "*xmin, *xmax = %f %f\n", *xmin, *xmax); */
}

void get_limits(float* x, float* y, int n, plot_params* params)
{
    /* Note: y is actually a 2D array, with ny rows, each of length n,
     * stored in a rectangular grid of total length arr_size. */

    int i, j;

    /* Determine the lower and upper limits on the data. */

    if (params->lims.xmin_from_data) params->lims.xmin = INFINITY;
    if (params->lims.xmax_from_data) params->lims.xmax = -INFINITY;
    if (params->lims.ymin_from_data) params->lims.ymin = INFINITY;
    if (params->lims.ymax_from_data) params->lims.ymax = -INFINITY;

    for (i = 0; i < n; i++) {
	if (params->lims.xmin_from_data
	     && x[i] < params->lims.xmin) params->lims.xmin = x[i];
	if (params->lims.xmax_from_data
	     && x[i] > params->lims.xmax) params->lims.xmax = x[i];
	for (j = 0; j < params->ny; j++) {
	    float yy = *(y+j*arr_size+i);
	    if (params->lims.ymin_from_data
		 && yy < params->lims.ymin) params->lims.ymin = yy;
	    if (params->lims.ymax_from_data
		 && yy > params->lims.ymax) params->lims.ymax = yy;
	}
    }

    /* Force the limits to something sensible. */

    force_limits(&(params->lims.xmin), params->lims.xmin_from_data,
		 &(params->lims.xmax), params->lims.xmax_from_data);
    force_limits(&(params->lims.ymin), params->lims.ymin_from_data,
		 &(params->lims.ymax), params->lims.ymax_from_data);

    if (!params->quiet) {
	fprintf(stderr, "Selected x limits %f to %f,  ",
		params->lims.xmin, params->lims.xmax);
	fprintf(stderr, "y limits %f to %f\n",
		params->lims.ymin, params->lims.ymax);
    }
}

Window set_up_window(plot_params* params)
{
    /* Open an X-window. */
    
    params->xwin = lux_openwin(params->xorigin, params->yorigin,
			       params->xsize, params->ysize);

    if (params->xwin > 0) set_window_parameters(params);

    return params->xwin;
}

void swapi(int* i, int* j)
{
    int k;

    k = *i;
    *i = *j;
    *j = k;
}

void sort_columns(plot_params* params,
		  int* read_list, int* nread, int* index)
{
    int i, j, index1[NYMAX+2];

    /* List all columns and establish connections with actual data.
     *
     *		index1 = 0	 ==>  xcol
     *		index1 = 1 -- ny ==>  ycol
     *		index1 = ny+1	 ==>  zcol
     */

    *nread = 1;
    index1[0] = 0;
    read_list[0] = params->xcol;
    
    for (j = 0; j < params->ny; j++) {
	index1[(*nread)] = j+1;
	read_list[(*nread)++] = params->ycol[j];
    }

    if (params->zcol > 0) {
	index1[(*nread)] = params->ny + 1;
	read_list[(*nread)++] = params->zcol;
    }

    /* Sort the columns, retaining connections. */

    for (i = 0; i < *nread; i++)
	for (j = i + 1; j < *nread; j++)
	    if (read_list[j] < read_list[i]) {
		swapi(read_list+i, read_list+j);
		swapi(index1+i, index1+j);
	    }

    /* Finally, invert index1 to obtain index. */

    for (i = 0; i < *nread; i++) index[index1[i]] = i;
}

#define MAX_LINE 1024

int readxyz(plot_params* params, float* x, float* y, float* z,
	    int output)
{
    /* Read x, y (array), and possibly z, from stdin.  Return 1 if
       successful, 0 otherwise.  Also allow special processing of
       control lines embedded in the data. */

    static int read_list[NYMAX+2], index[NYMAX+2], nread, col_init = -1;
    int start[NYMAX+2];

    float temp[NYMAX+2];
    int i, j, last, p, len;

    /* If an EOF is encountered, return 0, which terminates the read.
     *
     * If a line cannot be interpreted as numeric data, attempt to read
     * an inline command from it, act on the command, if possible, then
     * return 0 to force all data up to this point to be plotted and all
     * buffers to be flushed.  Graphics changes will take effect after
     * the NEXT read.
     *
     * The "modify" variable controls when changes can be made.
     * It is 0 if the buffers contain current data, 1 otherwise.
     * Only allow changes if modify = 1.
     *
     * Otherwise, return 1 to continue the read on the next line.
     */

    /* Problems with reading and interpreting data from stdin...
     * Read in an entire line and handle it as a string.
     */

    char input_line[MAX_LINE], *c;

    /* Rules:	fgets > 0 is the first character read.
     *		fgets = 0 ==> error
     *		fgets < 0 ==> EOF (= -1, usually)
     *
     * Return for fgets <= 0.  Otherwise, input_line is a null-terminated
     * string read containing the input line (including the '\n'). 
     */

    if (fgets(input_line, MAX_LINE, stdin) <= 0) return 0;

    /* We have read a valid line from stdin.  Check for embedded
     * commands, then read the data.
     *
     * Embedded commands always start with "-", just like command-line
     * arguments.
     */

    while (input_line[0] == '-' && input_line[1] >= 'A') {

	/* Line appears to contain an embedded inline command.  Read and
	 * decode it if and only if ignore_inline = 0, limits_from_data = 0,
	 * and modify = 1.
	 */ 

	if (!params->ignore_inline &&
	    !params->limits_from_data &&
	    !params->modify) {

	    /* Force buffers to be flushed before rereading line */

	    params->modify = 1;
	    return 0;
	}

	/* fprintf(stderr, "line length = %d\n", strlen(input_line)); */

	/* Need to strip any trailing newline... */

	i = strlen(input_line);
	if (i > 0 && input_line[i-1] == '\n') input_line[i-1] = '\0';

	process_inline_command(input_line, params, &col_init);

	/* Check for expose/resize events and refresh the display. */

	if (params->xwin > 0) win_checkevent(params->xwin);

	/* Read the next line. */

	if (fgets(input_line, MAX_LINE, stdin) <= 0) return 0;
    }

    /* Deal with unwanted punctuation. */

    for (c = input_line; *c != '\0'; c++)
	if (*c == ',') *c = ' ';

    /* Line apparently contains data.  Read it. */

    /* Maintain two static arrays to help in reading data:
     * 
     *	(1) read_list, an ordered list of columns to read (length nread),
     *	(2) index, to assign each column to the correct variable.
     *
     * Most of the complication here is because we may want to read
     * several columns, and that xcol, ycol and zcol are not in any
     * particular order.
     *
     * Columns to be read are read_list[0] < read_list[1] < ...
     *
     * Temp column index[0]	   ==> xcol
     *                  [1 -- ny]  ==> ycol
     *                  [ny+1]     ==> zcol
     */

    /* (Re)sort xcol, ycol (and possibly zcol), if necessary. */

    if (col_init < 0) {
	sort_columns(params, read_list, &nread, index);
	col_init = 1;
    }

    /* Note that we use FORTRAN numbering in specifying columns (in xcol,
     * ycol, etc.): the first column is 1, NOT 0!  We will continue this
     * convention here as we read in the temp array.
     */

    /* Read the data.   Start by locating the starting points of the
     * nread columns we want (avoids extra scanfs, note).  Note that
     * we don't check for end-of-line after every character.
     */

    len = strlen(input_line);

    /* Move to the start of the first column. */

    p = 0;
    while (input_line[p] <= ' ') p++;

    last = 0;
    for (j = 0; j < nread; j++) {

	if (p > len) return 0;

	for (i = last+1; i < read_list[j]; i++) {    /* skip unwanted data */
	    while (input_line[p] > ' ') p++;
	    while (input_line[p] <= ' ') p++;
	}

	/* Store this p for future use. */

	start[j] = p;

	/* Move to the start of the next column. */

	while (input_line[p] > ' ') p++;
	while (input_line[p] <= ' ') p++;

	last = read_list[j];
    }

    /* Now actually read the data. */

    for (j = 0; j < nread; j++)
	if (sscanf(input_line+start[j], "%f", temp+j) <= 0) return 0;

    /* Optionally print out the line only if we were successful in
     * reading it.
     */

    if (output) {
	printf("%s", input_line);
	fflush(stdout);
    }

    /* Redistribute the data. */

    *x = temp[index[0]];

    for (j = 0; j < params->ny; j++) *(y+j*arr_size) = temp[index[j+1]];

    if (params->zcol > 0) *z = temp[index[params->ny+2]];

    params->modify = 0;

    return 1;
}

void get_data(plot_params* params,
	      float* x, float* y, float* z, int* n)
{
    /* Read the input data. Note that the remainder of each input
     * line is discarded once the desired columns are read in.
     * The current read is terminated on encountering EOF or an
     * unreadable line (readxyz returns 0).
     */

    /* To avoid wasting too much space, we probably should buffer the
     * input and use malloc to increase the array sizes in the case
     * arr_size = NMAX (limits_from_data = 1)...
     *
     * Note also: y is actually a 2D array, with ny rows, each of length,
     * n, stored in a rectangular grid of total declared length arr_size.
     */

    params->modify = 1;
    *n = 0;

    while (*n < arr_size && readxyz(params, x + *n, y + *n, z + *n,
				    params->output)) {

#if 0
	/* This version would just echo the columns being plotted.
	   Probably better to echo everything, do plot_data can be
	   used in a pipeline. */

	if (params->output) {
	    int j;
	    printf("%f", *(x+*n));
	    for (j = 0; j < params->ny; j++) printf(" %f", *(y+j*arr_size+*n));
	    if (params->zcol > 0) printf(" %f", *(z+*n));
	    printf("\n");
	}
#endif

	(*n)++;
    }
}

void crop(float*x, float xmin, float xmax)
{
    if (*x < xmin) *x = xmin;
    if (*x > xmax) *x = xmax;
}

void crop_data(plot_params* params, float* x, float* y, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	if (params->crop_x) crop(x+i, params->lims.xmin, params->lims.xmax);
	if (params->crop_y) crop(y+i, params->lims.ymin, params->lims.ymax);
    }
}

void wrap(float*x, float xmin, float xmax)
{
    float dx = xmax - xmin;
    if (dx <= 0) return;

    while (*x < xmin) *x += dx;
    while (*x > xmax) *x -= dx;
}

void wrap_data(plot_params* params, float* x, float* y, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	if (params->wrap_x) wrap(x+i, params->lims.xmin, params->lims.xmax);
	if (params->wrap_y) wrap(y+i, params->lims.ymin, params->lims.ymax);
    }
}

int decode_color(Window xwin, float z)
{
    /* Completely arbitrary color scheme:

       		0 = white
		1 = black
		2 = red
		3 = green
		4 = blue
		5 = yellow
		6 = orange
		7 = purple
		8 = cyan
		9 = grey
    */

    if (z < 0) z = -z;
    while (z > 9.5) z -= 10;

    if (z < 0.5)      return lux_lookup_color(xwin, "white");
    else if (z < 1.5) return lux_lookup_color(xwin, "black");
    else if (z < 2.5) return lux_lookup_color(xwin, "red");
    else if (z < 3.5) return lux_lookup_color(xwin, "green");
    else if (z < 4.5) return lux_lookup_color(xwin, "blue");
    else if (z < 5.5) return lux_lookup_color(xwin, "yellow");
    else if (z < 6.5) return lux_lookup_color(xwin, "orange");
    else if (z < 7.5) return lux_lookup_color(xwin, "purple");
    else if (z < 8.5) return lux_lookup_color(xwin, "cyan");
    else if (z < 9.5) return lux_lookup_color(xwin, "grey");
    else              return lux_lookup_color(xwin, "black");
}

void plot_points(plot_params* params, float* x, float* y, int n)
{
    /* Plot one or more points in the desired style.
     * Note: point_size here is point size in x units.
     */

    if (params->point_size <= 0.0)

	    lux_draw_pointsf(params->xwin, x, y, n, 0);

    else {

	int i;
	for (i = 0; i < n; i++)
		lux_draw_arcf(params->xwin,
			      x[i] - params->point_size/2,
			      y[i] - params->point_size/2, 
			      params->point_size,
			      params->point_size * params->yscale,
			      0.0, 360.0);
    }

/*  lux_flush(win); */		/* Appears unnecessary, so long as the
    				   program generating the data forces
				   output and flushes the data properly. */

}

#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>

/* In general, UNIX systems are very vague about what value is assigned
   to RAND_MAX.  For the old rand() function, it was 32767.  However, this
   been superseded by random(), which should have RAND_MAX = 2147483647.
   Unfortunately, even in systems with the new random(), the old value of
   RAND_MAX is often still defined (e.g. DEC UNIX, Solaris 2,...).  For
   our purposes here, it doesn't much matter, but it *is* important that
   a properly normalized final result be obtained. */

float myrandom()
{
    int i = random() % RAND_MAX;
    return ((float)i) / RAND_MAX;
}

#define XOR_MIN		75
#define XOR_RANGE	100
#define YOR_MIN		75
#define YOR_RANGE	100

void randomize_origin(plot_params* params)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    srandom((unsigned int)(tv.tv_usec * tv.tv_usec));

    if (params->xorigin == -1)
	params->xorigin = XOR_MIN + (int)(XOR_RANGE*myrandom());

    if (params->yorigin == -1)
	params->yorigin = YOR_MIN + (int)(YOR_RANGE*myrandom());
}

void initialize_params(int argc, char** argv, plot_params* params)
{
   /* Set default plotting parameters: */

    params->xwin = 0;

    params->xcol = 1;
    params->nmax = NMAX;
    params->ycol[0] = 2;
    params->ny = 1;			/* Only 1 y column	      */
    params->zcol = 0;			/* No z-coloring	      */

    params->lines = 1;			/* Plot lines	 	      */
    params->point_size = 0.0;		/* Point size (0 ==> 1 pixel) */
    params->ntrail = 0;			/* No trails		      */
    params->xsave = NULL;
    params->ysave = NULL;
    params->output = 0;			/* Don't pass data to stdout  */

    params->xlabel[0] = '\0';
    params->xlabel_set = 0;
    params->ylabel[0] = '\0';
    params->ylabel_set = 0;
    params->header[0] = '\0';		/* No overall label 	      */

    strcpy(params->color[0], "black");
    params->icolor[0] = 1;
    params->ibox = 1;
    params->iback = 0;			/* Always 0 since I can't change it! */
    params->curr_color = 1;

    params->limits_from_data = 1;	/* Get limits from input data */

    params->lims.xmin_from_data = 1;
    params->lims.xmax_from_data = 1;
    params->lims.ymin_from_data = 1;
    params->lims.ymax_from_data = 1;

    params->lims.xmin = 0;
    params->lims.xmax = 1;
    params->lims.ymin = 0;
    params->lims.ymax = 1;

    params->crop_data = 0;
    params->crop_x = 0;
    params->crop_y = 0;

    params->wrap_data = 0;
    params->wrap_x = 0;
    params->wrap_y = 0;

    params->yscale = 0.0;

    params->xorigin = -1;
    params->yorigin = -1;
    params->xsize = 400;
    params->ysize = 400;

    params->skip = 0;
    params->ignore_inline = 0;
    params->quiet = 0;
    params->modify = 1;
    params->pause = 0;

    params->echo = 0;

    parse_command_line(argc, argv, params);
    
    if (params->xorigin == -1 || params->yorigin == -1)
	randomize_origin(params);
}

main(int argc, char* argv[])
{
    int n, nt = 0;
    float *x, *y, *z;
    float tempx[2], tempy[NYMAX][2], tempz;   /* Used in connecting segments */
    int i, j;
    int erase = 0, save = 0;
    int nread = NBUFFER;

    /* Initialize all parameters. */

    plot_params params;
    initialize_params(argc, argv, &params);

    /* Echo some input data. */

    if (!params.quiet) {
	fprintf(stderr, "xcol = %d\n", params.xcol);
	fprintf(stderr, "ycol =");
	for (j = 0; j < params.ny; j++) fprintf(stderr, " %d (%s) ",
					       params.ycol[j], params.color[j]);
	fprintf(stderr, "\n");
	fprintf(stderr, "%s, ", (params.lines ? "line mode" : "point mode"));
	if (params.ntrail)
	    fprintf(stderr, "ntrail = %d\n", params.ntrail);
	else
	    fprintf(stderr, "no trails\n");
	if (params.limits_from_data)
	  fprintf(stderr, "Getting (some) plot limits from input data...\n");
    }

    if (params.ntrail > 0) {
	if (params.limits_from_data)
	    err_exit("Inconsistent parameters limits_from_data and ntrail.");

	allocate_save(&params);
    }

    if (params.xsave) nread = 1;	/* Plot data point by point */

    if (params.limits_from_data)
	arr_size = params.nmax;		/* Probably should reduce this    */
					/* and expand buffer as needed... */
    else
	arr_size = nread;

    init_arr_size = arr_size;

    /* Establish sufficient storage for all data arrays. */

    if ((x = (float*)malloc(arr_size*sizeof(float))) == NULL)
	err_exit("Can't allocate x-storage space.");
    if ((y = (float*)malloc(NYMAX*arr_size*sizeof(float))) == NULL)
	err_exit("Can't allocate y-storage space.");
    if ((z = (float*)malloc(arr_size*sizeof(float))) == NULL)
	err_exit("Can't allocate z-storage space.");

    /* NOTE that we will have to use ugly explicit pointer notation
       below because y and ysave are not declared as 2-D arrays... */

    /* Establish default labels if none specified. */

    create_labels(&params);

    if (params.skip) {
	skip_lines(params.skip, params.output);
	params.skip = 0;
    }

    if (params.echo) echo_parameters(&params);

    if (params.limits_from_data) {

	/* Obtain limits from data in advance of plotting anything.
	 * Note that, in this case, inline control statements are ignored.
	 */

	get_data(&params, x, y, z, &n);
	if (n <= 0) exit(1);

	if (!params.quiet)
	    fprintf(stderr, "Read %d points from standard input\n", n);

	get_limits(x, y, n, &params);

	if (params.lims.xmin >= params.lims.xmax
	    || params.lims.ymin >= params.lims.ymax) {
	    fprintf(stderr, "Inconsistent limits: %f %f %f %f\n",
		    params.lims.xmin,params.lims.xmax,
		    params.lims.ymin, params.lims.ymax);
	    exit(1);
	}
    }

    in_line = 1;

    params.yscale = (params.lims.ymax - params.lims.ymin)
			/ (params.lims.xmax - params.lims.xmin);

    /* Plot the data. */

    if (set_up_window(&params) <= 0)

	fprintf(stderr, "Error opening X-window!\n");

    else {

	/* Establish colors, now that X is properly initialized. */

	params.ibox = lux_lookup_color(params.xwin, "black");
	params.iback = lux_lookup_color(params.xwin, "white");
	for (j = 0; j < params.ny; j++)
	    params.icolor[j] = lux_lookup_color(params.xwin, params.color[j]);

	lux_set_color(params.xwin, params.curr_color);

	/* Two-in-one loop:  Read and plot, or just plot, the data. */

	if (!params.limits_from_data) n = 1;	/* Just to initialize... */

	while (n > 0) {

	    /* Read the next block of data, if we don't already have it.
	     * Note that normally the data are buffered for efficiency.
	     * In the case of "trailing" points, we read the data point
	     * by point.
	     */

	    if (!params.limits_from_data) {
		get_data(&params, x, y, z, &n);
		if (params.crop_data)
		    crop_data(&params, x, y, n);
		else if (params.wrap_data)
		    wrap_data(&params, x, y, n);
	    }

	    if (n > 0) {

		if (params.lines) {

		    /* ----------  Drawing lines  ---------- */

		    if (params.ntrail > 0 && erase) {

			/* Erase a previous line segment. */

			params.curr_color = params.iback;
			lux_set_color(params.xwin, params.curr_color);

			if (save == params.ntrail - 1) {

			    float erasex[2], erasey[2];

			    erasex[0] = *(params.xsave+save);
			    erasex[1] = *(params.xsave);
			    for (j = 0; j < params.ny; j++) {
				erasey[0] = *(params.ysave
					      + j*params.ntrail+save);
				erasey[1] = *(params.ysave
					      + j*params.ntrail);
				lux_draw_linesf(params.xwin,
						erasex, erasey, 2, 0);
			    }

			} else

			    for (j = 0; j < params.ny; j++)
				lux_draw_linesf(params.xwin, params.xsave+save,
						params.ysave
						  + j*params.ntrail+save,
						2, 0);
		    }

		    if (nt == 2) {

			/* Connect to the ends of the previous lines, setting
			   color from the z-array if necessary. */

			if (params.zcol > 0) {
			    params.curr_color = decode_color(params.xwin,
							     tempz);
			    lux_set_color(params.xwin, params.curr_color);
			}

			tempx[1] = x[0];
			for (j = 0; j < params.ny; j++) {

			    if (params.zcol <= 0
				&& params.icolor[j] != params.curr_color) {
				params.curr_color = params.icolor[j];
				lux_set_color(params.xwin, params.curr_color);
			    }

			    tempy[j][1] = *(y+j*arr_size);
			    lux_draw_linesf(params.xwin,
					    tempx, &tempy[j][0], nt, 0);
			}
		    }

		    /* Plot the new line segments. */

		    if (params.zcol <= 0)

			for (j = 0; j < params.ny; j++) {

			    if (params.icolor[j] != params.curr_color) {
				params.curr_color = params.icolor[j];
				lux_set_color(params.xwin, params.curr_color);
			    }

			    lux_draw_linesf(params.xwin, x, y+j*arr_size,
					    n, 0);
			}

		    else {

			/* Plot point by point, selecting color
			   from the z array. */

			for (i = 1; i < n; i++) {
			    int ic;

			    if ((ic  = decode_color(params.xwin, z[i-1]))
				  != params.curr_color) {
				params.curr_color = ic;
				lux_set_color(params.xwin, params.curr_color);
			    }

			    for (j = 0; j < params.ny; j++)
				lux_draw_linesf(params.xwin,
						x+i-1, y+j*arr_size+i-1, 2, 0);
			}
		    }

		} else {

		    /* ----------  Drawing points only  ---------- */

		    if (params.ntrail > 0 && erase) {

			/* Erase a previous point. */

			params.curr_color = params.iback;
			lux_set_color(params.xwin, params.curr_color);

			for (j = 0; j < params.ny; j++)
			    plot_points(&params, params.xsave+save,
					params.ysave+j*params.ntrail+save, 1);
		    }

		    if (params.zcol <= 0)

			for (j = 0; j < params.ny; j++) {

			    if (params.icolor[j] != params.curr_color) {
				params.curr_color = params.icolor[j];
				lux_set_color(params.xwin, params.curr_color);
			    }

			    plot_points(&params, x, y + j*arr_size, n);
			}

		    else {

			/* Plot point by point, selecting color
			   from the z array. */

			for (i = 0; i < n; i++) {
			    int ic;

			    if ((ic  = decode_color(params.xwin, z[i]))
				  != params.curr_color) {
				params.curr_color = ic;
				lux_set_color(params.xwin, params.curr_color);
			    }

			    for (j = 0; j < params.ny; j++)
				plot_points(&params, x + i,
					    y + j*arr_size + i, 1);
			}
		    }
		}

		/* Check for expose/resize events and refresh the display. */

		win_checkevent(params.xwin);

		if (params.ntrail > 0) {
			
		    /* Save past data in a circular buffer. */

		    *(params.xsave+save) = x[0];
		    for (j = 0; j < params.ny; j++)
			*(params.ysave+j*params.ntrail+save) = *(y+j*arr_size);

		    if (++save >= params.ntrail) {
			erase = 1;
			save = 0;
		    }
		}

	    }	/* End of "if (n > 0)..." */

	    if (params.limits_from_data) {

		/* Time to quit. */

		n = 0;		/* This will force exit from the while loop */

	    } else if (n > 0) {

		/* Save last point for connection to next segment. */

		nt = 2;
		tempx[0] = x[n-1];
		for (j = 0; j < params.ny; j++)
		    tempy[j][0] = *(y+j*arr_size+n-1);
		tempz = z[n-1];
	    }
	}

	/* Clean up only at end (quite expensive!).
	   Redraw the box in case it was partly erased. */

	if (params.curr_color != params.ibox)
	    lux_set_color(params.xwin, params.ibox);

	lux_draw_axis(params.xwin);

	if (params.curr_color != params.ibox)
	    lux_set_color(params.xwin, params.curr_color);

	if (params.output) fclose(stdout);

	/* Enter idle mode before quitting. */

	fprintf(stderr,
		"Press any key in display window to quit current plot\n");
	while(!win_getkey(params.xwin));

    }
}
