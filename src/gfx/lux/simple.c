
/* simple.c:	Simple C interface onto the lux libraries.
 *		 -- reinventing the wheel yet again...
 *
 *		Accessible functions:
 *
 *			draw_box
 *			add_label
 *			set_color
 *			plot
 *			move_to
 *			draw_to
 *			point
 *			get_mouse
 *			pause
 *                      clear_graphics
 *			exit_graphics
 *
 *		Less commonly used:
 *
 *			crop_gfx
 *			nocrop_gfx
 * 			wait_for_mouse
 * 			check_mouse
 *
 *		All positions are specified in "data" units, so the
 *		"draw" functions can only be used after the box has
 *		been constructed.
 */

#include <stdio.h>

typedef unsigned long Window;

static int	print	= 1;

static Window	xwin	= 0;
static int	crop	= 1;

void crop_gfx()   {crop = 1;}
void nocrop_gfx() {crop = 0;}

static int	xorigin	= 100;
static int	yorigin	= 100;
static int	xsize	= 400;
static int	ysize	= 400;

static float	xmin	= 0;
static float	xmax	= 0;
static float	ymin	= 0;
static float	ymax	= 0;

static float	xlast[2] = {0, 0};
static float	ylast[2] = {0, 0};

static void initialize_graphics()
{
    if ((xwin = lux_openwin(xorigin, yorigin, xsize, ysize)) <= 0)
	fprintf(stderr, "Can't initialize graphics\n");
}

void set_box(int xo, int yo, int xs, int ys)
{
    if (xo >= 0) xorigin = xo;
    if (yo >= 0) yorigin = yo;
    if (xs > 0) xsize = xs;
    if (ys > 0) ysize = ys;
}

void Draw_box(float xmin_in, float xmax_in, char *xlabel,
	      float ymin_in, float ymax_in, char *ylabel,
	      int draw)
{
    if (!xwin) initialize_graphics();

    /* Draw a box, with tick marks. */

    xmin = xmin_in;
    xmax = xmax_in;
    ymin = ymin_in;
    ymax = ymax_in;

    lux_setup_region(xwin, 1.5, 1.5, 7.5, 7.5);
    lux_setup_axis(xwin, xmin, xmax, ymin, ymax);

    if (draw) {

	lux_draw_axis(xwin);

	/* Add optional labels. */

	if (xlabel && *xlabel != '\0')
	    lux_draw_string(xwin,
			    0.5*(xmin + xmax), ymin - 0.075*(ymax - ymin),
			    -1., xlabel, 0);

	if (ylabel && *ylabel != '\0')
	    lux_draw_vstring(xwin,
			     xmin - 0.077*(xmax - xmin), 0.5*(ymin + ymax),
			     0., ylabel, 0);
    }
}

void draw_box(float xmin_in, float xmax_in, char *xlabel,
	      float ymin_in, float ymax_in, char *ylabel)
{
    Draw_box(xmin_in, xmax_in, xlabel, ymin_in, ymax_in, ylabel, 1);
}

void nodraw_box()
{
    Draw_box(0, 1, "", 0, 1, "", 0);
}

void add_label(char *label)
{
    if (xmax <= xmin) Draw_box(0, 1, "", 0, 1, "", 0);

    lux_draw_string(xwin,
		    0.5*(xmin + xmax), ymax + 0.05*(ymax - ymin),
		    0., label, 0);
}

void set_color(char *color)
{
    if (!xwin) initialize_graphics();

    if (color && *color != '\0')
	lux_set_color(xwin, lux_lookup_color(xwin, color));
}

static float xinterp(float x1, float y1, float x2, float y2, float y)
{
    if (y2 == y1)
	return x1;
    else
	return x1 + (y - y1) * (x2 - x1) / (y2 - y1);}

static float yinterp(float x1, float y1, float x2, float y2, float x)
{
    if (x2 == x1)
	return y1;
    else
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

static void plot_inside(float x1, float y1, float x2, float y2)
{
    /* Draw the portion of the line from (x1, y1)to (x2, y2) that
     * lies within the box. */

    float xpl[2], ypl[2];

    if (x1 >= xmin && x1 <= xmax && y1 >= ymin && y1 <= ymax) {
	xpl[0] = x1;
	ypl[0] = y1;
    } else if (x2 >= xmin && x2 <= xmax && y2 >= ymin && y2 <= ymax) {
	xpl[0] = x2;
	ypl[0] = y2;
	x2 = x1;
	y2 = y1;
	x1 = xpl[0];
	y1 = ypl[0];
    } else
	return;

    /* Plot from (x1, y1) = (xpl[0], ypl[0]) inside the box to (x2, y2). */

    if (x2 < xmin) {
	float yb = yinterp(x1, y1, x2, y2, xmin);
	if (yb < ymin) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymin);
	    ypl[1] = ymin;
	} else if (yb > ymax) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymax);
	    ypl[1] = ymax;
	} else {
	    xpl[1] = xmin;
	    ypl[1] = yb;
	}
    } else if (x2 > xmax) {
	float yb = yinterp(x1, y1, x2, y2, xmax);
	if (yb < ymin) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymin);
	    ypl[1] = ymin;
	} else if (yb > ymax) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymax);
	    ypl[1] = ymax;
	} else {
	    xpl[1] = xmax;
	    ypl[1] = yb;
	}
    } else {
	if (y2 < ymin) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymin);
	    ypl[1] = ymin;
	} else if (y2 > ymax) {
	    xpl[1] = xinterp(x1, y1, x2, y2, ymax);
	    ypl[1] = ymax;
	} else {
	    xpl[1] = x2;
	    ypl[1] = y2;
	}
    }

    lux_draw_linesf(xwin, xpl, ypl, 2, 0);

}

void plot(float *x, float *y, int n)
{
    if (!xwin) {

	if (print)
	    fprintf(stderr, "plot: must draw box before plotting data\n");
	print = 0;

    } else {

	if (crop) {

	    /* Apply cropping as we go. */

	    int i1 = 0, i2;

	    while (i1 < n) {

		while (i1 < n && (x[i1] < xmin || x[i1] > xmax
				  || y[i1] < ymin || y[i1] > ymax)) i1++;

		if (i1 >= n) return;

		/* Point #i1 is the first to lie inside the box. */

		i2 = i1;
		while (i2 < n && x[i2] >= xmin && x[i2] <= xmax
		       && y[i2] >= ymin && y[i2] <= ymax) i2++;

		/* Point #i2 is the first to lie outside the box. */

		/* Special treatment of line from i1-1 to i1. */

		if (i1 > 0)
		    plot_inside(x[i1-1], y[i1-1], x[i1], y[i1]);

		/* Draw from i1 to i2-1. */

		if (i2 > i1+1) lux_draw_linesf(xwin, x+i1, y+i1, i2-i1, 0);

		/* Special treatment of line from i2-1 to i2. */

		if (i2 < n)
		    plot_inside(x[i2-1], y[i2-1], x[i2], y[i2]);

		i1 = i2;
	    }

	} else

	    lux_draw_linesf(xwin, x, y, n, 0);
    }
}

void move_to(float x, float y)
{
    if (!xwin) {
	if (print)
	    fprintf(stderr, "plot: must draw box before plotting data\n");
	print = 0;
    } else {
	xlast[0] = x;
	ylast[0] = y;
    }
}

void draw_to(float x, float y)
{
    if (!xwin) {
	if (print)
	    fprintf(stderr, "plot: must draw box before plotting data\n");
	print = 0;
    } else {
	xlast[1] = x;
	ylast[1] = y;

	if (crop)
	    plot_inside(xlast[0], ylast[0], xlast[1], ylast[1]);
	else
	    lux_draw_linesf(xwin, xlast, ylast, 2, 0);

	xlast[0] = x;
	ylast[0] = y;
    }
}

void point(float x, float y, float point_size)
{
    if (!xwin) {
	if (print)
	    fprintf(stderr, "plot: must draw box before plotting data\n");
	print = 0;
    } else {
	lux_fill_arcf(xwin,
		      x - point_size/2,
		      y - point_size/2, 
		      point_size,
		      point_size * (ymax-ymin) / (xmax-xmin),
		      0.0, 360.0);
    }
}

int check_mouse()
{
    /* Return 0 if no mouse button has been pressed; return the number
       of the button (1, 2, 3) otherwise. */

    int i = lux_check_buttonpress(xwin);

    if (i < 0 || i > 3)
	return 0;
    else
	return i;
}

static int mouse_print = 1;

void get_mouse(float *x, float *y)
{

    /* Return cursor position on right mouse click. */

    if (mouse_print) {
	fprintf(stderr,
		"use left button to get cursor position information,\n");
	fprintf(stderr,
		"use right button to capture cursor coordinates\n");
	mouse_print = 0;
    }

    get_mouse_position(xwin, x, y);
    /* fprintf(stderr, "x = %f, y = %f\n", *x, *y); */
}

void wait_for_mouse()
{
    /* Wait for mouse press before continuing. */

    float xdum, ydum;

    if (!xwin) return;

    fprintf(stderr,
	    "\a\npress any mouse button in display window to continue\n");
    get_mouse_position(xwin, &xdum, &ydum);
}

void pause(int time)
{
    lux_pause(time);
}

void clear_graphics()
{
    /* Clear the graphics window. */

    lux_reset_window(xwin);
}

void exit_graphics()
{
    /* Enter idle mode before quitting.  Exit on any keypress. */

    if (!xwin) return;

    fprintf(stderr, "\a\npress any key in display window to exit\n");
    while(!win_getkey(xwin));
}


#ifdef TEST

/* Run a standard test program... */

#define N 6

main()
{
    float x[N] = {0.1, 0.9, 0.9, 0.5, 0.1, 0.1};
    float y[N] = {0.1, 0.1, 0.9, 2.2, 0.9, 0.1};

    set_color("blue");
    draw_box(0, 1, "x-axis", 0, 2, "y-axis");

    nocrop_gfx();
    set_color("black");
    plot(x, y, N);

    crop_gfx();
    set_color("yellow");
    plot(x, y, N);

    set_color("green");
    move_to(0.2, 0.5);
    draw_to(0.2, 0.6);
    draw_to(0.3, 0.6);
    draw_to(0.3, 0.5);
    draw_to(0.2, 0.5);

    set_color("red");
    add_label("hello, this is a test");

    set_color("purple");
    point(0.9, 1.0, .1);

    get_mouse(x, y);

    clear_graphics();
    set_color("pink");
    draw_box(0, (*x > 0 ? *x : 1), "x-axis",
	     0, (*y > 0 ? *y : 1), "y-axis");

    get_mouse(x, y);
    get_mouse(x, y);

    while(1) {
	int i = check_mouse();
	fprintf(stderr, "mouse = %d\n", i);
	if (i > 0) break;
	pause(100000);
    }

    exit_graphics();
}

#endif
