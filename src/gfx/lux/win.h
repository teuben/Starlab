/* win.h
 * Biao Lu                     biao@eagle.drexel.edu
 */

#ifndef LUX_WIN_H
#define LUX_WIN_H

#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <math.h>




#define register 




#define BITMAPDEPTH 1

#define SMALL 1
#define OK 0

#define GRAPH_WINDOW    1
#define POPUP_WINDOW    2
#define DIALOG_WINDOW   3
#define BUTTON_WINDOW   4
#define INPUT_WINDOW    5
#define TEXT_WINDOW     6    /* Actualy not */

#define NO_TYPE         0
#define OK_BUTTON       1
#define CANCEL_BUTTON   2
#define CHECK_BUTTON    3
#define OK_KEEP_BUTTON  4

/* For axis style */
#define BOXED           1
#define FRAMED          2
#define NORMAL          3

#define NOOP          	      0
#define DRAW_LINE     	      1
#define DRAW_LINES    	      2
#define DRAW_PONIT    	      3
#define DRAW_POINTS   	      4
#define SET_LINE_WIDTH	      5
#define SETUP_REGION  	      6
#define SETUP_AXIS    	      7
#define DRAW_AXIS     	      8
#define DRAW_LINE_F   	      9
#define DRAW_STRING           10
#define DRAW_VSTRING          11
#define DRAW_SEGMENTS_F       12
#define SET_COLOR             13
#define DRAW_POINT_F          14
#define CLEAR_CURRENT_REGION  15
#define CLEAR_WINDOW          16
#define DRAW_IMAGE_STRING     17
#define SET_BG_COLOR          18
#define UPDATE_FOREGROUND     19
#define DRAW_RECTANGLES_F     20
#define DRAW_LINES_F          21
#define DRAW_POINTS_F         22
#define DRAW_RECTANGLE_F      23
#define DRAW_ARCS_F           24
#define DRAW_ARC_F            25
#define FILL_ARC_F            26
#define FILL_ARCS_F           27
#define FILL_RECTANGLE_F      28
#define FILL_RECTANGLES_F     29
#define FILL_POLYGON_F        30
#define SET_LINE_STYLE        31
#define SET_WINDOW_BG_COLOR   32
#define SET_UPDATE            33
#define SET_NO_UPDATE         34
#define SETUP_AXIS_STYLE      35


typedef struct _lux_data {
  unsigned int type;
  union {
    char          *b;
    short         *s;
    int           *i;
    unsigned int  *u;
    long          *l;
    unsigned long *ul;
    float         *f;
    double        *d;
  }  data;
  struct  _lux_data  *next;
} lux_data;

typedef struct {

  int            type;   /* This is for checking the type of window */
  int            subtype;

  Display       *display;
  unsigned int   screen;

  unsigned int   display_width;
  unsigned int   display_height;

  unsigned int   window_depth;

  int            update_flag;    /* For animation  1 for update*/
  int            discard_flag;   /* For save data or not  0 for save*/

  int            x;
  int            y;
  unsigned int   width;
  unsigned int   height;
  unsigned int   old_width;      /* To determine if window size is changed*/  
  unsigned int   old_height;     /* or not */
  unsigned int   user_width;     /* This is to honor the size user specified */
  unsigned int   user_height;    /* The purpose is for resize */
  float          xresizefactor;
  float          yresizefactor;

  char          *window_name;
  char          *icon_name;
  int            window_size;    /* This is for checking the window size */

  Window         parent;
  Window         serial;         /* For item number, dialog attached window */
  Window         window;
  Pixmap         pixmap;
  Pixmap         icon_pixmap;

  XSizeHints     size_hints;
  GC             gc;

  Cursor         cursor;
  Colormap       colormap;
  char          *colormapfile;
  int            lux_colormap;

  lux_data      *data;
  lux_data      *currentdata;
  char           msg[255];          /* msg got after return from dialog */

/* The following are used for drawing axis */
  float          fxsize, fysize;   /* This for the size of axis 
				    * 10.0 same as window */
  float          fxorg;   /* This specify (x,y) origin in [0.0,10.0] rect */
  float          fyorg;   /* for axis box */
  float          xmin, xmax;
  float          ymin, ymax;
  float          xfactor, yfactor;  /* related to real data */
  int            xorg, yorg;        /* Real orgin for axis box */
  int            xsize, ysize;      /* size of axis in window */

  int            lnx, lny;          /* Log axis */

  int            buttonflag;  /* for after get button posion to leave or not */
  float          xbutton, ybutton;   /* for button position*/

} lux_win;

typedef struct _lux_wins {
  lux_win    win;
  struct _lux_wins  *next;
} lux_wins;

#endif  /*LUX_WIN_H*/
