/* win.c
 * Biao Lu                     biao@eagle.drexel.edu
 */

#include "win.h"
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#include <string.h>

/* Global variables: only 1 of each allowed in lux! */

Display      *display;
Visual       *visual;
char         *display_name = NULL;
unsigned int  display_width, display_height, default_depth;
Pixmap        default_icon_pixmap;
int           pmpallocflag = 1;
int           screen;
lux_wins     *windows = (lux_wins *)NULL;
XStandardColormap  map_info;
int                colormap_size;
char         *colorfile = (char *)NULL;
Colormap      colormap;
Window        popup;
Atom          protocols;   /* for delete */

XFontStruct  *font_info;
GC            defaultgc, redrawgc;
Cursor        cursor;
int           lux_colors = 0;
int           lux_colormap;


extern        redraw(), lux_setup_region(), lux_setup_axis(), 
              lux_reconvert_coord();

extern  Colormap   lux_setup_colormap();
extern  unsigned   long lux_rgb_pixel();
extern  unsigned   long lux_lookup_color();


#define REPORT_EVENT 0
#define DOUBLECLICKLENGTH 300L


/* Used for McDraw */

static  int  mcdrawflag = 0;

/* Used for McDraw */

static char* default_font = "fixed";

set_default_font(char* font)
{
    default_font = malloc(strlen(font)+1);
    strcpy(default_font, font);
}

/* The following are for input focus */
Window        focus;
int           revert;
/* The above     are for input focus */

long      XMAXREQUESTSIZE;

#define X_CreatePixmap 53   /* Got this from XErrorDB */

int lux_handler(display, lux_err)
Display      *display;
XErrorEvent  *lux_err;
{
    char msg[80];
    XGetErrorText(display, lux_err->error_code, msg, 80);
    fprintf(stderr, "luxwins: Error %s\n Request code %d\n", 
	    msg, lux_err->request_code);
    if (lux_err->request_code == X_CreatePixmap) pmpallocflag = 0;
}


load_font(font_info)
XFontStruct **font_info;
{
    char *fontname = default_font;

    /* Access font */

    if ((*font_info = XLoadQueryFont(display, fontname)) == NULL) {

	fprintf(stderr, "LUX: Cannot open font \"%s\"\n", fontname);
	fprintf(stderr, "     Trying font \"fixed\"\n");

	if ((*font_info = XLoadQueryFont(display, "fixed")) == NULL) {
	    if (mcdrawflag)
		return 0;
	    else
		exit(-1);
	}

    }

    /* fprintf(stderr, "font loaded\n");*/

    return 1;
}


lux_init_gc()
{
  static int  status = 0;
  XGCValues            values;
  unsigned long        valuemask   = 0; /* ignore XGCvalues and use defaults */

  values.graphics_exposures = False;

  defaultgc = XCreateGC(display, RootWindow(display, screen), 
			valuemask, &values);
  redrawgc  = XCreateGC(display, RootWindow(display, screen), 
			valuemask, &values);

  if (load_font(&font_info))
    XSetFont(display, defaultgc, font_info->fid);

  XSetForeground(display, defaultgc, BlackPixel(display, screen));
  XSetBackground(display, defaultgc, WhitePixel(display, screen));
    
  XChangeGC(display, defaultgc, GCGraphicsExposures, &values);
  return 0;
}

lux_get_wingc(win)
lux_win *win;
{
    XGCValues            values;
    unsigned long        valuemask = 0; /* ignore XGCvalues and use defaults */

    win->gc = XCreateGC(win->display, RootWindow(win->display, win->screen), 
			valuemask, &values);

    XCopyGC(win->display, defaultgc, (long)0x003FFFFF, win->gc);

}

lux_create_popup()
{
  
    XSetWindowAttributes   attrib;
    unsigned     long      valuemask;

    if (colormap != DefaultColormap(display, screen)) {
      attrib.colormap = colormap;
      attrib.background_pixel = WhitePixel(display, screen);
      attrib.border_pixel = BlackPixel(display, screen);
      attrib.save_under = True;
      valuemask = CWColormap | CWBackPixel | CWBorderPixel | CWSaveUnder ;
      popup = XCreateWindow(display, 
			  RootWindow(display, screen), 
			  0, 0, 40*font_info->max_bounds.width, 
			  font_info->ascent+font_info->descent, 2, 
			  default_depth, InputOutput, visual, valuemask, 
		          &attrib);
    }
    else {	     
      popup = XCreateSimpleWindow(display, 
				  RootWindow(display, screen), 
				  0, 0, 40*font_info->max_bounds.width, 
				  font_info->ascent+font_info->descent, 2, 
				  BlackPixel(display, screen), 
				  WhitePixel(display, screen));
      attrib.save_under = True;
      XChangeWindowAttributes(display, popup, CWSaveUnder, 
			      &attrib);

    }

    if (REPORT_EVENT) fprintf(stderr, "Popup window id = %ld \n", popup);
}

lux_xinit()
{
    static int status = 0;

    if( !status ) {

      XSetErrorHandler(lux_handler);
      if ( (display=XOpenDisplay(display_name)) == NULL )
	{
	  (void) fprintf( stderr, 
			 "Cannot connect to X server %s\n", 
			 XDisplayName(display_name));
	  if (mcdrawflag) return 1;
	  else exit( -1 );
	}
 
      screen = DefaultScreen(display);
      display_width  = DisplayWidth(display, screen);
      display_height = DisplayHeight(display, screen);
      default_depth  = DefaultDepth(display, screen);

      cursor = XCreateFontCursor(display, XC_crosshair);
      {
	/* To get the window ID where the program start from, not for sure */
	XGetInputFocus(display, &focus, &revert);
/*	fprintf(stderr, "Window focus in %d \n", focus);*/
      }
      XMAXREQUESTSIZE = XMaxRequestSize(display);
/*      fprintf(stderr, "XMaxRequestSize = %d \n", XMaxRequestSize(display);*/

      colormap = lux_setup_colormap(display, screen, &visual);
      if (colormap == DefaultColormap(display, screen)) {
	visual = DefaultVisual(display, screen);
/*	fprintf(stderr, "Use default visual\n");*/
      }
      protocols = XInternAtom(display, "WM_DELETE_WINDOW", True);
      lux_init_gc();
      lux_create_popup();
      status = 1;
     }
/*    fprintf(stderr, "Initialized\n");*/
    return 0;
}

/* Convenient accessors: */

int lux_get_display_width()
{
    return display_width;
}

int lux_get_display_height()
{
    return display_height;
}

Bool predproc(display, event, arg)
Display  *display;
XEvent   *event;
char     *arg;
{
    register Window win;

    win = *(Window *)arg;

    switch (event->type) {
        case Expose:
          if (event->xexpose.window != win) return False;
	  return True;
	  break;
	default:
	  return False;
	  break;
    }
}

lux_wait_windowexpose(win)
Window win;
{
    XEvent report;
    Window window;

    window = win;   /* Very tricky */
    XPeekIfEvent(display, &report, predproc, (char *)(&window));
}

Bool keypredproc(display, event, arg)
Display     *display;
XEvent   *event;
char        *arg;
{
    Window              win;
    char                buffer[255];
    KeySym              keysym;
    int                 buffersize;
    int                 charcount;
    XComposeStatus      compose;
    char                key;

    buffersize = 20;

    win = *(Window *)arg;
    key = arg[sizeof(Window)];

    switch (event->type) {
        case KeyPress:
          if (event->xkey.window != win) break;
	  charcount = XLookupString(&(event->xkey), buffer, buffersize, 
				    &keysym, &compose);
	  if (charcount == 0) break;
	  if (((keysym >= XK_KP_Space) && (keysym <= XK_KP_9))
	      || ((keysym >= XK_space) && (keysym <= XK_asciitilde))) {
	    buffer[charcount] = 0;
/*	    fprintf(stderr, "get %d  %s keypress in win %u\n", charcount, 
		   buffer, win);*/
	    if (buffer[0] == key) return True;
	  }
	  break;
	default:
	  return False;
	  break;
    }
    return False;
}

Bool lux_check_keypress(win, key)
Window win;
char key;
{

    XEvent report;
    char   s[20];

    ((Window *)s)[0] = win;
    s[sizeof(Window)] = key;
    return XCheckIfEvent(display, &report, keypredproc, s);
}


Bool lux_next_keypress(win, key, function, shift, control)
Window win;
char *key, *shift, *control, function[20];
{
    char                buffer[255];
    KeySym              keysym;
    int                 buffersize = 20;
    int                 charcount;
    XComposeStatus      compose;
    XEvent              report;

    *key = 0;
    *shift = 0;
    *control = 0;
    function[0] = 0;

    /* Really should take care of expose/resize events while waiting... */

    if (XCheckTypedWindowEvent(display, win, KeyPress, &report)) {

	charcount = XLookupString(&(report.xkey), buffer, buffersize, 
				  &keysym, &compose);
/*
	fprintf(stderr, "Charcount=%d, Keystate=%d, keysym=%d, key=%c\n", 
		charcount, report.xkey.state, keysym, buffer[0]);
*/
	if (charcount != 0 ) {
	    *key = buffer[0];
	    return 1;
	} else {
	    if      (keysym == XK_Up)    strcat(function, "Up");
	    else if (keysym == XK_Down)  strcat(function, "Down");
	    else if (keysym == XK_Right) strcat(function, "Right");
	    else if (keysym == XK_Left)  strcat(function, "Left");
	    else if (keysym == XK_R7)    strcat(function, "Home");
	    else if (keysym == XK_R9)    strcat(function, "PgUp");
	    else if (keysym == XK_R11)   strcat(function, "R11");
	    else if (keysym == XK_R15)   strcat(function, "PgDn");
	    else if (keysym == XK_End)   strcat(function, "End");
	    else return 0; /* Other keys are not yet handled */

	    if (report.xkey.state & ControlMask) *control = 1;
	    if (report.xkey.state & ShiftMask) *shift = 1;
	}
	return 1;
    }
    return 0;
}

Bool buttonpredproc(display, event, arg)
Display  *display;
XEvent   *event;
char     *arg;
{
    Window  win;

    win = *((Window *)arg);

    switch (event->type) {
        case ButtonPress:
          if (event->xbutton.window != win) break;
	  ((unsigned int *)arg)[0] = event->xbutton.button;
	  return True;
	  break;
	default:
	  return False;
	  break;
    }
    return False;
}

lux_check_buttonpress(win)
Window win;
{
    XEvent report;
    char   s[32];

    ((Window *)s)[0] = win;
    if (XCheckIfEvent(display, &report, buttonpredproc, s)) {
/*      fprintf(stderr, "Button pressed %d \n", 
	      ((unsigned int *)s)[0]);*/
      return (((unsigned int *)s)[0]);
    }
    else return -1;
}


Bool buttonpositionpred(display, event, arg)
Display  *display;
XEvent   *event;
char     *arg;
{
    Window  win;

    win = *((Window *)arg);

    switch (event->type) {
        case ButtonPress:
          if (event->xbutton.window != win) break;
	  ((unsigned int *)arg)[0] = event->xbutton.button;
	  ((int *)(&(arg[sizeof(Window)])))[0] = event->xbutton.x;
	  ((int *)(&(arg[sizeof(Window)])))[1] = event->xbutton.y;
	  return True;
	  break;
	default:
	  return False;
	  break;
    }
    return False;
}

lux_check_buttonposition(win, x, y)
Window win;
int *x, *y;
{

    XEvent report;
    char   s[sizeof(Window)+sizeof(int)*2];

    ((Window *)s)[0] = win;
    if (XCheckIfEvent(display, &report, buttonpositionpred, s)) {
/*      fprintf(stderr, "Button pressed %d \n", 
	      ((unsigned int *)s)[0]);*/
	*x = ((int *)(&(s[sizeof(Window)])))[0];
	*y = ((int *)(&(s[sizeof(Window)])))[1];
      return (((unsigned int *)s)[0]);
    }
    else return 0;
}

lux_clear_pixmap(display, pixmap, gc, x, y, width, height)
Display *display;
Drawable pixmap;
GC       gc;
unsigned int x, y, width, height;
{
    XGCValues values;

    XGetGCValues(display, gc, GCForeground | GCBackground, &values);

    XSetForeground(display, gc, values.background);
    XFillRectangle(display, pixmap, gc, x, y, width, height);
    XSetForeground(display, gc, values.foreground);
}


#define icon_bitmap_width 40
#define icon_bitmap_height 40

lux_createwin(win)
lux_win *win;
{
    
    static char default_win_name[255];
    static char default_icon_name[255];

    static unsigned char default_icon_bitmap_bits[] = {
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
	0x00, 0x00, 0x60, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x08, 0x00, 0x60, 0x00, 
	0x00, 0x30, 0x00, 0x60, 0x00, 0x00, 0xe0, 0x00, 0x40, 0x00, 0x00, 0xc0, 
	0x00, 0x00, 0x00, 0x00, 0xc0, 0x05, 0x00, 0xc0, 0x00, 0x80, 0x03, 0x00, 
	0x00, 0x00, 0x80, 0x13, 0x00, 0x00, 0x00, 0x00, 0x8f, 0x00, 0x00, 0x00, 
	0x00, 0x9e, 0x00, 0x00, 0x02, 0x00, 0x7c, 0x48, 0x00, 0x00, 0x00, 0xf8, 
	0xfe, 0x00, 0x00, 0x00, 0xf0, 0xf8, 0x01, 0x00, 0x00, 0xf0, 0xf3, 0x03, 
	0x00, 0x08, 0xe0, 0xe3, 0x03, 0x00, 0x00, 0xc0, 0xd7, 0x03, 0x00, 0x10, 
	0xe0, 0x8f, 0x07, 0x00, 0x10, 0xa0, 0x5f, 0x06, 0x00, 0x00, 0xe0, 0x3f, 
	0x04, 0x10, 0x00, 0xc0, 0xff, 0x01, 0x00, 0x00, 0x40, 0xff, 0x05, 0x00, 
	0x00, 0x80, 0xff, 0x03, 0x40, 0x00, 0x80, 0xfe, 0x0b, 0x00, 0x00, 0x00, 
	0xfa, 0x03, 0x00, 0x00, 0x00, 0xa8, 0x1f, 0x00, 0x00, 0x00, 0x00, 0x3e, 
	0x00, 0x00, 0x00, 0x00, 0x2c, 0x00, 0x00, 0x00, 0x00, 0x70, 0x00, 0x18, 
	0x00, 0x00, 0x60, 0x40, 0x0c, 0x00, 0x00, 0xa0, 0x00, 0x08, 0x40, 0x00, 
	0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x05, 0x00, 0x80, 0x00, 0x00, 0x0c, 
	0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00, 0x30, 0x00, 0x00, 
	0x00, 0x00, 0x40, 0x00, 0x10, 0x00, 0x00, 0x00};
    static XSizeHints default_size_hints = {
      USPosition | USSize | PPosition | PSize, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    static int        status = 0, num = 0;
    XSetWindowAttributes   attrib;
    unsigned     long      valuemask;
    Window root;      unsigned int border_width;

    if ( !status ) {
      default_icon_pixmap
	  = (Pixmap) XCreateBitmapFromData(display, 
					   RootWindow(display, screen), 
					   (char*) default_icon_bitmap_bits, 
					   icon_bitmap_width, 
					   icon_bitmap_height);
      sprintf(default_win_name, " Window (%d) ", num);
      sprintf(default_icon_name, " Icon (%d) ", num);
      status = 1;
    }
    else {
      if (win->type == GRAPH_WINDOW || win->type == DIALOG_WINDOW) {
      num++;
      sprintf(default_win_name, " Window (%d) ", num);
      sprintf(default_icon_name, " Icon (%d) ", num);
      }
    }

    switch(win->type) {
      case GRAPH_WINDOW:
      case DIALOG_WINDOW:
        win->parent = RootWindow(win->display, win->screen);
	break;
      case BUTTON_WINDOW:
      case INPUT_WINDOW:
	break;
      default:
        win->parent = RootWindow(win->display, win->screen);
	break;
    }

    if (colormap != DefaultColormap(win->display, win->screen)) {
      attrib.colormap = colormap;
      attrib.background_pixel=WhitePixel(display, screen);
      attrib.border_pixel = BlackPixel(display, screen);
      valuemask = CWColormap | CWBackPixel | CWBorderPixel;
      win->window = XCreateWindow(win->display, win->parent, 
			  win->x, win->y, win->width, win->height, 1, 
			  win->window_depth, InputOutput, visual, valuemask, 
		          &attrib);
      win->colormap = colormap;
    }
    else {	     
      win->window = XCreateSimpleWindow(
		    win->display, win->parent, 
		    win->x, win->y, win->width, win->height, 1, 
                    BlackPixel(win->display, win->screen), 
		    WhitePixel(win->display, win->screen));
      win->colormap = DefaultColormap(win->display, win->screen);
    }

    if ( win->window_name == (char *)NULL ) win->window_name = default_win_name;
    if ( win->icon_name   == (char *)NULL ) win->icon_name   = default_icon_name;
    if ( win->icon_pixmap == (Pixmap)NULL ) win->icon_pixmap = default_icon_pixmap;


    if ( win->size_hints.flags == (long)NULL ) {
      win->size_hints.flags = default_size_hints.flags;
      win->size_hints.x = win->x;
      win->size_hints.y = win->y;
      win->size_hints.width  = win->width;
      win->size_hints.height = win->height;
      win->size_hints.min_width  = 0;
      win->size_hints.min_height = 0;
      switch(win->type) {
        case DIALOG_WINDOW:
        case BUTTON_WINDOW:
        case INPUT_WINDOW:
	  win->size_hints.flags = PAllHints | USPosition | USSize;
	  win->size_hints.min_width  = win->size_hints.max_width  = win->width;
	  win->size_hints.min_height = win->size_hints.max_height = win->height;
	  win->size_hints.width_inc  = win->size_hints.height_inc = 0;
	  break;
        case GRAPH_WINDOW:
        default:
	  break;
	}
    }

    lux_get_wingc(win);

    XSetStandardProperties(win->display, win->window, win->window_name, 
			   win->icon_name, win->icon_pixmap, 
			   (char **)NULL, 0, &win->size_hints);

    switch(win->type) {
      case GRAPH_WINDOW:
      case DIALOG_WINDOW:
        valuemask = ExposureMask | KeyPressMask | ButtonPressMask | 
	            ButtonReleaseMask | StructureNotifyMask |
		    ButtonMotionMask | EnterWindowMask | LeaveWindowMask 
		    /*| FocusChangeMask*/;
	break;
      case BUTTON_WINDOW:
	valuemask = ExposureMask | ButtonPressMask | 
	            EnterWindowMask | LeaveWindowMask;
	break;
      case INPUT_WINDOW:
	valuemask = ExposureMask | KeyPressMask | 
	            EnterWindowMask | LeaveWindowMask;
	break;
      default:
        valuemask = ExposureMask | KeyPressMask | ButtonPressMask | 
	            ButtonReleaseMask | StructureNotifyMask |
		    ButtonMotionMask | EnterWindowMask | LeaveWindowMask 
		    /*| FocusChangeMask*/;
	break;
      }
    XSelectInput(win->display, win->window, valuemask);

    XSetWMProtocols(win->display, win->window, &protocols, 1);

    switch(win->type) {
      case BUTTON_WINDOW:
      case INPUT_WINDOW:
        XMapWindow(win->display, win->window);
	break;
      case GRAPH_WINDOW:

	XMapWindow(win->display, win->window);
	XFlush(win->display);

	lux_wait_windowexpose(win->window);

	/*
	 * This is for reconfiguring the window in case it was changed by
	 * the WM.  The x, y will not be the values you provide, and even
	 * not relative to the root, but to its parent which WM creates
	 * for you
	 */

	XGetGeometry(win->display, win->window, &root, &win->x, &win->y, 
		     &win->width, &win->height, &border_width,
		     &win->window_depth);

	if (win->width != win->old_width) 
	  win->xresizefactor = (float)win->width/(float)win->user_width;
	if (win->height != win->old_height) 
	  win->yresizefactor = (float)win->height/(float)win->user_height;


	win->pixmap = XCreatePixmap(win->display, 
				    win->window, 
				    win->width, win->height, 
				    win->window_depth);

	XSync(win->display, False);/* Wait for operation to finish in order
				      for error to be handled */
	if (!pmpallocflag) {win->pixmap = (Pixmap)NULL; pmpallocflag = 1;}

	if (win->pixmap)         /* Here is to clear the pixmap  */
	  lux_clear_pixmap(win->display, win->pixmap, win->gc, 
			   0, 0, win->width, win->height);
	break;
      case DIALOG_WINDOW:
/*        XMapWindow(win->display, win->window);*/
	break;
      default:
	fprintf(stderr, "No window type specified !");
	XMapWindow(win->display, win->window);
	break;
      }
}

lux_wins *lux_allocwin()
{

     register lux_wins *newwin;

     if ( windows == (lux_wins *)NULL ) {
       if ( (newwin = windows = (lux_wins *)malloc(sizeof(lux_wins))) == NULL )
	              return (lux_wins *)NULL;
       windows->next = windows;
     }
     else {
       if ( (newwin = (lux_wins *)malloc(sizeof(lux_wins))) == NULL ) 
	              return (lux_wins *)NULL;
       newwin->next  = windows->next;
       windows->next = newwin;
     }
     
     /* Setup default values for window */
     newwin->win.display = display;
     newwin->win.screen  = screen;
     newwin->win.update_flag  = 1;      /* default to no animation    */
     newwin->win.discard_flag    = 0;      /* default to save everything */
     newwin->win.pixmap  = (Pixmap)NULL;
     newwin->win.window_depth = default_depth;
     newwin->win.window_name  = (char *)NULL;
     newwin->win.icon_name    = (char *)NULL;
     newwin->win.icon_pixmap  = (Pixmap)NULL;
     newwin->win.size_hints.flags   =   (long)NULL;
     newwin->win.data = (lux_data *)NULL;
     newwin->win.currentdata = (lux_data *)NULL;
     newwin->win.window_size = 0;      /* set to ok */
     newwin->win.msg[0] = (char)NULL;
     if (colorfile == NULL) newwin->win.colormapfile = NULL;
     else {
	 newwin->win.colormapfile = (char *)malloc(strlen(colorfile)+1);
	 newwin->win.colormapfile[0] = 0;
	 strcat(newwin->win.colormapfile, colorfile);
     }
     newwin->win.lux_colormap = lux_colormap;

     return newwin;
}

Window lux_openwin(x, y, width, height)	  /* Open a new graphics window */
unsigned int x, y, width, height;
{
    register lux_wins *mywin;
    static   int       status = 0;

    if ( !status ) {
      if(lux_xinit()) return 0;
      status = 1;
    }

    if ((mywin = lux_allocwin()) == (lux_wins *)NULL ) return (Window)NULL;
    
    mywin->win.type = GRAPH_WINDOW;

    mywin->win.x = x;
    mywin->win.y = y;
    mywin->win.width  = width;
    mywin->win.height = height;
    mywin->win.old_width   = width;
    mywin->win.old_height  = height;
    mywin->win.user_width  = width;
    mywin->win.user_height = height;
    mywin->win.xresizefactor = 1.0;
    mywin->win.yresizefactor = 1.0;

    mywin->win.buttonflag = 0;

    mywin->win.fxorg = 0.0;
    mywin->win.fyorg = 0.0;
    mywin->win.fxsize = 10.0;
    mywin->win.fysize = 10.0;
    mywin->win.xfactor = 1.0;
    mywin->win.yfactor = 1.0;
    mywin->win.xmin = 1.0;
    mywin->win.xmax = width;
    mywin->win.ymin = 1.0;
    mywin->win.ymax = height;
    mywin->win.xsize = width;
    mywin->win.ysize = height;
    mywin->win.xorg = 0;
    mywin->win.yorg = 0;
    
    mywin->win.lnx = 0;    /* Default to ordinary axis */
    mywin->win.lny = 0;

    lux_createwin(&mywin->win);
 
/* Make sure that all this setup to default when start to redraw */   
    lux_setup_region(mywin->win.window, 
		     mywin->win.fxorg, mywin->win.fyorg, 
		     mywin->win.fxsize, mywin->win.fysize);
    lux_setup_axis_style(mywin->win.window, mywin->win.lnx, mywin->win.lny);
    lux_setup_axis(mywin->win.window, 
		   mywin->win.xmin, mywin->win.xmax, 
		   mywin->win.ymin, mywin->win.ymax);
/* Maybe we should add clear screen to set the window clear to default */

    return mywin->win.window;
}

lux_wins *get_currentwin(window)
Window    window;
{
    register lux_wins *start;
    
    if ( windows == (lux_wins *)NULL ) 
      fprintf(stderr, "lux_err: NULL window!\n");
    start = windows;
    while (start->win.window != window ) { 
      start = start->next;
      if (start == windows && start->win.window != window) {
	fprintf(stderr, "lux_err: Window id = %d not found!\n", window);
	break;
      }
    }
    windows = start;
    return start;    
}

lux_freedata(win)
Window win;
{
    register lux_wins *current;
    register int      *ptr;

    current = get_currentwin(win);

    if ( current->win.data == (lux_data *)NULL ) return;
    
    while(1) {
      current->win.currentdata = current->win.data;
      if (current->win.data == (current->win.data)->next) {
	if((current->win.currentdata)->data.b != (char *)NULL) {
	    if ((current->win.currentdata)->type == DRAW_VSTRING) {
		ptr = (int *)(&(current->win.currentdata)->data.f[3]);
		XDestroyImage(((XImage **)(&ptr[2]))[0]);
	    }
	    free((current->win.currentdata)->data.b);
	}
	free(current->win.currentdata);
	break;
      }

      current->win.data = (current->win.data)->next;
      if((current->win.currentdata)->data.b != (char *)NULL) {
	  if ((current->win.currentdata)->type == DRAW_VSTRING) {
	      ptr = (int *)(&(current->win.currentdata)->data.f[3]);
	      XDestroyImage(((XImage **)(&ptr[2]))[0]);
	  }
	  free((current->win.currentdata)->data.b);
      }
      free(current->win.currentdata);
    }

    current->win.currentdata = current->win.data = (lux_data *)NULL;
}

lux_freewins()
{
    register lux_wins *current, *start, *tmp;
    XEvent             report;

    if ( windows == (lux_wins *)NULL ) return;

    current = start = windows;

    XDestroyWindow(current->win.display, popup);
    lux_create_popup();

    do {
      lux_freedata(current->win.window);
      current = current->next;
    } while(current != start);
  
    do {
      tmp = current;
      current = current->next;
      if ((tmp->win.icon_pixmap != default_icon_pixmap) && 
	  tmp->win.icon_pixmap) {
	XFreePixmap(tmp->win.display, tmp->win.icon_pixmap);
	tmp->win.icon_pixmap = (Pixmap)NULL;
      }
      if (tmp->win.pixmap) {
	XFreePixmap(tmp->win.display, tmp->win.pixmap);
	tmp->win.pixmap = (Pixmap)NULL;
      }
      XUnmapWindow(tmp->win.display, tmp->win.window);
      XDestroyWindow(tmp->win.display, tmp->win.window);
      XSync(tmp->win.display, False);
      while(XCheckWindowEvent(tmp->win.display, 
			      tmp->win.window, 
			      0xFFFFFFFF, &report));
      XFreeGC(tmp->win.display, tmp->win.gc);
      free(tmp);
    } while(current != start);

    windows = (lux_wins *)NULL;
}

lux_reset_window(win)    /* Clear everything and the default values are random*/
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);
  
    XClearWindow(current->win.display, current->win.window);
    if (current->win.pixmap) 
      lux_clear_pixmap(current->win.display, current->win.pixmap, 
		       defaultgc, 0, 0, 
		       current->win.width, current->win.height);
    XFlush(current->win.display);
    lux_freedata(current->win.window);
}

lux_freewin(win)	/* Delete the specified window */
Window win;
{
    register lux_wins *current, *tmp;
    XEvent             report;

    if ( windows == (lux_wins *)NULL ) return;

    current = get_currentwin(win);

    XDestroyWindow(current->win.display, popup);
    lux_create_popup();

    lux_freedata(current->win.window);

    tmp = current;
    if ((tmp->win.icon_pixmap != default_icon_pixmap) && 
	tmp->win.icon_pixmap) {
      XFreePixmap(tmp->win.display, tmp->win.icon_pixmap);
      tmp->win.icon_pixmap = (Pixmap)NULL;
    }
    if (tmp->win.pixmap) {
      XFreePixmap(tmp->win.display, tmp->win.pixmap);
      tmp->win.pixmap = (Pixmap)NULL;
    }
    XUnmapWindow(tmp->win.display, tmp->win.window);
    XDestroyWindow(tmp->win.display, tmp->win.window);
    XSync(tmp->win.display, False);
    while(XCheckWindowEvent(tmp->win.display, 
			    tmp->win.window, 
			    0xFFFFFFFF, &report));
    XFreeGC(tmp->win.display, tmp->win.gc);
    while(tmp->next != current) tmp = tmp->next;
    tmp->next = current->next;
    if (current == current->next)
      windows = (lux_wins*)NULL;
    else  windows = current->next;
    free(current);
}

lux_exit()
{
  if (display == (Display *)NULL) {
    fprintf(stderr, "No open window!\n");
    exit(0);
  }

  lux_freewins();
  XFreePixmap(display, default_icon_pixmap);
  XUnloadFont(display, font_info->fid);
  XFreeGC(display, defaultgc);
  XFreeGC(display, redrawgc);
  XCloseDisplay(display);
  exit(1);
}

lux_quick_exit()
{
  if (display == (Display *)NULL) {
/*    fprintf(stderr, "No open window!\n");*/
   if (mcdrawflag) return 0;
   else exit(0);
  }

  XFreePixmap(display, default_icon_pixmap);
  XUnloadFont(display, font_info->fid);
  XFreeGC(display, defaultgc);
  XFreeGC(display, redrawgc);
  XCloseDisplay(display);
  if (mcdrawflag) return 0;
  else exit(1);
}

static int wait_time = 0;
void set_wait_time(i)
int i;
{
    wait_time = i;
}
int get_wait_time()
{
    return wait_time;
}

char win_getkey(Window win)

/*
 * See if a key has been pressed in the specified window, and respond
 * to a resize event in any window.  Return after the first event with
 * NULL or the key that was pressed.
 */

{
    XEvent              report;
    register lux_wins  *current;
    char                buffer[255];
    KeySym              keysym;
    int                 buffersize = 20;
    int                 charcount;
    XComposeStatus      compose;

    buffer[0] = buffer[1] = 0;

    if (windows == (lux_wins *)NULL) {
	if (REPORT_EVENT) fprintf(stderr, "No open window!\n");
	return 0;
    }

    /* Grab the next event from the queue. */

    XNextEvent(display, &report);

    current = get_currentwin(report.xany.window);

    /* (The current window ID, as returned by lux_openwin, is
     *
     *		current->win.window
     *
     *  which is of type Window, aka unsigned long.)
     */

    /* Expose and ConfigureNotify events must be acted upon
     * regardless of the window they occur in.
     *
     * DISCARD all other events, except a keypress in the
     * specified window.
     */

    if (report.type == Expose) {	/* Window exposed (or resized?) */

	if (REPORT_EVENT) fprintf(stderr, "Get Expose event in win %u\n", 
				  current->win.window);

	if (!report.xexpose.count) {

	    if (current->win.window_size == SMALL)
		/* do something for small? */
		;

	    redraw(current->win.window, 0);
	}

    } else if (report.type == ConfigureNotify) {	/* Window resized */

	if (REPORT_EVENT) fprintf(stderr,
				  "Get ConfigureNotify in win %u\n", 
				  current->win.window);

	current->win.x = report.xconfigure.x;
	current->win.y = report.xconfigure.y;
	current->win.width  = report.xconfigure.width;
	current->win.height = report.xconfigure.height;

	if (current->win.old_width  != current->win.width) 
	    current->win.xresizefactor
		= (float)current->win.width /
		    (float)current->win.user_width;

	if (current->win.old_height != current->win.height) 
	    current->win.yresizefactor
		= (float)current->win.height /
		    (float)current->win.user_height;

	if ((current->win.width
	     < current->win.size_hints.min_width) ||
	    (current->win.height
	     < current->win.size_hints.min_height))
	    current->win.window_size = SMALL;
	else
	    current->win.window_size = OK;

    } else if (current->win.window == win && report.type == KeyPress)

	if (XLookupString(&(report.xkey), 
			  buffer, buffersize, 
			  &keysym, &compose)) return buffer[0];

    return 0;
}

void win_checkevent(Window win)

/*
 * Respond to all pending expose or resize events in the specified
 * window.  Return when the pending stack is exhausted.
 */

{
    XEvent              report;
    register lux_wins  *current;
    char                buffer[255];
    KeySym              keysym;
    int                 buffersize = 20;
    int                 charcount;
    XComposeStatus      compose;

    buffer[0] = buffer[1] = 0;

    if (windows == (lux_wins *)NULL)
	if (REPORT_EVENT) fprintf(stderr, "No open window!\n");

    current = get_currentwin(win);

    /* (The current window ID, as returned by lux_openwin, is
     *
     *		current->win.window
     *
     *  which is of type Window, aka unsigned long.)
     */

    /* Loop through all pending events for this window. */

    while (XCheckWindowEvent(current->win.display, 
			    current->win.window, 
			    0xFFFFFFFF, &report)) {

	/* Act upon Expose and ConfigureNotify events.  DISCARD all others. */

	if (report.type == Expose) {	/* Window exposed (or resized?) */

	    if (REPORT_EVENT) fprintf(stderr, "Get Expose event in win %u\n", 
				      current->win.window);

	    if (!report.xexpose.count) {

		if (current->win.window_size == SMALL)
		    /* do something for small? */
		    ;

		redraw(current->win.window, 0);
	}

	} else if (report.type == ConfigureNotify) {	/* Window resized */

	    if (REPORT_EVENT) fprintf(stderr,
				      "Get ConfigureNotify in win %u\n", 
				      current->win.window);

	    current->win.x = report.xconfigure.x;
	    current->win.y = report.xconfigure.y;
	    current->win.width  = report.xconfigure.width;
	    current->win.height = report.xconfigure.height;

	    if (current->win.old_width  != current->win.width) 
		current->win.xresizefactor
		    = (float)current->win.width /
			(float)current->win.user_width;

	    if (current->win.old_height != current->win.height) 
		current->win.yresizefactor
		    = (float)current->win.height /
			(float)current->win.user_height;

	    if ((current->win.width
		 < current->win.size_hints.min_width) ||
		(current->win.height
		 < current->win.size_hints.min_height))
		current->win.window_size = SMALL;
	    else
		current->win.window_size = OK;

	}
    }
}

int win_getevent(Window win, int win_type, int hide_on_cr)

/* Expands on functionality of old lux_getevent.
 * Comments and modifications added 10/95, SLWM.
 */

{
    XEvent              report;
    register lux_wins  *current;
    char                buffer[255];
    KeySym              keysym;
    int                 buffersize = 20;
    int                 charcount;
    XComposeStatus      compose;
    Time                buttondown;
    unsigned int        button, xbutton, ybutton;
    int                 click = 0;
    Window              wbutton;
    int                 pop = 0;

    buffer[0] = buffer[1] = 0;

    if (windows == (lux_wins *)NULL) {
	if (REPORT_EVENT) fprintf(stderr, "No open window!\n");
	return 0;
    }

    /* Enter infinite event loop. */

    while (1)  {

	if (wait_time > 0) lux_pause(wait_time);  /* Unit = microseconds */

	if (mcdrawflag) {
	    if (XPending(display) == 0) {
		if (get_timeout()) return 7;
		continue;
	    }
	}

	/* Grab the next event from the queue. */

	XNextEvent(display, &report);

	current = get_currentwin(report.xany.window);

	/* (The current window ID, as returned by lux_openwin, is
	 *
	 *		current->win.window
	 *
	 *  which is of type Window, aka unsigned long.)
	 */

/*
	fprintf(stderr,
                "win_getevent: target win = %ld, type = %d\n",
		win, win_type, current->win.window);
	fprintf(stderr,
	        "              event win = %ld, type = %d, parent = %ld\n",
		current->win.type, current->win.parent);
*/

	/* Note: Expose and ConfigureNotify events must be acted upon
	 * regardless of the window they occur in.  Others may be
	 * limited depending on the input preferences.
	 */

	if (report.type == Expose) {	/* Window exposed (or resized?) */

	    if (REPORT_EVENT) fprintf(stderr, "Get Expose event in win %u\n", 
				      current->win.window);

	    if (!report.xexpose.count) {

		if (current->win.window_size == SMALL)
		    /* do something for small? */
		    ;

		redraw(current->win.window, 0);
	    }

	} else if (report.type == ConfigureNotify) {	/* Window resized */

	    if (REPORT_EVENT) fprintf(stderr,
				      "Get ConfigureNotify in win %u\n", 
				      current->win.window);

	    current->win.x = report.xconfigure.x;
	    current->win.y = report.xconfigure.y;
	    current->win.width  = report.xconfigure.width;
	    current->win.height = report.xconfigure.height;

	    if (current->win.old_width  != current->win.width) 
		current->win.xresizefactor
		    = (float)current->win.width /
			(float)current->win.user_width;

	    if (current->win.old_height != current->win.height) 
		current->win.yresizefactor
		    = (float)current->win.height /
			(float)current->win.user_height;

	    if ((current->win.width
		 < current->win.size_hints.min_width) ||
		(current->win.height
		 < current->win.size_hints.min_height))
		current->win.window_size = SMALL;
	    else
		current->win.window_size = OK;

	} else if ( (win <= 0 && win_type <= 0)
		     || current->win.window == win
		     || current->win.parent == win
		     || current->win.type == win_type) {

	    /* Other events may be restricted in scope. */

	    switch (report.type) {		/* Branch on event type */

		case ButtonPress:		/* A button has been pressed */

		  if (REPORT_EVENT) {
		      fprintf(stderr, "get button press in win %u\n", 
			      current->win.window);
		      fprintf(stderr, 
			      "Button %d, position x= %d, y= %d, time= %u\n", 
			      report.xbutton.button, 
			      report.xbutton.x, 
			      report.xbutton.y, 
			      report.xbutton.time);
		  }

		  switch(current->win.type) {	/* Action depends on window */

		    case DIALOG_WINDOW:
		    case GRAPH_WINDOW:

		      if (report.xbutton.button == Button1
			   || report.xbutton.button == Button3) {

			  /* Left or right mouse button pressed: pop up a
			   * window showing the current cursor position. */

			  XReparentWindow(current->win.display, 
					  popup, 
					  current->win.window, 
					  0, 
					  (report.xbutton.y
					    > current->win.height/2) ?
					  0 : (current->win.height
					       - (font_info->ascent
						  + font_info->descent) - 4) );
			  XMapWindow(current->win.display, popup);
			  pop = 1;
			  XFlush(current->win.display);

			  {
			      float fx, fy;
			      lux_reconvert_coord(current->win.window, 
						  report.xmotion.x, 
						  report.xmotion.y, 
						  &fx, &fy);
			      sprintf(buffer, "(%d, %d) (%f, %f)        ", 
				      report.xmotion.x, report.xmotion.y, 
				      fx, fy);
			      XDrawImageString(display, popup, defaultgc, 
					       0, font_info->ascent, 
					       buffer, strlen(buffer));
			      XFlush(current->win.display);
			  }

			  XDefineCursor(current->win.display, 
					current->win.window, cursor);
			  XFlush(current->win.display);

			  /* Special treatment of button 1: */

			  if (report.xbutton.button == Button1) {
			      if (wbutton != report.xbutton.window) click = 0;
			      if (report.xbutton.time - buttondown
				    >= DOUBLECLICKLENGTH)
				  click = 0; 
			      if (!click) {
				  buttondown = report.xbutton.time;
				  wbutton = report.xbutton.window;
			      }
			  }
		      }
		      break;

		    case BUTTON_WINDOW:

		      /* Return a value of 3 or 4 if an "exit"
		       * button has been pressed. */

		      if (current->win.subtype == OK_BUTTON) {

			  lux_ok_data(current->win.parent);
			  lux_hide_dialog(current->win.parent);

			  return 3;   /* OK data */

		      } else if (current->win.subtype == OK_KEEP_BUTTON) {

			  lux_ok_data(current->win.parent);

			  return 3;   /* OK data */

		      } else if (current->win.subtype == CANCEL_BUTTON) {

			  lux_cancel_data(current->win.parent);
			  lux_hide_dialog(current->win.parent);

			  return 4;   /* Cancel data */

		      } else if (current->win.subtype == CHECK_BUTTON) {

			  /* Otherwise, modify a check button. */

			  if (current->win.msg[0] == 1) {

			      /* Erase the box contents. */

			      lux_reset_window(current->win.window);
			      current->win.msg[0] = 0;

			  } else if (current->win.msg[0] == 0) {

			      /* Place an "X" in the box. */

			      lux_draw_line(current->win.window, 0, 0, 
					    current->win.width, 
					    current->win.height);
			      lux_draw_line(current->win.window, 0, 
					    current->win.height, 
					    current->win.width, 0);
			      current->win.msg[0] = 1;
			  }
		      }  
		      break;

		    case INPUT_WINDOW:
		      default:
		      break;
		  }
		  break;

		case ButtonRelease:		/* A button was released */

		  if (REPORT_EVENT) {
		      fprintf(stderr, "get button release in win %u\n", 
			      current->win.window);
		      fprintf(stderr, 
			      "Button %d, position x= %d, y= %d, time= %u\n", 
			      report.xbutton.button, 
			      report.xbutton.x, 
			      report.xbutton.y, 
			      report.xbutton.time);
		  }

		  if (report.xbutton.button == Button1) {

		      /* Left mouse button: remove popup and deal
		       * with double click. */

		      XUnmapWindow(display, popup);
		      pop = 0;

		      if (report.xbutton.time - buttondown
			    < DOUBLECLICKLENGTH) {

			  if (click) {

			      if (REPORT_EVENT) fprintf(stderr,
							"Double Click\n");

			      lux_reconvert_coord(current->win.window, 
						  report.xbutton.x, 
						  report.xbutton.y, 
						  &current->win.xbutton, 
						  &current->win.ybutton);
			      if (REPORT_EVENT)
				  fprintf(stderr, 
					  "Button position %f %f\n", 
					  current->win.xbutton, 
					  current->win.ybutton);

			      XBell(current->win.display, 100);
			      click = 0;
			      XUndefineCursor(current->win.display, 
					      current->win.window);

			      if (current->win.buttonflag) return 0;

			  } else {

			      if (REPORT_EVENT) fprintf(stderr, "One Click\n");

			      click = 1;
			      XUndefineCursor(current->win.display, 
					      current->win.window);
			  }

		      } else {

			  click = 0;
			  XUndefineCursor(current->win.display, 
					  current->win.window);
		      }

		  } else if (report.xbutton.button == Button3) {

		      /* Right mouse button: just remove popup */

		      XUnmapWindow(display, popup);
		      pop = 0;
		      lux_reconvert_coord(current->win.window, 
					  report.xbutton.x, report.xbutton.y, 
					  &current->win.xbutton, 
					  &current->win.ybutton);
		      XUndefineCursor(current->win.display, 
				      current->win.window);

		      if (current->win.buttonflag) return 0;

		  }
		  break;

		case MotionNotify:		/* The cursor has moved */

		  if (pop) {

		      /* Update the popup display. */

		      {
			  float fx, fy;
			  lux_reconvert_coord(current->win.window, 
					      report.xmotion.x, 
					      report.xmotion.y, 
					      &fx, &fy);	
			  sprintf(buffer, 
			  "(%d, %d) (%g, %g)                                  ", 
				  report.xmotion.x, report.xmotion.y, fx, fy);
			  XDrawImageString(display, popup, defaultgc, 
					   0, font_info->ascent, 
					   buffer, strlen(buffer));
		      }
		      XFlush(current->win.display);
		  } 
		  break;

		case KeyPress:			/* Keyboard input */

		  switch(current->win.type) {	/* Action depends on window */

		    case DIALOG_WINDOW:	/* Dialog, outside of input regions */

		      charcount = XLookupString(&(report.xkey), 
						buffer, buffersize, 
						&keysym, &compose);

		      if (   ((keysym >= XK_KP_Space) && (keysym <= XK_KP_9))
			  || ((keysym >= XK_space)
			        && (keysym <= XK_asciitilde))
			  || (keysym == XK_Return
			      || keysym == XK_Linefeed
			      || keysym == XK_KP_Enter) ) {

			  /* Keypad ' ' - '9', <CR>, or nonwhite
			     ASCII character. */

			  /* Terminate the input string. */

			  buffer[charcount] = 0;   /* A little bit stupid */
			  if (charcount == 0) break;

			  /* o/O in dialog window ==> OK / OK and hide
			   * 				(<CR> <==> 'o' or 'O')
			   * c/C in dialog window ==> cancel / cancel and hide
			   * q/Q in dialog window ==> quit program
			   */

			  if (buffer[0] == 'o' || buffer[0] == 'O'
			      || keysym == XK_Return
			      || keysym == XK_Linefeed
			      || keysym == XK_KP_Enter) {

			      /* Return, use dialog data, optionally
			       * hide window. */

			      lux_ok_data(current->win.window);
			      if (buffer[0] == 'O'
				  || (buffer[0] != 'o' && hide_on_cr))
				  lux_hide_dialog(current->win.window);

			      return 3;   /* OK data */

			  } else if (buffer[0] == 'c' || buffer[0] == 'C') {

			      /* Return, cancel dialog data, optionally
			       * hide window. */

			      lux_cancel_data(current->win.window);
			      if (buffer[0] == 'C')
				  lux_hide_dialog(current->win.window);

			      return 4;   /* Cancel data */

			  } else if (buffer[0] == 'q' || buffer[0] == 'Q') {

			      /* Quit program. */

			      lux_quick_exit();
			  }
		      }
		      break;

		    case GRAPH_WINDOW:

		      if (mcdrawflag) {

			  charcount = XLookupString(&(report.xkey), 
						    buffer, buffersize, 
						    &keysym, &compose);
			  if (REPORT_EVENT)
			      fprintf(stderr, 
				      "get %d %d %s keypress in win %u\n", 
				      charcount, buffer[0], 
				      buffer, current->win.window);

			  if (charcount == 0) break;

			  if (keysym == XK_Return
			      || keysym == XK_Linefeed
			      || keysym == XK_KP_Enter)
			      *buffer = '\n';

			  if (keysym == XK_BackSpace || keysym == XK_Delete)
			      *buffer = '\b';

			  if (process_char(*buffer)) return 6;

		      } else {

			  /* Allow limited input via the graph window. */

			  charcount = XLookupString(&(report.xkey), 
						    buffer, buffersize, 
						    &keysym, &compose);

			  if (((keysym >= XK_KP_Space)
			         && (keysym <= XK_KP_9))
			      || ((keysym >= XK_space)
				    && (keysym <= XK_asciitilde))) {

			      /* Keypad ' ' - '9' or nonwhite character. */

			      /* Terminate the input string. */

			      buffer[charcount] = 0;  /* A little bit stupid */

			      if (REPORT_EVENT)
				  fprintf(stderr, 
					  "get %d %s keypress in win %u\n", 
					  charcount, buffer, 
					  current->win.window);

			      if (charcount == 0) break;

			      /* q ==> quit
			       * d ==> show dialog (what if more than 1?)
			       * c ==> return, leave dialog unread
			       * C ==> clear window
			       * r ==> redraw window
			       */

			      if (buffer[0] == 'q')

				  lux_quick_exit();

			      else if (buffer[0] == 'd')

				  lux_show_dialog(current->win.window);

			      else if (buffer[0] == 'c') {

				  /* A 'c' in the graph window will cause an
				     immediate return.  Any dialog windows
				     will remain open, but will not be
				     read. */

				  return 1;

			      } else if (buffer[0] == 'C')

				  lux_reset_window(current->win.window);

			      else if (buffer[0] == 'r')

				  redraw(current->win.window, 1);
			  }
		      }

		      break;

		    case INPUT_WINDOW:		/* Modify an input window */

		      charcount = XLookupString(&(report.xkey), 
						buffer, buffersize, 
						&keysym, &compose);
		      if (charcount == 0) break;

		      if (REPORT_EVENT)
			  fprintf(stderr, 
			      "get %d %s keypress in window %u, keysym %x\n", 
				  charcount, buffer, current->win.window,
				  keysym);

		      /* Terminate the string. */

		      buffer[charcount] = 0; 

		      {
			  char *tmp;
			  tmp = &(current->win.data->data.b[3*sizeof(float)
							     + 1]); /* wow! */

			  /* Note: tmp is a pointer to the first character of
			     the string associated with this input window. */

			  if ((keysym == XK_Return)
			      || (keysym == XK_KP_Enter)
			      || (keysym == XK_Linefeed))  {

			      lux_ok_data(current->win.parent);

			      if (hide_on_cr)
				  lux_hide_dialog(current->win.parent);

			      return 3;

			  } else if (((keysym >= XK_KP_Space)
				      && (keysym <= XK_KP_9))
				   || ((keysym >= XK_space)
				         && (keysym <= XK_asciitilde))) {

			      /* Keypad ' ' - '9' or nonwhite character. */

			      /* Build up the string. */

			      if ((int)(strlen(tmp) + strlen(buffer)) > 230) 
				  XBell(current->win.display, 100);
			      else
				  strcat(tmp, buffer);

			  } else if ((keysym == XK_BackSpace)
				       || (keysym == XK_Delete)) {

			      int length;
			      if ((length = strlen(tmp)) > 0)
				  tmp[length-1] = (char)NULL;

			  } else if (keysym == XK_Escape) {

			      /* ESCAPE deletes the entire string. */

			      tmp[0] = (char)NULL;

			  }
		      }
		      redraw(current->win.window, 1);
		      break;

		    case BUTTON_WINDOW:
		      default:
		      break;
		  }
		  buffer[1]=(char)NULL;
		  break;

		case EnterNotify:		/* Highlight input or button
						   window on entry */

		  if (REPORT_EVENT) fprintf(stderr, "EnterNotify in win %u\n", 
					    current->win.window);

		  /* XSetInputFocus(current->win.display, focus, 
				    RevertToPointerRoot, CurrentTime); */

		  if (current->win.type == INPUT_WINDOW || 
		      current->win.type == BUTTON_WINDOW) 
		      lux_highlight(current->win.window);

		  break;

		case LeaveNotify:		/* Unhighlight input or button
						   window on exit */

		  if (REPORT_EVENT) fprintf(stderr, "LeaveNotify in win %u\n", 
					    current->win.window);

		  /* XSetInputFocus(current->win.display, PointerRoot, 
				    RevertToPointerRoot, CurrentTime); */

		  if (current->win.type == INPUT_WINDOW || 
		      current->win.type == BUTTON_WINDOW) 
		      lux_unhighlight(current->win.window);

		  break;

		case FocusIn:

		  if (REPORT_EVENT) fprintf(stderr, "FocusIn in win %u\n", 
					    current->win.window);

		  /* XSetInputFocus(current->win.display, focus, 
				    RevertToPointerRoot, CurrentTime); */

		  break;

		case FocusOut:

		  if (REPORT_EVENT) fprintf(stderr, "FocusOut in win %u\n", 
					    current->win.window);

		  break;

		case DestroyNotify:	  /* Window has been deleted via WM */

		  if (REPORT_EVENT) fprintf(stderr, "Destroy win %u\n", 
					    current->win.window);

		  lux_freewin(current->win.window);
		  break;

		case ClientMessage:

		  if (REPORT_EVENT) fprintf(stderr, 
					    "Got a ClientMessage in win %u\n", 
				      current->win.window);

		  if (current->win.type != GRAPH_WINDOW)

		      XBell(current->win.display, 100);

		  else {

		      lux_wins *tmp; tmp = current->next;
		      while(tmp->win.type != GRAPH_WINDOW) tmp = tmp->next;

		      if (tmp == current) {

			  XBell(current->win.display, 100);
			  fprintf(stderr, 
				  "Ambiguous command! Type q to quit..\n");

		      } else if (report.xclient.message_type == 
				 XInternAtom(display, "WM_PROTOCOLS", True)) {

/*
			  XDestroyWindow(current->win.display, popup);
			  lux_create_popup();
			  XUnmapWindow(current->win.display, 
				       current->win.window);
			  XDestroyWindow(current->win.display, 
					 current->win.window);
			  XSync(current->win.display, False);
			  while(XCheckWindowEvent(current->win.display, 
						  current->win.window, 
						  0xFFFFFFFF, &report));
*/
			  lux_freewin(current->win.window);
		      }
		  }
		  break;

		case MappingNotify:

		  if (REPORT_EVENT) fprintf(stderr, 
					    "Got a MappingNotify in win %u\n", 
					    current->win.window);

		  /* XRefreshKeyboardMapping(&report); */
		  break;

		default:
		  break;

	    }	/* switch (report.type)	*/

	}	/* else if (win...	*/

    }		/* while(1)		*/
}

int lux_getevent()
{
    return win_getevent(0, 0, 0);
}

get_mouse_position(win, fx, fy)
Window win;
float *fx, *fy;
{
  register lux_wins *current;

  current = get_currentwin(win);
  
  current->win.buttonflag = 1;

  lux_getevent();

  current->win.buttonflag = 0;

  *fx = current->win.xbutton;
  *fy = current->win.ybutton;

}

lux_set_window_name(win, name)
Window win;
char  *name;
{
  register lux_wins *current;

  current = get_currentwin(win);

  XStoreName(current->win.display, current->win.window, name);
  XSetIconName(current->win.display, current->win.window, name);

}


lux_display(width, height)
unsigned int *width, *height;
{
    Display *display;

    if ( (display=XOpenDisplay(display_name)) == NULL ) {
	(void) fprintf( stderr, 
		       "Cannot connect to X server %s\n", 
		       XDisplayName(display_name));
	if (mcdrawflag) return 1;
	else exit( -1 );
    }
 
    screen = DefaultScreen(display);
    *width  = DisplayWidth(display, screen);
    *height = DisplayHeight(display, screen);
    XCloseDisplay(display);
    return 0;
}


lux_use_mcdraw()
{
    mcdrawflag = 1;
}

lux_iconify_window(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    XIconifyWindow(current->win.display, current->win.window, current->win.screen);
}


lux_uniconify_window(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

/*    XUnmapWindow(current->win.display, current->win.window);*/
    XMapWindow(current->win.display, current->win.window);
}

lux_raise_window(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    XRaiseWindow(current->win.display, current->win.window);
}

lux_clear_keyboard_buffer()
{
    XEvent report;
    while(XCheckTypedEvent(display, KeyPress, &report));
}
