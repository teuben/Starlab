/* draw.c
 * Biao Lu                     biao@eagle.drexel.edu
 */

#include "win.h"

extern GC        redrawgc, defaultgc;
extern int       pmpallocflag;
extern XFontStruct  *font_info;
extern lux_wins *get_currentwin();
extern unsigned  long lux_get_bgcolor();

#define RESIZED              1
#define ONE_MORE_EXPOSE      2
#define ABORT_REDRAW         3


#define DEBUG 0


Bool lux_expose_pred(display, event, msg)
Display  *display;
XEvent   *event;
char     *msg;
{
    register Window *win;

    win = (Window *)msg;

    switch (event->type) {
        case ConfigureNotify:
          if (win[0] == event->xany.window) win[1] = 1;
	  break;
        case Expose:
          if (win[0] == event->xany.window) win[2] = 2;
	  break;
	case ButtonPress:   
	  if ((win[0] == event->xany.window) && 
	      (event->xbutton.button == Button2)) {
	    win[1] = win[2] = 3;
	    return(True);
	  }
	  break;
	default:
	  break;
    }
    return(False);
}

lux_peek_expose(display,win)
Display *display;
Window   win;
{
    XEvent   report;
    Window  *w;
    char     msg[3*sizeof(Window)];
    
    w = (Window *)msg;
    w[0] = win; 

    w[1] = w[2] = 0;
    XCheckIfEvent(display,&report,lux_expose_pred,msg);

    if (w[2] == 2) {
/*      fprintf(stderr,"lux_msg: One more expose in win %u\n", win);*/
      return(ONE_MORE_EXPOSE);
    }
    if (w[1] == 3 || w[2] == 3) {
/*      fprintf(stderr,"lux_msg: Abort in redraw in win %u\n", win);*/
      return(ABORT_REDRAW);
    }

    return 0;
}

lux_peek_expose_config(display,win)
Display *display;
Window   win;
{
    XEvent   report;
    Window  *w;
    char     msg[3*sizeof(Window)];
    
    w = (Window *)msg;
    w[0] = win; 

    w[1] = w[2] = 0;
    XCheckIfEvent(display,&report,lux_expose_pred,msg);

    if (w[1] == 1) {
/*      fprintf(stderr,"lux_msg: Size changed in win %u\n", win);*/
      return(RESIZED);
    }
    if (w[2] == 2) {
/*      fprintf(stderr,"lux_msg: One more expose in win %u\n", win);*/
      return(ONE_MORE_EXPOSE);
    }
    if (w[1] == 3 || w[2] == 3) {
/*      fprintf(stderr,"lux_msg: Abort in redraw in win %u\n", win);*/
      return(ABORT_REDRAW);
    }

    return 0;
}

lux_data *get_newdata(wins)
lux_wins *wins;
{
    register lux_data *newdata, *current;

    if (wins->win.data == NULL) {

      if ( (wins->win.data = (lux_data *)malloc(sizeof(lux_data))) == 
	    (lux_data *)NULL ) return((lux_data *)NULL);

      (wins->win.data)->next = wins->win.data;
      wins->win.currentdata = wins->win.data;
      wins->win.data->data.d = (double *)NULL;
      newdata = wins->win.data;

    } else {

      newdata = (lux_data *)malloc(sizeof(lux_data));

      if ( newdata == (lux_data *)NULL )
	  return((lux_data *)NULL);

      newdata->data.d = (double *)NULL;
      current = wins->win.currentdata;
      while(current->next != current) current = current->next;
      current->next = newdata;
      newdata->next = newdata;
      wins->win.currentdata = newdata;

    }

    return(newdata);
}

lux_clear_window(win)    /* Just clear the graph not data */
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);
  
    if (current->win.pixmap) 
      lux_clear_pixmap(current->win.display,current->win.pixmap,
		       current->win.gc,0,0,
		       current->win.width,current->win.height);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      lux_clear_pixmap(current->win.display,current->win.window,
		       current->win.gc,0,0,
		       current->win.width,current->win.height);
      XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't clear window -- not enough memory\n"); 
	return 0;}

    newdata->type = CLEAR_WINDOW;


}

lux_reclear_window(win)    /* Just clear the graph not data */
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);
  
    if (current->win.pixmap) 
      lux_clear_pixmap(current->win.display,current->win.pixmap,
		       redrawgc,0,0,
		       current->win.width,current->win.height);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      lux_clear_pixmap(current->win.display,current->win.window,
		       redrawgc,0,0,
		       current->win.width,current->win.height);
      XFlush(current->win.display);
    }

}

lux_clear_current_region(win)
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;
    int x,y,width,height;

    current = get_currentwin(win);

    x = current->win.xorg*current->win.xresizefactor+0.5;
    y = current->win.yorg*current->win.yresizefactor+0.5;
    width  = current->win.xsize*current->win.xresizefactor+0.5;
    height = current->win.ysize*current->win.yresizefactor+0.5;

    if (current->win.pixmap) {
      lux_clear_pixmap(current->win.display, current->win.pixmap, 
		       current->win.gc, 
		       x+1, y+1, width-2, height-2);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  
      lux_clear_pixmap(current->win.display, current->win.window, 
		       current->win.gc, 
		       x+1, y+1, width-2, height-2);

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't clear display -- not enough memory\n"); 
	return 0;}

    newdata->type = CLEAR_CURRENT_REGION;

    return 1;
}

lux_reclear_current_region(win)
Window win;
{
    register lux_wins *current;
    int x,y,width,height;

    current = get_currentwin(win);

    x = current->win.xorg*current->win.xresizefactor+0.5;
    y = current->win.yorg*current->win.yresizefactor+0.5;
    width  = current->win.xsize*current->win.xresizefactor+0.5;
    height = current->win.ysize*current->win.yresizefactor+0.5;

    if (current->win.pixmap) {
/*
	XCopyArea(current->win.display, current->win.pixmap,
		  current->win.window,
		  redrawgc, 
		  0, 0, current->win.width, current->win.height, 0, 0);
	XFlush(current->win.display);
*/
      lux_clear_pixmap(current->win.display, current->win.pixmap, redrawgc, 
		       x+1, y+1, width-2, height-2);
    }

    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  
      lux_clear_pixmap(current->win.display, current->win.window, redrawgc, 
		       x+1, y+1, width-2, height-2);

}

lux_update_fg(win)
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    if (current->win.update_flag) return;

    if (current->win.pixmap) {
      XCopyArea(current->win.display, current->win.pixmap, current->win.window,
		current->win.gc, 
		0, 0, current->win.width, current->win.height, 0, 0);
      XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't clear display -- not enough memory\n"); 
	return 0;}

    newdata->type = UPDATE_FOREGROUND;

    return 1;
}

lux_reupdate_fg(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    if (current->win.pixmap) {
      XCopyArea(current->win.display, current->win.pixmap, current->win.window,
		current->win.gc, 
		0, 0, current->win.width, current->win.height, 0, 0);
      XFlush(current->win.display);
    }

    return 1;
}

lux_setup_region(win, fxorg, fyorg, fxsize, fysize)
Window win;
float  fxorg, fyorg, fxsize, fysize;
{
    register float    *ptr;
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    current->win.fxorg  = fxorg;
    current->win.fyorg  = fyorg;
    current->win.fxsize = fxsize;
    current->win.fysize = fysize;

    current->win.xsize = (int)(fxsize*current->win.user_width/10.0);
    current->win.ysize = (int)(fysize*current->win.user_height/10.0);
    current->win.xorg  = (int)(current->win.fxorg*current->win.user_width/10.0);
    current->win.yorg  = (int)((10.0-current->win.fyorg)*
			       current->win.user_height/10.0)-
			       current->win.ysize;

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't set up coord -- not enough memory\n"); 
	return 0;}

    newdata->type = SETUP_REGION;
    ptr = newdata->data.f = (float *)malloc(4*sizeof(float));
    
    if (ptr == (float *)NULL) {
      fprintf(stderr,"Can't set up coord -- not enough memory\n"); 
      fprintf(stderr,"Error may occur!!\n");
      newdata->type = NOOP;
      return 0;
    }

    ptr[0] = fxorg;
    ptr[1] = fyorg;
    ptr[2] = fxsize;
    ptr[3] = fysize;
    return 1;
}

lux_resetup_region(win)
Window win;
{
    register float *ptr;
    register lux_wins *current;
  
    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.f;

    current->win.fxorg = ptr[0];
    current->win.fyorg = ptr[1];
    current->win.fxsize = ptr[2];
    current->win.fysize = ptr[3];

    current->win.xsize = (int)(current->win.fxsize*
			       current->win.user_width/10.0+0.5);
    current->win.ysize = (int)(current->win.fysize*
			       current->win.user_height/10.0+0.5);
    current->win.xorg  = (int)(current->win.fxorg*
			       current->win.user_width/10.0+0.5);
    current->win.yorg  = (int)((10.0-current->win.fyorg)*
			       current->win.user_height/10.0+0.5)-
                               current->win.ysize;

}

lux_setup_axis(win, xmin, xmax, ymin, ymax)
Window win;
float  xmin, xmax, ymin, ymax;
{
    register float    *ptr;
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    current->win.xmin = xmin;
    current->win.xmax = xmax;
    current->win.ymin = ymin;
    current->win.ymax = ymax;

    if (current->win.lnx == 0) 
      current->win.xfactor = (float)current->win.xsize/(xmax-xmin);
    else {
	if ( xmin <= 0.0 || xmax <= 0.0) {
	    fprintf(stderr,
		    "For log axis, the limit must be greater than zero\n");
	    return 0;
	}
	current->win.xfactor = (float)current->win.xsize
				    /(log10(xmax)-log10(xmin));
    }
    if (current->win.lny == 0) 
	current->win.yfactor = (float)current->win.ysize/(ymax-ymin);
    else {
	if ( ymin <= 0.0 || ymax <= 0.0) {
	    fprintf(stderr,
		    "For log axis, the limit must be greater than zero\n");
	    return 0;
	}
	current->win.yfactor = (float)current->win.ysize
				    /(log10(ymax)-log10(ymin));
    }

    if (current->win.discard_flag) return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) {
	fprintf(stderr, "Can't set up axis -- not enough memory\n"); 
	return 0;
    }

    newdata->type = SETUP_AXIS;

    ptr = newdata->data.f = (float *)malloc(4*sizeof(float));

    if ( newdata->data.f == (float *)NULL) {
	fprintf(stderr,"Can't set up axis -- not enough memory\n"); 
	fprintf(stderr,"Error may occur!!\n"); 
	newdata->type = NOOP;
	return 0;
    }
    ptr[0] = xmin;
    ptr[1] = xmax;
    ptr[2] = ymin;
    ptr[3] = ymax;

    return 1;
}

lux_resetup_axis(win)
Window win;
{
    register float *ptr;
    register lux_wins *current;

    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.f;

    current->win.xmin = ptr[0];
    current->win.xmax = ptr[1];
    current->win.ymin = ptr[2];
    current->win.ymax = ptr[3];

    if (current->win.lnx == 0) 
	current->win.xfactor = (float)current->win.xsize/(ptr[1]-ptr[0]);
    else 
	current->win.xfactor = (float)current->win.xsize
				    /(log10(ptr[1])-log10(ptr[0]));
    if (current->win.lny == 0) 
	current->win.yfactor = (float)current->win.ysize/(ptr[3]-ptr[2]);
    else 
	current->win.yfactor = (float)current->win.ysize
				    /(log10(ptr[3])-log10(ptr[2]));
}

lux_setup_axis_style(win, lnx, lny)
Window win;
int lnx, lny;
{
    register lux_wins *current;
    register lux_data *newdata;
    int               *ptr;

    current = get_currentwin(win);

    current->win.lnx = lnx;
    current->win.lny = lny;

    if (current->win.lnx == 0) 
	current->win.xfactor=(float)current->win.xsize
				  /(current->win.xmax-current->win.xmin);
    else {
	if ( current->win.xmin <= 0.0 || current->win.xmax <= 0.0) {
	    fprintf(stderr,
		    "For log axis, the limit must be greater than zero\n");
	    return 0;
	}
	current->win.xfactor=(float)current->win.xsize
		    /(log10(current->win.xmax)-log10(current->win.xmin));
    } 
    if (current->win.lny == 0) 
	current->win.yfactor=(float)current->win.ysize
		    /(current->win.ymax-current->win.ymin);      
    else {
	if ( current->win.ymin <= 0.0 || current->win.ymax <= 0.0) {
	    fprintf(stderr,
		    "For log axis, the limit must be greater than zero\n");
	    return 0;
	}
	current->win.yfactor=(float)current->win.ysize
		    /(log10(current->win.ymax)-log10(current->win.ymin));
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) {
	fprintf(stderr,"Can't set up axis style -- not enough memory\n"); 
	return 0;
    }

    newdata->type = SETUP_AXIS_STYLE;

    ptr = newdata->data.i = (int *)malloc(2*sizeof(int));

    if ( newdata->data.i == (int *)NULL) {
	fprintf(stderr,"Can't set up axis style -- not enough memory\n"); 
	fprintf(stderr,"Error may occur!!\n"); 
	newdata->type = NOOP;
	return 0;
    }
    ptr[0] = lnx;
    ptr[1] = lny;

    return 1;     
}

lux_resetup_axis_style(win)
Window win;
{
    register int *ptr;
    register lux_wins *current;

    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.i;

    current->win.lnx = ptr[0];
    current->win.lny = ptr[1];

    if (current->win.lnx == 0)
	current->win.xfactor=(float)current->win.xsize
			  /(current->win.xmax-current->win.xmin);
    else
	current->win.xfactor=(float)current->win.xsize
		    /(log10(current->win.xmax)-log10(current->win.xmin));
    if (current->win.lny == 0)
	current->win.yfactor=(float)current->win.ysize
		    /(current->win.ymax-current->win.ymin);      
    else
	current->win.yfactor=(float)current->win.ysize
		    /(log10(current->win.ymax)-log10(current->win.ymin));
}

lux_draw_line(win, x1, y1, x2, y2)
Window win;
int x1, y1, x2, y2;
{
    register lux_wins *current;
    register lux_data *newdata;
    register int      *ptr;
    int xx1,yy1,xx2,yy2;

    current = get_currentwin(win);

    xx1 = (int)(x1*current->win.xresizefactor+0.5);
    yy1 = (int)(y1*current->win.yresizefactor+0.5);
    xx2 = (int)(x2*current->win.xresizefactor+0.5);
    yy2 = (int)(y2*current->win.yresizefactor+0.5);

    if (current->win.pixmap) 
	XDrawLine(current->win.display,current->win.pixmap,
		  current->win.gc,xx1,yy1,xx2,yy2);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawLine(current->win.display,current->win.window,
		  current->win.gc,xx1,yy1,xx2,yy2);
	XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) {
	fprintf(stderr,"Can't draw line -- not enough memory\n"); 
	return 0;
    }

    newdata->type = DRAW_LINE;
    ptr = newdata->data.i = (int *)malloc(4*sizeof(int));

    if (ptr == (int *)NULL) {
	fprintf(stderr,"Can't draw line -- not enough memory\n"); 
	fprintf(stderr,"Error may occur!!\n"); 
	newdata->type = NOOP;
	return 0;
    }

    ptr[0] = x1;
    ptr[1] = y1;
    ptr[2] = x2;
    ptr[3] = y2;

    return 1;
}

lux_redraw_line(win)
Window win;
{
    register int *ptr;
    register lux_wins *current;
    int      x1, x2, y1, y2;
  
    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.i;
    x1  = ptr[0]*current->win.xresizefactor+0.5;
    y1  = ptr[1]*current->win.yresizefactor+0.5;
    x2  = ptr[2]*current->win.xresizefactor+0.5;
    y2  = ptr[3]*current->win.yresizefactor+0.5;

    if (current->win.pixmap) 
	XDrawLine(current->win.display, current->win.pixmap,
		  redrawgc, x1, y1, x2, y2);
   
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawLine(current->win.display,current->win.window,
		  redrawgc, x1, y1, x2, y2);

/*      XFlush(current->win.display); */

    }
}

lux_convert_coords(win,fx,fy,x,y,n)
Window   win;
float   *fx,*fy;
int     *x, *y;
unsigned int n;
{
    register unsigned i;
    register lux_wins *current;
    float    logxmin, logxmax, logymin, logymax;

    current = get_currentwin(win);

    if (current->win.lnx != 0) {
	logxmin = log10(current->win.xmin);
	logxmax = log10(current->win.xmax);
    }
    if (current->win.lny != 0) {
	logymin = log10(current->win.ymin);
	logymax = log10(current->win.ymax);
    }

    if (fx == fy)  /* x == y too */
      for(i=0;i<2*n;i+=2) {
	  if ( current->win.lnx == 0)
	      x[i] = (int)(((fx[i] - current->win.xmin)*current->win.xfactor +
			  current->win.xorg)*current->win.xresizefactor + 0.5);
	  else 
	      x[i] = (int)(((log10(fx[i]) - logxmin)*current->win.xfactor +
			  current->win.xorg)*current->win.xresizefactor + 0.5);
	  if ( current->win.lny == 0)
	      y[i+1] = (int)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			    current->win.yfactor + current->win.yorg)*
			   current->win.yresizefactor + 0.5);
	  else
	      y[i+1] = (int)((current->win.ysize - (log10(fy[i+1])-logymin)*
			    current->win.yfactor + current->win.yorg)*
			    current->win.yresizefactor + 0.5);
      }  
    else if (x == y) 
      for(i=0;i<n;i++) {
	  if ( current->win.lnx == 0)
	    x[i*2] = (int)(((fx[i] - current->win.xmin)*current->win.xfactor +
			    current->win.xorg)*current->win.xresizefactor
			   						+ 0.5);
	  else 
	    x[i*2] = (int)(((log10(fx[i]) - logxmin)*current->win.xfactor +
			    current->win.xorg)*current->win.xresizefactor
			   						+ 0.5);
	  if ( current->win.lny == 0)
	    x[i*2+1] = (int)((current->win.ysize - (fy[i]-current->win.ymin)*
			      current->win.yfactor + current->win.yorg)*
			     current->win.yresizefactor + 0.5);
	  else
	    x[i*2+1] = (int)((current->win.ysize - (log10(fy[i])-logymin)*
			      current->win.yfactor + current->win.yorg)*
			     current->win.yresizefactor + 0.5);
	
      }        
    else 
      for(i=0;i<n;i++) {
	  if ( current->win.lnx == 0)
	    x[i] = (int)(((fx[i]-current->win.xmin)*current->win.xfactor+
			  current->win.xorg)*current->win.xresizefactor+0.5);
	  else 
	    x[i] = (int)(((log10(fx[i])-logxmin)*current->win.xfactor+
			  current->win.xorg)*current->win.xresizefactor+0.5);
	  if ( current->win.lny == 0)
	    y[i] = (int)((current->win.ysize-(fy[i]-current->win.ymin)*
			  current->win.yfactor+current->win.yorg)*
			 current->win.yresizefactor+0.5);
	  else
	    y[i] = (int)((current->win.ysize-(log10(fy[i])-logymin)*
			  current->win.yfactor+current->win.yorg)*
			 current->win.yresizefactor+0.5);
      }
    
}

lux_convert_coord(win,fx,fy,x,y)
Window   win;
float    fx,fy;
int     *x, *y;
{
    register lux_wins *current;
    float    logxmin, logxmax, logymin, logymax;

    current = get_currentwin(win);

    if (current->win.lnx != 0) {
      logxmin = log10(current->win.xmin);
      logxmax = log10(current->win.xmax);
    }
    if (current->win.lny != 0) {
      logymin = log10(current->win.ymin);
      logymax = log10(current->win.ymax);
    }

    if ( current->win.lnx == 0)
      *x = (int)(((fx-current->win.xmin)*current->win.xfactor+
		  current->win.xorg)*current->win.xresizefactor+0.5);
    else
      *x = (int)(((log10(fx)-logxmin)*current->win.xfactor+
		  current->win.xorg)*current->win.xresizefactor+0.5);
    if ( current->win.lny == 0)
      *y = (int)((current->win.ysize-(fy-current->win.ymin)*
		  current->win.yfactor+current->win.yorg)*
		 current->win.yresizefactor+0.5);
    else
      *y = (int)((current->win.ysize-(log10(fy)-logymin)*
		  current->win.yfactor+current->win.yorg)*
		 current->win.yresizefactor+0.5);
      
}      

lux_reconvert_coord(win,ix,iy,fx,fy)
Window   win;
int      ix, iy;
float    *fx,*fy;
{
    register long i;
    register lux_wins *current;

    current = get_currentwin(win);
/*
    *x = (int)(((fx-current->win.xmin)*current->win.xfactor+
		 current->win.xorg)*current->win.xresizefactor+0.5);
    *y = (int)((current->win.ysize-(fy-current->win.ymin)*
		current->win.yfactor+current->win.yorg)*
		current->win.yresizefactor+0.5);
*/
    if ( current->win.lnx == 0)
      *fx = (((float)ix)/current->win.xresizefactor-current->win.xorg)/
	      current->win.xfactor+current->win.xmin;
    else
      *fx = pow(10.0,((((float)ix)/current->win.xresizefactor
		       -current->win.xorg)/
	      current->win.xfactor+log10(current->win.xmin)));
    if ( current->win.lny == 0)
      *fy = (current->win.ysize-((float)iy)/current->win.yresizefactor+
	     current->win.yorg)/current->win.yfactor+current->win.ymin;
    else
      *fy = pow(10.0,((current->win.ysize-((float)iy)
		       /current->win.yresizefactor+
		       current->win.yorg)/current->win.yfactor
		      +log10(current->win.ymin)));
        
}
      
lux_reconvert_rcoord(win,ix,iy,fx,fy)/*input coord are relative to the old one*/
Window   win;
int      ix, iy;
float    *fx,*fy;
{
    register long i;
    register lux_wins *current;

    current = get_currentwin(win);

    if ( current->win.lnx == 0)
      *fx = (((float)ix)-current->win.xorg)/
		  current->win.xfactor+current->win.xmin;
    else
      *fx = pow(10.0,((((float)ix)-current->win.xorg)/
		  current->win.xfactor+log10(current->win.xmin)));
    if ( current->win.lny == 0)
      *fy = (current->win.ysize-((float)iy)+current->win.yorg)/
	    current->win.yfactor+current->win.ymin;
    else
      *fy = pow(10.0,((current->win.ysize-((float)iy)+current->win.yorg)/
	    current->win.yfactor+log10(current->win.ymin)));

}      


lux_draw_linef(win, x1, y1, x2, y2)
Window win;
float x1, y1, x2, y2;
{
    register lux_wins *current;
    register lux_data *newdata;
    register float    *ptr;
    int      xx1,yy1,xx2,yy2;

    current = get_currentwin(win);

    lux_convert_coord(win,x1,y1,&xx1,&yy1);
    lux_convert_coord(win,x2,y2,&xx2,&yy2);

    if (current->win.pixmap) 
      XDrawLine(current->win.display,current->win.pixmap,
		current->win.gc,xx1,yy1,xx2,yy2);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawLine(current->win.display,current->win.window,
		current->win.gc,xx1,yy1,xx2,yy2);
      XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't draw line -- not enough memory\n"); 
	return 0;}

    newdata->type = DRAW_LINE_F;
    ptr = newdata->data.f = (float *)malloc(4*sizeof(float));

    if (ptr == (float *)NULL) {
      fprintf(stderr,"Can't draw line -- not enough memory\n"); 
      fprintf(stderr,"Error may occur!!\n"); 
      newdata->type = NOOP;
      return 0;
    }


    ptr[0] = x1;
    ptr[1] = y1;
    ptr[2] = x2;
    ptr[3] = y2;

    return 1;
}

lux_redraw_linef(win)
Window win;
{
    register float *ptr;
    register lux_wins *current;
    int      x1, x2, y1, y2;
    int      data[4];
  
    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.f;

    lux_convert_coords(win,ptr,ptr,data,data,2);

    if (current->win.pixmap) 
      XDrawLine(current->win.display, current->win.pixmap,
		redrawgc, data[0], data[1], data[2], data[3]);
   
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawLine(current->win.display,current->win.window,
		redrawgc, data[0], data[1], data[2], data[3]);
 /*     XFlush(current->win.display);*/
    }

}

lux_draw_pointf(win, x1, y1)
Window win;
float x1, y1;
{
    register lux_wins *current;
    register lux_data *newdata;
    register float    *ptr;
    int      xx1,yy1;

    current = get_currentwin(win);

    lux_convert_coord(win,x1,y1,&xx1,&yy1);

    if (current->win.pixmap) 
      XDrawPoint(current->win.display,current->win.pixmap,
		current->win.gc,xx1,yy1);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawPoint(current->win.display,current->win.window,
		 current->win.gc,xx1,yy1);
      XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't draw line -- not enough memory\n"); 
	return 0;}

    newdata->type = DRAW_POINT_F;
    ptr = newdata->data.f = (float *)malloc(4*sizeof(float));

    if (ptr == (float *)NULL) {
      fprintf(stderr,"Can't draw line -- not enough memory\n"); 
      fprintf(stderr,"Error may occur!!\n"); 
      newdata->type = NOOP;
      return 0;
    }

    ptr[0] = x1;
    ptr[1] = y1;
    return 1;
}

lux_redraw_pointf(win)
Window win;
{
    register float *ptr;
    register lux_wins *current;
    int      x1, x2, y1, y2;
    int      data[2];
  
    current = get_currentwin(win);

    ptr = (current->win.currentdata)->data.f;

    lux_convert_coords(win,ptr,ptr,data,data,1);

    if (current->win.pixmap) 
      XDrawPoint(current->win.display, current->win.pixmap,
		redrawgc, data[0], data[1]);
   
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawPoint(current->win.display,current->win.window,
		redrawgc, data[0], data[1]);
 /*     XFlush(current->win.display);*/
    }

}

lux_set_linewidth(win,width)
Window win;
unsigned int width;
{
    register lux_wins  *current;
    register lux_data  *newdata;
    XGCValues  values;

    current = get_currentwin(win);

    XGetGCValues(current->win.display, current->win.gc, 
		 GCLineStyle | GCCapStyle | GCJoinStyle,
		 &values);

    XSetLineAttributes(current->win.display, current->win.gc, width, 
		      values.line_style,values.cap_style,values.join_style);

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't set line width -- not enough memory\n"); 
	return 0;}

    newdata->type = SET_LINE_WIDTH;
    newdata->data.u = (unsigned int *)malloc(sizeof(int));
    if ( newdata->data.u == (unsigned int *)NULL) {
      fprintf(stderr,"Can't setlinewidth -- not enough memory\n"); 
      fprintf(stderr,"Error may occur!!\n"); 
      newdata->type = NOOP;
      return 0;
    }

    newdata->data.u[0] = width;
    return 1;
}

lux_reset_linewidth(win)
Window win;
{
    register lux_wins  *current;
    unsigned int width;
    XGCValues  values;

    current = get_currentwin(win);

    width = (current->win.currentdata)->data.u[0];

    XGetGCValues(current->win.display, redrawgc, 
		 GCLineStyle | GCCapStyle | GCJoinStyle,
		 &values);

    XSetLineAttributes(current->win.display, redrawgc, width, 
		      values.line_style,values.cap_style,values.join_style);
}

#define MAX_LINE_STYLE           11 

#define SHORT_DASHED_LIST_LENGTH 2
#define LONG_DASHED_LIST_LENGTH  2
#define DOTTED_LIST_LENGTH       2
#define DOT_DASHED_LIST_LENGTH   4
#define ODD_DASHED_LIST_LENGTH   3

static int dash_list_length[] = {
    SHORT_DASHED_LIST_LENGTH,
    LONG_DASHED_LIST_LENGTH,
    DOTTED_LIST_LENGTH,
    DOT_DASHED_LIST_LENGTH,
    ODD_DASHED_LIST_LENGTH
};


/* must be at least one element in each list */

static char short_dashed[SHORT_DASHED_LIST_LENGTH] = {4, 4};
static char long_dashed[LONG_DASHED_LIST_LENGTH] = {4, 7};
static char dotted[DOTTED_LIST_LENGTH] = {3, 1};
static char dot_dashed[DOT_DASHED_LIST_LENGTH] = {3, 4, 3, 1};
static char odd_dashed[ODD_DASHED_LIST_LENGTH] = {1,2,3};

static char *dash_list[] = {
    short_dashed,
    long_dashed,
    dotted,
    dot_dashed,
    odd_dashed,
};

static int dash_offset = 0;

lux_set_linestyle(win,style)
Window win;
int    style; /* 0 Solid line, 1-5 OnOffDash, 6-10 DoubleDash(with background */
{
    register lux_wins  *current;
    register lux_data  *newdata;
    XGCValues  values;
    int      line_style;

    current = get_currentwin(win);

    XGetGCValues(current->win.display, current->win.gc, 
		 GCLineWidth | GCCapStyle | GCJoinStyle,
		 &values);

    if (style>=MAX_LINE_STYLE) style = (style>0)?style%(MAX_LINE_STYLE):(-style)%(MAX_LINE_STYLE);

    line_style = style;

    if (style>5) {
	XSetDashes(current->win.display, current->win.gc, dash_offset, 
		   dash_list[style-6], dash_list_length[style-6]);
	style = 2; /* Set to LineDoubleDash (With background color)*/
    } else if (style>0) {
	XSetDashes(current->win.display, current->win.gc, dash_offset, 
		   dash_list[style-1], dash_list_length[style-1]);
	style = 1; /* Set to LineOnOffDash */
    }	
      

    XSetLineAttributes(current->win.display,current->win.gc, values.line_width, 
		       style,values.cap_style,values.join_style);

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr,"Can't set line style -- not enough memory\n"); 
	return 0;}

    newdata->type = SET_LINE_STYLE;
    newdata->data.u = (unsigned int *)malloc(sizeof(int));
    if ( newdata->data.u == (unsigned int *)NULL) {
      fprintf(stderr,"Can't setlinestyle -- not enough memory\n"); 
      fprintf(stderr,"Error may occur!!\n"); 
      newdata->type = NOOP;
      return 0;
    }

    newdata->data.i[0] = line_style;  /* Save the old value */
    return 1;
}

lux_reset_linestyle(win)
Window win;
{
    register lux_wins  *current;
    int      style;
    XGCValues  values;

    current = get_currentwin(win);

    style = (current->win.currentdata)->data.i[0];

    XGetGCValues(current->win.display, redrawgc, 
		 GCLineWidth | GCCapStyle | GCJoinStyle,
		 &values);

    if (style>5) {
	XSetDashes(current->win.display, redrawgc, dash_offset, 
		   dash_list[style-6], dash_list_length[style-6]);
	style = 2; /* Set to LineDoubleDash (With background color)*/
    } else if (style>0) {
	XSetDashes(current->win.display, redrawgc, dash_offset, 
		   dash_list[style-1], dash_list_length[style-1]);
	style = 1; /* Set to LineOnOffDash */
    }	
      
    XSetLineAttributes(current->win.display, redrawgc, values.line_width, 
		       style,values.cap_style,values.join_style);
}

lux_string_dimensions(char* s, int* length, int* height) /* Results in PIXELS */
{
    *length = XTextWidth(font_info, s, strlen(s));
    *height = font_info->ascent + font_info->descent;
}



static int buffer_gfx = 1;

void lux_set_nobuffer()
{
    /* Suppress buffering for *next* lux_draw*string call.
       Flag buffer_gfx is ALWAYS reset to 1 in lux_draw*string. */

    buffer_gfx = 0;
}



lux_draw_string(win, fx, fy, fh, s, position)
Window  win;
float   fx,fy,fh;     /* x,y position, h is offset in units of font height */ 
char   *s, position;  /* -1 -> left, 0 -> center, 1->right */
{
    int ix,iy;
    lux_wins *current;
    lux_data *newdata;
    char    *string;
    int      stringwidth, slength;

    current = get_currentwin(win);
    
    lux_convert_coord(win, fx, fy, &ix, &iy);

    slength = strlen(s);    
    stringwidth = XTextWidth(font_info, s, slength);

    if (position == 0) ix = ix - (stringwidth+1)/2;
    else if (position == 1) ix = ix-stringwidth;

    iy -= (int)(0.5 + fh * (font_info->ascent+font_info->descent)); 

    if (current->win.pixmap) 
      XDrawString(current->win.display, current->win.pixmap, current->win.gc,
		  ix, iy, s, slength);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawString(current->win.display, current->win.window, current->win.gc,
		    ix, iy, s, slength);
	XFlush(current->win.display);
    }

    if (current->win.discard_flag)  return 1;

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) {
	fprintf(stderr,"Can't set line width -- not enough memory\n"); 
	return 0;
    }

    if (buffer_gfx) {

	newdata->type = DRAW_STRING;
	newdata->data.f=
	    (float *)malloc((slength+2)*sizeof(char)+3*sizeof(float));

	if ( newdata->data.f == (float *)NULL) {
	    fprintf(stderr,"Can't drawstring -- not enough memory\n"); 
	    fprintf(stderr,"Error may occur!!\n"); 
	    newdata->type = NOOP;
	    return 0;
	}

	string = &newdata->data.b[3*sizeof(float)];

	newdata->data.f[0] = fx;
	newdata->data.f[1] = fy;
	newdata->data.f[2] = fh;
	string[0] = position;
	strcpy(&string[1],s);
    }

    buffer_gfx = 1;
    return 1;
} 

lux_redraw_string(win)
Window  win;
{
    int ix,iy;
    char *string;
    register lux_wins *current;
    register lux_data *newdata;
    int      stringwidth, slength;
    char     position;

    current = get_currentwin(win);

    lux_convert_coord(win, (current->win.currentdata)->data.f[0],
		      (current->win.currentdata)->data.f[1],&ix, &iy);

    string = &((current->win.currentdata)->data.b[3*sizeof(float)]);

    slength = strlen(&string[1]);    
    stringwidth = XTextWidth(font_info, &string[1], slength);

    position = string[0];
    if (position == 0) ix = ix - (stringwidth+1)/2;
    else if (position == 1) ix = ix-stringwidth;

    iy = iy - (int)(0.5+(current->win.currentdata)->data.f[2]
                        *(font_info->ascent+font_info->descent));

    if (current->win.pixmap) 
      XDrawString(current->win.display, current->win.pixmap, redrawgc,
		  ix, iy, &string[1], slength);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XDrawString(current->win.display, current->win.window, redrawgc,
		  ix, iy, &string[1], slength);

/*      XFlush(current->win.display); */
    
}

lux_draw_image_string(win, fx, fy, fh, s, position)
Window  win;
float   fx,fy,fh;     /* x,y position, h is times of font height */ 
char   *s, position;  /* -1 -> left, 0 -> center, 1->right */
{
    int ix,iy;
    lux_wins *current;
    lux_data *newdata;
    char     *string;
    int      stringwidth, slength;

    current = get_currentwin(win);

    lux_convert_coord(win, fx, fy, &ix, &iy);

    slength = strlen(s);    
    stringwidth = XTextWidth(font_info, s, slength);

    if (position == 0) ix = ix - (stringwidth+1)/2;
    else if (position == 1) ix = ix-stringwidth;

    iy = iy - (int)(0.5+fh*(font_info->ascent+font_info->descent));

    if (current->win.pixmap) 
      XDrawImageString(current->win.display, current->win.pixmap, 
		       current->win.gc, ix, iy, s, slength);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawImageString(current->win.display, current->win.window, 
		       current->win.gc, ix, iy, s, slength);
      XFlush(current->win.display);
    }

    if (current->win.discard_flag)
	return 1;

    if (buffer_gfx) {

	/* Buffer graphics operations for replay. */

	newdata = get_newdata(current);

	if ( newdata == (lux_data *)NULL ) {
	    fprintf(stderr,"Can't set line width -- not enough memory\n"); 
	    return 0;
	}

	newdata->type = DRAW_IMAGE_STRING;
	newdata->data.f=(float *)malloc((slength+2)*sizeof(char)
					 + 3*sizeof(float));

	if ( newdata->data.f == (float *)NULL) {
	    fprintf(stderr,"Can't draw image string -- not enough memory\n"); 
	    fprintf(stderr,"Error may occur!!\n"); 
	    newdata->type = NOOP;
	    return 0;
	}
	string = &newdata->data.b[3*sizeof(float)];

	newdata->data.f[0] = fx;
	newdata->data.f[1] = fy;
	newdata->data.f[2] = fh;
	string[0] = position;
	strcpy(&string[1], s);
    }

    buffer_gfx = 1;
    return 1;
} 

lux_redraw_image_string(win)
Window  win;
{
    int ix,iy;
    char *string;
    register lux_wins *current;
    register lux_data *newdata;
    int      stringwidth, slength;
    char     position;

    current = get_currentwin(win);

    lux_convert_coord(win, (current->win.currentdata)->data.f[0],
		      (current->win.currentdata)->data.f[1],&ix, &iy);

    string = &((current->win.currentdata)->data.b[3*sizeof(float)]);

    slength = strlen(&string[1]);    
    stringwidth = XTextWidth(font_info, &string[1], slength);

    position = string[0];
    if (position == 0) ix = ix - (stringwidth+1)/2;
    else if (position == 1) ix = ix-stringwidth;

    iy = iy - (int)(0.5+(current->win.currentdata)->data.f[2]
                        *(font_info->ascent+font_info->descent));

    if (current->win.pixmap) 
      XDrawImageString(current->win.display, current->win.pixmap, redrawgc,
		  ix, iy, &string[1], slength);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XDrawImageString(current->win.display, current->win.window, redrawgc,
		  ix, iy, &string[1], slength);
/*      XFlush(current->win.display);*/
    
} 

lux_draw_vstring(win, fx, fy, fh, s, position) 
Window  win;
float   fx,fy,fh;
char   *s, position;  /* -1 -> left, 0 -> center, 1->right */
{
    int ix,iy;
    register lux_wins *current;
    register lux_data *newdata;
    int      stringwidth, slength, charheight, height, *ptr;
    Pixmap   pixmap, pixmapwhite;
    XImage   *in, *out;
    register int i,j;
    unsigned long c;
    XGCValues value;

    current = get_currentwin(win);
    
    lux_convert_coord(current->win.window, fx, fy, &ix, &iy);

    slength = strlen(s);    
    stringwidth = XTextWidth(font_info, s, slength);
    charheight = font_info->ascent+font_info->descent;
    pixmap = XCreatePixmap(current->win.display, current->win.window, 
		  stringwidth, charheight,  current->win.window_depth);
    if (!pmpallocflag) {
      fprintf(stderr, "In draw_vstring, not enough video memory\n");
      pixmap=(Pixmap)NULL; pmpallocflag=1;return 0;}

    XDrawImageString(current->win.display, pixmap, current->win.gc,
		0, font_info->ascent, s, slength);
/*    XCopyArea(current->win.display, pixmap, current->win.pixmap, 
		current->win.gc, 0, 0, stringwidth, charheight, 0, 0);
*/
    in = XGetImage(current->win.display, pixmap, 
		      0, 0, stringwidth, charheight, AllPlanes, XYPixmap);

/*    XPutImage(current->win.display, current->win.pixmap, current->win.gc,
	      in, 0, 0, 0, 0, stringwidth, charheight);*/ /* Looks Bad ! */

    XFreePixmap(current->win.display, pixmap);

    pixmap = XCreatePixmap(current->win.display, current->win.window, 
			   charheight, stringwidth, current->win.window_depth);
    if (!pmpallocflag) {
      fprintf(stderr, "In draw_vstring, not enough video memory\n");
      pixmap=(Pixmap)NULL; pmpallocflag=1;
      XDestroyImage(in);
      return 0;}
    pixmapwhite = XCreatePixmap(current->win.display, current->win.window, 
			   charheight, stringwidth, current->win.window_depth);
    if (!pmpallocflag) {
      fprintf(stderr, "In draw_vstring, not enough video memory\n");
      pixmap=(Pixmap)NULL; pmpallocflag=1;
      XFreePixmap(current->win.display,pixmap);
      XDestroyImage(in);
      return 0;}

    XGetGCValues(current->win.display, current->win.gc,
		 GCBackground | GCForeground, &value);
    c = value.foreground;
    XSetForeground(current->win.display, current->win.gc, value.background);
    lux_clear_pixmap(current->win.display, pixmapwhite, current->win.gc,
		     0, 0, charheight, stringwidth);
    XSetForeground(current->win.display, current->win.gc, c);

    out = XGetImage(current->win.display, pixmap, 
		      0, 0, charheight, stringwidth, AllPlanes, XYPixmap);

    for(i=0;i<stringwidth;i++)
      for(j=0;j<charheight;j++) 
	XPutPixel(out, j, stringwidth-i-1, XGetPixel(in, i, j));

    XPutImage(current->win.display, pixmap, current->win.gc,
	      out, 0, 0, 0, 0, charheight, stringwidth); /* Looks Bad ! */

    if (position == (char)(-1)) {ix=ix-charheight;iy=iy-stringwidth;}
    else if (position == 0) {ix=ix-charheight;iy=iy-(stringwidth+1)/2;}
    else if (position == 1) {ix=ix-charheight;iy=iy-stringwidth;}

    ix = ix - (int)(0.5 + fh*charheight);

    XGetGCValues(current->win.display, current->win.gc, GCFunction, &value);
    XSetFunction(current->win.display, current->win.gc, GXxor);
    XCopyArea(current->win.display, pixmapwhite, pixmap, 
		current->win.gc, 0, 0, charheight, stringwidth, 0, 0);

    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmapwhite, current->win.pixmap, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XCopyArea(current->win.display, pixmapwhite, current->win.window, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);
    XSetFunction(current->win.display, current->win.gc, GXor);

    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmap, current->win.pixmap, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))
      XCopyArea(current->win.display, pixmap, current->win.window, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);

    XSetFunction(current->win.display, current->win.gc, GXxor);

    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmapwhite, current->win.pixmap, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XCopyArea(current->win.display, pixmapwhite, current->win.window, 
		current->win.gc, 0, 0, charheight, stringwidth, ix, iy);
      XFlush(current->win.display);
    }

    XSetFunction(current->win.display, current->win.gc, value.function);

    XFreePixmap(current->win.display, pixmap);
    XFreePixmap(current->win.display, pixmapwhite);
    XDestroyImage(in);

    if (current->win.discard_flag)  {XDestroyImage(out);return 1;}

    if (buffer_gfx) {

	if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) {
	    fprintf(stderr,"Can't draw_vstring -- not enough memory\n"); 
	    return 0;
	}

	newdata->type = DRAW_VSTRING;
	newdata->data.b = (char *)malloc(3*sizeof(float) + 2*sizeof(int)
					 + sizeof(XImage*)+1);

	if (newdata->data.b == (char *)NULL) {
	    fprintf(stderr,"Can't draw_vstring -- not enough memory\n"); 
	    fprintf(stderr,"Error may occur!!\n"); 
	    newdata->type = NOOP;
	    return 0;
	}

	newdata->data.f[0] = fx;
	newdata->data.f[1] = fy;
	newdata->data.f[2] = fh;
	ptr = (int *)&newdata->data.f[3];
	ptr[0] = charheight;
	ptr[1] = stringwidth;


	/* NO idea what Biao was doing here...  SLWM 4/98.
	   However, it is necessary when redrawing windows... */

	*((XImage **)(&ptr[2])) = out;	/* Causes "unaligned access" warning
					   in DEC UNIX... */


	newdata->data.b[3*sizeof(float)+2*sizeof(int)+sizeof(XImage*)]
								    = position;
    }

    buffer_gfx = 1;
    return 1;
} 

lux_redraw_vstring(win)
Window  win;
{
    int       ix,iy,*ptr;
    register  lux_wins *current;
    int       stringwidth, charheight;
    unsigned  long c;
    Pixmap    pixmap, pixmapwhite;
    XImage   *out;
    XGCValues value;
    float     fx,fy,fh;
    char      position;


    current = get_currentwin(win);

    fx = (current->win.currentdata)->data.f[0];
    fy = (current->win.currentdata)->data.f[1];
    fh = (current->win.currentdata)->data.f[2];
    ptr = (int*)(&(current->win.currentdata)->data.f[3]);
    charheight = ptr[0];
    stringwidth = ptr[1];
    out = ((XImage **)(&ptr[2]))[0];
    position = (current->win.currentdata)->data.b[3*sizeof(float)+
				     2*sizeof(int)+sizeof(XImage*)];

    lux_convert_coord(current->win.window, fx, fy, &ix, &iy);

    pixmap = XCreatePixmap(current->win.display, current->win.window, 
			   charheight, stringwidth, current->win.window_depth);
    if (!pmpallocflag) {pixmap=(Pixmap)NULL; pmpallocflag=1;
      fprintf(stderr, "In redraw_vstring, not enough video memory\n");
      return 0;
    }
    pixmapwhite = XCreatePixmap(current->win.display, current->win.window, 
			   charheight, stringwidth, current->win.window_depth);
    if (!pmpallocflag) {pixmap=(Pixmap)NULL; pmpallocflag=1;
      fprintf(stderr, "In redraw_vstring, not enough video memory\n");
      XFreePixmap(current->win.display, pixmap);
      return 0;
    }

    XGetGCValues(current->win.display, redrawgc,
		 GCBackground | GCForeground, &value);
    c = value.foreground;
    XSetForeground(current->win.display, redrawgc, value.background);
    lux_clear_pixmap(current->win.display, pixmapwhite, redrawgc,
		     0, 0, charheight, stringwidth);
    XSetForeground(current->win.display, redrawgc, c);

    XPutImage(current->win.display, pixmap, redrawgc,
	      out, 0, 0, 0, 0, charheight, stringwidth);

    if (position == (char)(-1)) {ix=ix-charheight;iy=iy-stringwidth;}
    else if (position == 0) {ix=ix-charheight;iy=iy-(stringwidth+1)/2;}
    else if (position == 1) {ix=ix-charheight;iy=iy-stringwidth;}

    ix = ix - (int)(0.5 + fh*charheight);

    XGetGCValues(current->win.display, redrawgc, GCFunction, &value);
    XSetFunction(current->win.display, redrawgc, GXxor);
    XCopyArea(current->win.display, pixmapwhite, pixmap, 
		redrawgc, 0, 0, charheight, stringwidth, 0, 0);
    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmapwhite, current->win.pixmap, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XCopyArea(current->win.display, pixmapwhite, current->win.window, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);
    XSetFunction(current->win.display, redrawgc, GXor);
    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmap, current->win.pixmap, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XCopyArea(current->win.display, pixmap, current->win.window, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);
    XSetFunction(current->win.display, redrawgc, GXxor);
    if (current->win.pixmap) 
      XCopyArea(current->win.display, pixmapwhite, current->win.pixmap, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) 
      XCopyArea(current->win.display, pixmapwhite, current->win.window, 
		redrawgc, 0, 0, charheight, stringwidth, ix, iy);

    XSetFunction(current->win.display, redrawgc, value.function);
    XFreePixmap(current->win.display, pixmap);
    XFreePixmap(current->win.display, pixmapwhite);
    return 1;
} 

lux_draw_axis1(win)
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;
    unsigned int       x,y,width,height,i;
    char msg[20];
    float fx,fy;

    current = get_currentwin(win);

    x = current->win.xorg*current->win.xresizefactor+0.5;
    y = current->win.yorg*current->win.yresizefactor+0.5;
    width  = current->win.xsize*current->win.xresizefactor+0.5;
    height = current->win.ysize*current->win.yresizefactor+0.5;

    if (current->win.pixmap) 
      XDrawRectangle(current->win.display, current->win.pixmap, 
		     current->win.gc, x, y, width, height);

    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawRectangle(current->win.display, current->win.window, 
		     current->win.gc, x, y, width, height);
      XFlush(current->win.display);
    }

    fy = current->win.ymin;
    for(i=0;i<5;i++) {
      fx = current->win.xmin+(current->win.xmax-current->win.xmin)*i/4.0;
      sprintf(msg,"%6.2f",fx);
      lux_draw_string(current->win.window,fx,fy,-1.1,msg,0);
    }
    fx = current->win.xmin;
    for(i=0;i<5;i++) {
      fy = current->win.ymin+(current->win.ymax-current->win.ymin)*i/4.0;
      sprintf(msg,"%6.2f",fy);
      lux_draw_vstring(current->win.window,fx,fy,0.1,msg,0);
    }

    if (current->win.discard_flag) return 1;

    if ((newdata = get_newdata(current)) == (lux_data *)NULL)
      {fprintf(stderr,"Can't draw axis -- not enough memory\n"); 
       return 0;}

    newdata->type = DRAW_AXIS;
}

lux_redraw_axis1(win)
Window win;
{
    register lux_wins *current;
    unsigned int       x,y,width,height;

    current = get_currentwin(win);

    x = current->win.xorg*current->win.xresizefactor+0.5;
    y = current->win.yorg*current->win.yresizefactor+0.5;
    width  = current->win.xsize*current->win.xresizefactor+0.5;
    height = current->win.ysize*current->win.yresizefactor+0.5;

    if (current->win.pixmap) 
      XDrawRectangle(current->win.display, current->win.pixmap, redrawgc,
		     x, y, width, height);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawRectangle(current->win.display, current->win.window, redrawgc,
		     x, y, width, height);
/*      XFlush(current->win.display);*/
    }
}

redraw(win, status)
Window win;
int status;
{

    register lux_wins *current;

    current = get_currentwin(win);

    if (current->win.data == (lux_data *)NULL) return;

/* Important for redraw if you have large amount of points */
    XSync(current->win.display,False);

    if (lux_peek_expose(current->win.display,current->win.window)) return;

    if ( status == 0 ) { /* added for movie */
	if (current->win.old_width  == current->win.width  &&
	    current->win.old_height == current->win.height &&
	    current->win.pixmap) {
	    XCopyArea(current->win.display, current->win.pixmap, current->win.window,
		      current->win.gc, 0, 0, current->win.width, current->win.height,
		      0, 0);
	    return;
	}            

	if (current->win.pixmap) {
	    XFreePixmap(current->win.display, current->win.pixmap);

/*If you remove comment it will insist to alloc pixmap everytime you resize*
	}
	{
 *If you remove comment it will insist to alloc pixmap everytime you resize*/

	    current->win.pixmap = XCreatePixmap(current->win.display, 
						current->win.window,
						current->win.width, current->win.height,
						current->win.window_depth);

	    XSync(current->win.display,False); 
	
	    if (!pmpallocflag) {current->win.pixmap=(Pixmap)NULL; pmpallocflag=1;}

	    if (current->win.pixmap)         /* Here is to clear the pixmap  */
	      lux_clear_pixmap(current->win.display, current->win.pixmap, defaultgc,
			       0, 0, current->win.width, current->win.height);
	}

	if ((current->win.old_width  != current->win.width  ||
	     current->win.old_height != current->win.height) /*&&
	    current->win.pixmap == (Pixmap)NULL*/)
	  XClearWindow(current->win.display,current->win.window);
    }
    else {
	if (current->win.pixmap)         /* Here is to clear the pixmap  */
	  lux_clear_pixmap(current->win.display, current->win.pixmap, defaultgc,
			   0, 0, current->win.width, current->win.height);
        XClearWindow(current->win.display,current->win.window);
    }

    XCopyGC(current->win.display, defaultgc, (long)0x003FFFFF, redrawgc);

    current->win.currentdata = current->win.data;

    while(1) {

      switch((current->win.currentdata)->type) {
      case SETUP_REGION:
	lux_resetup_region(current->win.window);
	break;
      case SETUP_AXIS:
	lux_resetup_axis(current->win.window);
	break;
      case SETUP_AXIS_STYLE:
	lux_resetup_axis_style(current->win.window);
	break;
      case DRAW_AXIS:
	lux_redraw_axis1(current->win.window);
	break;
      case DRAW_LINE:
	lux_redraw_line(current->win.window);
	break;
      case DRAW_LINE_F:
	lux_redraw_linef(current->win.window);
	break;
      case DRAW_LINES_F:
	lux_redraw_linesf(current->win.window);
	break;
      case SET_LINE_WIDTH:
	lux_reset_linewidth(current->win.window);
	break;
      case SET_LINE_STYLE:
	lux_reset_linestyle(current->win.window);
	break;
      case SET_COLOR:
	lux_reset_color(current->win.window);
	break;
      case SET_BG_COLOR:
	lux_reset_bgcolor(current->win.window);
	break;
      case DRAW_STRING:
	lux_redraw_string(current->win.window);
	break;
      case DRAW_IMAGE_STRING:
	lux_redraw_image_string(current->win.window);
	break;
      case DRAW_VSTRING:
	lux_redraw_vstring(current->win.window);
	break;
      case DRAW_SEGMENTS_F:
	lux_redraw_segmentsf(current->win.window);
	break;
      case DRAW_RECTANGLES_F:
	lux_redraw_rectanglesf(current->win.window);
	break;
      case DRAW_RECTANGLE_F:
	lux_redraw_rectanglef(current->win.window);
	break;
      case DRAW_ARCS_F:
	lux_redraw_arcsf(current->win.window);
	break;
      case DRAW_ARC_F:
	lux_redraw_arcf(current->win.window);
	break;
      case FILL_RECTANGLES_F:
	lux_refill_rectanglesf(current->win.window);
	break;
      case FILL_RECTANGLE_F:
	lux_refill_rectanglef(current->win.window);
	break;
      case FILL_ARCS_F:
	lux_refill_arcsf(current->win.window);
	break;
      case FILL_ARC_F:
	lux_refill_arcf(current->win.window);
	break;
      case FILL_POLYGON_F:
	lux_refill_polygonf(current->win.window);
	break;
      case DRAW_POINTS_F:
	lux_redraw_pointsf(current->win.window);
	break;
      case DRAW_POINT_F:
	lux_redraw_pointf(current->win.window);
	break;
      case SET_UPDATE:
	lux_reset_update(current->win.window);
	break;
      case SET_NO_UPDATE:
	lux_reset_noupdate(current->win.window);
	break;
      case UPDATE_FOREGROUND:
	lux_reupdate_fg(current->win.window);
/* Added for pause in animation */
	if (lux_check_keypress(current->win.window, 'p')) {
	    XBell(current->win.display,50);
	    while(!lux_check_keypress(current->win.window, 'c'));
	    XBell(current->win.display,50);
	}
	break;
      case CLEAR_WINDOW:
	lux_reclear_window(current->win.window);
	break;
      case CLEAR_CURRENT_REGION:   
	lux_reclear_current_region(current->win.window);
	break;
      default:
	break;
      }
      
      if ((current->win.currentdata)->next == current->win.currentdata) break;


      switch(lux_peek_expose_config(current->win.display,current->win.window)) {
	  case ABORT_REDRAW:
/* You must hit r to finish redraw  for these two lines*/
	  current->win.old_width  = current->win.width;
	  current->win.old_height = current->win.height;
/* You must hit r to finish redraw  for these two lines*/
	  if (current->win.pixmap) {
	      XCopyArea(current->win.display, current->win.pixmap, 
			current->win.window, redrawgc, 0, 0, 
			current->win.width, current->win.height, 0, 0);
	      XFlush(current->win.display);
	  }
	  return;
	  case ONE_MORE_EXPOSE:
	  if (current->win.pixmap == (Pixmap)NULL) {
/* This not good for catching the expose event */
	      if (current->win.old_width  == current->win.width &&
		  current->win.old_height == current->win.height) {
		  current->win.old_width  = current->win.width + 1;
		  current->win.old_height = current->win.height+ 1;
	      }
	      return;
	  }
	  break;
	  case RESIZED:
/* Check if it is really resized */
	  {
	      Window root;       
	      int x,y;
	      XEvent event; 
	      unsigned int width,height,border_width,depth;
	      XGetGeometry(current->win.display,current->win.window,
			   &root,&x,&y,&width,&height,&border_width,&depth);
	      if (current->win.width  != width ||
		  current->win.height != height) return;
/* Remove the useless events */
	      XCheckTypedWindowEvent(current->win.display,current->win.window,
				     ConfigureNotify,&event);
          }
	  break;
	  default:
	  break;
      }

      current->win.currentdata = (current->win.currentdata)->next;

    }
    current->win.old_width  = current->win.width;
    current->win.old_height = current->win.height;
    
    if (current->win.pixmap) 
      XCopyArea(current->win.display, current->win.pixmap, current->win.window,
		redrawgc, 0, 0, current->win.width, current->win.height, 0, 0);
    XFlush(current->win.display);
    
/*    fprintf(stderr,"Finish redraw in win %u\n",current->win.window);*/
}


