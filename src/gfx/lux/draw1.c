/* draw1.c
 * Biao Lu                   Drexel University
 *                           E-mail biao@eagle.drexel.edu
 */

#include "win.h"

#define MAXDRAWSIZE 8000  /* Not max size you can draw in the main program 
                             In main program the size is unlimited */

extern GC        redrawgc, defaultgc;
extern lux_wins *get_currentwin();
extern lux_data *get_newdata();
extern long      XMAXREQUESTSIZE;
static XSegment  segs[MAXDRAWSIZE];



lux_convert_coordsshort(win,fx,fy,x,y,n)
Window   win;
float   *fx,*fy;
short   *x, *y;
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
	    x[i] = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
			  current->win.xorg)*current->win.xresizefactor + 0.5);
	  else 
	    x[i] = (short)(((log10(fx[i]) - logxmin)*current->win.xfactor +
			  current->win.xorg)*current->win.xresizefactor + 0.5);
	  if ( current->win.lny == 0)
	    y[i+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			    current->win.yfactor + current->win.yorg)*
			   current->win.yresizefactor + 0.5);
	  else
	    y[i+1] = (short)((current->win.ysize - (log10(fy[i+1])-logymin)*
			    current->win.yfactor + current->win.yorg)*
			    current->win.yresizefactor + 0.5);
      }  
    else if (x == y) 
      for(i=0;i<n;i++) {
	  if ( current->win.lnx == 0)
	    x[i*2] = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
			    current->win.xorg)*current->win.xresizefactor + 0.5);
	  else 
	    x[i*2] = (short)(((log10(fx[i]) - logxmin)*current->win.xfactor +
			    current->win.xorg)*current->win.xresizefactor + 0.5);    
	  if ( current->win.lny == 0)
	    x[i*2+1] = (short)((current->win.ysize - (fy[i]-current->win.ymin)*
			      current->win.yfactor + current->win.yorg)*
			     current->win.yresizefactor + 0.5);
	  else
	    x[i*2+1] = (short)((current->win.ysize - (log10(fy[i])-logymin)*
			      current->win.yfactor + current->win.yorg)*
			     current->win.yresizefactor + 0.5);
	
      }        
    else 
      for(i=0;i<n;i++) {
	  if ( current->win.lnx == 0)
	    x[i] = (short)(((fx[i]-current->win.xmin)*current->win.xfactor+
			  current->win.xorg)*current->win.xresizefactor+0.5);
	  else 
	    x[i] = (short)(((log10(fx[i])-logxmin)*current->win.xfactor+
			  current->win.xorg)*current->win.xresizefactor+0.5);
	  if ( current->win.lny == 0)
	    y[i] = (short)((current->win.ysize-(fy[i]-current->win.ymin)*
			  current->win.yfactor+current->win.yorg)*
			 current->win.yresizefactor+0.5);
	  else
	    y[i] = (short)((current->win.ysize-(log10(fy[i])-logymin)*
			  current->win.yfactor+current->win.yorg)*
			 current->win.yresizefactor+0.5);
      }
    
}

lux_convert_coordshort(win,fx,fy,x,y)
Window   win;
float    fx,fy;
short     *x, *y;
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
      *x = (short)(((fx-current->win.xmin)*current->win.xfactor+
		  current->win.xorg)*current->win.xresizefactor+0.5);
    else
      *x = (short)(((log10(fx)-logxmin)*current->win.xfactor+
		  current->win.xorg)*current->win.xresizefactor+0.5);
    if ( current->win.lny == 0)
      *y = (short)((current->win.ysize-(fy-current->win.ymin)*
		  current->win.yfactor+current->win.yorg)*
		 current->win.yresizefactor+0.5);
    else
      *y = (short)((current->win.ysize-(log10(fy)-logymin)*
		  current->win.yfactor+current->win.yorg)*
		 current->win.yresizefactor+0.5);
      
}      

/*******
lux_convert_coordsshort(win,fx,fy,x,y,n)  // Support 3 kinds of format //
Window   win;
float   *fx,*fy;
short   *x, *y;
unsigned int n;
{
    register int i;
    register lux_wins *current;

    current = get_currentwin(win);

    if (fx == fy) // x == y too //
      for(i=0;i<2*n;i+=2) {
	x[i]   = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	y[i+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5);
      }  
    else if (x == y) 
      for(i=0;i<n;i++) {
	x[i*2] = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	x[i*2+1] = (short)((current->win.ysize - (fy[i]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5);
      }        
    else 
      for(i=0;i<n;i++) {
	x[i] = (short)(((fx[i]-current->win.xmin)*current->win.xfactor+
		      current->win.xorg)*current->win.xresizefactor+0.5);
	y[i] = (short)((current->win.ysize-(fy[i]-current->win.ymin)*
		      current->win.yfactor+current->win.yorg)*
		      current->win.yresizefactor+0.5);
      }
}

lux_convert_coordshort(win,fx,fy,x,y)
Window   win;
float    fx,fy;
short     *x, *y;
{
    register long i;
    register lux_wins *current;

    current = get_currentwin(win);

    *x = (short)(((fx-current->win.xmin)*current->win.xfactor+
		 current->win.xorg)*current->win.xresizefactor+0.5);
    *y = (short)((current->win.ysize-(fy-current->win.ymin)*
		current->win.yfactor+current->win.yorg)*
		current->win.yresizefactor+0.5);
}      
*******/

lux_convert_coordsshort_rects(win,fx,fy,w,h,x,y,n)
Window   win;
float   *fx,*fy;
float    w,h;
short   *x, *y;
unsigned int n;
{
    register int i;
    register lux_wins *current;
    short    width, height;

    current = get_currentwin(win);

    width  = (short)(w*current->win.xfactor*current->win.xresizefactor + 0.5);
    height = (short)(h*current->win.yfactor*current->win.yresizefactor + 0.5);

    if ((fx == fy) && ((w <0.0) || (h<0.0))) /* x == y too */
      for(i=0;i<4*n;i+=4) {
	x[i+2] = (short)(fx[i+2]*current->win.xfactor*current->win.xresizefactor + 0.5);
	y[i+3] = (short)(fy[i+3]*current->win.yfactor*current->win.yresizefactor + 0.5);
	x[i]   = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	y[i+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - y[i+3];
      }  
    else if (fx == fy) /* x == y too */
      for(i=0;i<2*n;i+=2) {
	x[i*2]   = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	x[i*2+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - height;
	x[i*2+2] = width;
	x[i*2+3] = height;
      }  
    else if (x == y) 
      for(i=0;i<n;i++) {
	x[i*4] = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	x[i*4+1] = (short)((current->win.ysize - (fy[i]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - height;
	x[i*4+2] = width;
	x[i*4+3] = height;
      }        
}

lux_convert_coord_rect(win,fx,fy,fw,fh,ix,iy,iw,ih)
Window   win;
float   fx,fy,fw,fh;
int    *ix, *iy, *iw, *ih;
{
    register lux_wins *current;

    current = get_currentwin(win);

    *iw  = (int)(fw*current->win.xfactor*current->win.xresizefactor + 0.5);
    *ih  = (int)(fh*current->win.yfactor*current->win.yresizefactor + 0.5);

    *ix  = (int)(((fx - current->win.xmin)*current->win.xfactor +
		  current->win.xorg)*current->win.xresizefactor + 0.5);
    *iy  = (int)((current->win.ysize - (fy-current->win.ymin)*
		  current->win.yfactor + current->win.yorg)*
		  current->win.yresizefactor + 0.5) - *ih;
}

lux_convert_coordsshort_arcs(win,fx,fy,w,h,a1,a2,x,y,n)
Window   win;
float   *fx,*fy;
float    w,h;
float    a1,a2;
short   *x, *y;
unsigned int n;
{
    register int i;
    register lux_wins *current;
    short    width, height;

    current = get_currentwin(win);

    width  = (short)(w*current->win.xfactor*current->win.xresizefactor + 0.5);
    height = (short)(h*current->win.yfactor*current->win.yresizefactor + 0.5);

    a1 = (short)(a1*64.0 + 0.5);
    a2 = (short)(a2*64.0 + 0.5);

    if ((fx==fy) && ((w<0.0)||(h<0.0))) /* x == y too */
      for(i=0;i<6*n;i+=6) {
	x[i+2] = (short)(fx[i+2]*current->win.xfactor*current->win.xresizefactor + 0.5);
	y[i+3] = (short)(fy[i+3]*current->win.yfactor*current->win.yresizefactor + 0.5);
	x[i]   = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	y[i+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - y[i+3];
	x[i+4] = (short)(fx[i+4]*64.0 + 0.5);
	y[i+5] = (short)(fy[i+5]*64.0 + 0.5);
      }  
    else if (fx == fy) /* x == y too */
      for(i=0;i<2*n;i+=2) {
	x[i*3]   = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	x[i*3+1] = (short)((current->win.ysize - (fy[i+1]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - height;
	x[i*3+2] = width;
	x[i*3+3] = height;
	x[i*3+4] = a1;
	x[i*3+5] = a2;
      }  
    else if (x == y) 
      for(i=0;i<n;i++) {
	x[i*6] = (short)(((fx[i] - current->win.xmin)*current->win.xfactor +
		      current->win.xorg)*current->win.xresizefactor + 0.5);
	x[i*6+1] = (short)((current->win.ysize - (fy[i]-current->win.ymin)*
			  current->win.yfactor + current->win.yorg)*
		      current->win.yresizefactor + 0.5) - height;
	x[i*6+2] = width;
	x[i*6+3] = height;
	x[i*6+4] = a1;
	x[i*6+5] = a2;
      }        
}

lux_convert_coord_arc(win,fx,fy,fw,fh,fa1,fa2,ix,iy,iw,ih,ia1,ia2)
Window   win;
float   fx,fy,fw,fh,fa1,fa2;
int    *ix, *iy, *iw, *ih, *ia1, *ia2;
{
    register lux_wins *current;

    current = get_currentwin(win);

    *iw  = (int)(fw*current->win.xfactor*current->win.xresizefactor + 0.5);
    *ih  = (int)(fh*current->win.yfactor*current->win.yresizefactor + 0.5);

    *ix  = (int)(((fx - current->win.xmin)*current->win.xfactor +
		  current->win.xorg)*current->win.xresizefactor + 0.5);
    *iy  = (int)((current->win.ysize - (fy-current->win.ymin)*
		  current->win.yfactor + current->win.yorg)*
		  current->win.yresizefactor + 0.5) - *ih;
    *ia1 = (int)(fa1*64.0 + 0.5);
    *ia2 = (int)(fa2*64.0 + 0.5);

}

lux_draw_pointsf(win,fx,fy,n,mode)  /* 2 stand for CoordModePrevious */
Window   win;
float   *fx, *fy;
unsigned int n;
int      mode;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int size;
    int      m;

    current = get_currentwin(win);

    size=(MAXDRAWSIZE*2>(XMAXREQUESTSIZE-3))?(XMAXREQUESTSIZE-3):MAXDRAWSIZE*2;

    if ((mode==2)&&(n>size)) {
	fprintf(stderr,"ERR: Two many points in draw points with mode 2\n     Maximum points you can draw at one time is %u \n", size);
	return(0);
    }


    if (mode == 2) m = CoordModePrevious;
    else           m = CoordModeOrigin;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	if (fx == fy) 
	 lux_draw_pointsf(win,fx+i*2,fy+i*2,size,mode);
	else
	 lux_draw_pointsf(win,fx+i,fy+i,size,mode);
      }
      if (fx == fy) 
	lux_draw_pointsf(win,fx+i*2,fy+i*2,j,mode);
      else
	lux_draw_pointsf(win,fx+i,fy+i,j,mode);
      return(1);
    }

    if ( current->win.discard_flag ) {
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XDrawPoints(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawPoints(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
      }
      return(1);
    }
    
    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XDrawPoints(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawPoints(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw pointsf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(2*sizeof(int)+2*n*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if ((n <= 10) || (mode == 2)) {
	  lux_convert_coordsshort(win, fx, fy, segs, segs, n);

	  if (current->win.pixmap)
	    XDrawPoints(current->win.display, current->win.pixmap,
		       current->win.gc, (XPoint*)segs, n, m);
	  if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	      XDrawPoints(current->win.display, current->win.window,
			 current->win.gc, (XPoint*)segs, n, m);
	      XFlush(current->win.display);
	  }
	  fprintf(stderr, "\n Not enough memory in draw_pointsf with mode %d\n",mode);
	  newdata->type = NOOP;
	  return(0);
      }
      else {
	j = (int)(n/2);
	if ( fx == fy ) {
	  lux_draw_pointsf(win,fx,fy,j,mode);
	  lux_draw_pointsf(win,fx+j*2,fy+j*2,n-j,mode);
	  return(1);
	}
	else {
	  lux_draw_pointsf(win,fx,fy,j,mode);
	  lux_draw_pointsf(win,fx+j,fy+j,n-j,mode);
	  return(1);
	}
      }
    }

    lux_convert_coordsshort(win, fx, fy, segs, segs, n);

    if (current->win.pixmap)
      XDrawPoints(current->win.display, current->win.pixmap,
		 current->win.gc, (XPoint*)segs, n, m);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawPoints(current->win.display, current->win.window,
		   current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
    }

    newdata->type = DRAW_POINTS_F;

    newdata->data.u[0] = n;
    newdata->data.i[1] = m;
    ptr = (float *)(&(newdata->data.i[2]));
    
    if (fx == fy) 
      for(i=0;i<2*n;i++) ptr[i] = fx[i];
    else
      for(i=0;i<2*n;i+=2) {
	ptr[i]   = fx[i/2];
	ptr[i+1] = fy[i/2];
      }
    return(1);
}

lux_redraw_pointsf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int  m; unsigned int n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    m = (current->win.currentdata)->data.i[1];
    ptrf = (float *)(&((current->win.currentdata)->data.i[2]));

    lux_convert_coordsshort(win, ptrf, ptrf, (short *)segs, (short *)segs, n);

    if (current->win.pixmap) {
      XDrawPoints(current->win.display, current->win.pixmap,
		    redrawgc, (XPoint*)segs, n, m);
    }
     if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawPoints(current->win.display, current->win.window,
		    redrawgc, (XPoint*)segs, n, m);
    }
    return(1);
}

lux_draw_linesf(win,fx,fy,n,mode)
Window   win;
float   *fx, *fy;
unsigned int n;
int      mode;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    int      size, m;

    current = get_currentwin(win);

    size=(MAXDRAWSIZE*2>(XMAXREQUESTSIZE-3))?(XMAXREQUESTSIZE-3):MAXDRAWSIZE*2;

    if ((mode==2)&&(n>size)) {
	fprintf(stderr,"luxwins: Two many points in draw lines with mode 2\n     Maximum points you can draw at one time is %u \n", size);
	return(0);
    }

    if (mode == 2) m = CoordModePrevious;
    else           m = CoordModeOrigin;

    if (n > size) {
      int m = 0;
      j = n+n/size;
      for(i=0;i<n+n/size-size;i+=size) {
	j = j - size;
	if (fx == fy) 
	  /* Make sure that after split the lines are still connected */
	 lux_draw_linesf(win,fx+(i*2-m*2),fy+(i*2-m*2),size,mode);
	else
	 lux_draw_linesf(win,fx+(i-m),fy+(i-m),size,mode);
	m = m + 1;
      }
      if (fx == fy)
	lux_draw_linesf(win,fx+i*2-m*2,fy+i*2-m*2,j,mode);
      else
	lux_draw_linesf(win,fx+i-m,fy+i-m,j,mode);
      return(1);
    }

    if ( current->win.discard_flag ) {
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XDrawLines(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawLines(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XDrawLines(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawLines(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw linesf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(2*sizeof(int)+2*n*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10 || mode == 2) {
	  lux_convert_coordsshort(win, fx, fy, segs, segs, n);

	  if (current->win.pixmap)
	    XDrawLines(current->win.display, current->win.pixmap,
		       current->win.gc, (XPoint*)segs, n, m);
	  if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	      XDrawLines(current->win.display, current->win.window,
			 current->win.gc, (XPoint*)segs, n, m);
	      XFlush(current->win.display);
	  }
	  fprintf(stderr, "\n Not enough memory in draw_linesf with mode %d\n", mode);
	  newdata->type = NOOP;
	  return(0);
      }
      else {
	j = (int)(n/2);
	if ( fx == fy ) {
	  lux_draw_linesf(win,fx,fy,j,mode);
	  lux_draw_linesf(win,fx+j*2-2,fy+j*2-2,n-j+1,mode);
	  return(1);
	}
	else {
	  lux_draw_linesf(win,fx,fy,j,mode);
	  lux_draw_linesf(win,fx+j-1,fy+j-1,n-j+1,mode);
	  return(1);
	}
      }
    }

    lux_convert_coordsshort(win, fx, fy, segs, segs, n);

    if (current->win.pixmap)
      XDrawLines(current->win.display, current->win.pixmap,
		 current->win.gc,(XPoint*)segs, n, m);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawLines(current->win.display, current->win.window,
		   current->win.gc, (XPoint*)segs, n, m);
	XFlush(current->win.display);
    }

    newdata->type = DRAW_LINES_F;

    newdata->data.u[0] = n;
    newdata->data.i[1] = m;
    ptr = (float *)(&(newdata->data.i[2]));
    
    if (fx == fy) 
      for(i=0;i<2*n;i++) ptr[i] = fx[i];
    else
      for(i=0;i<2*n;i+=2) {
	ptr[i]   = fx[i/2];
	ptr[i+1] = fy[i/2];
      }
    return(1);
}

lux_redraw_linesf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int  m; unsigned int n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    m = (current->win.currentdata)->data.i[1];
    ptrf = (float *)(&((current->win.currentdata)->data.i[2]));

    lux_convert_coordsshort(win, ptrf, ptrf, (short *)segs, (short *)segs, n);

    if (current->win.pixmap) {
      XDrawLines(current->win.display, current->win.pixmap,
		    redrawgc, (XPoint*)segs, n, m);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawLines(current->win.display, current->win.window,
		    redrawgc, (XPoint*)segs, n, m);
    }
    return(1);
}

lux_draw_segmentsf(win,fx1,fy1,fx2,fy2,n)
Window  win;
float   *fx1, *fy1, *fx2, *fy2;
unsigned int n;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int size;
   
    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3)/2)?(XMAXREQUESTSIZE-3)/2:MAXDRAWSIZE;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	if (fx1 == fy1) 
	  lux_draw_segmentsf(win,fx1+i*4,fy1+i*4,fx2+i*4,fy2+i*4,size);
	else if (fx1 == fx2)
	  lux_draw_segmentsf(win,fx1+i*2,fy1+i*2,fx2+i*2,fy2+i*2,size);
	else
	  lux_draw_segmentsf(win,fx1+i,fy1+i,fx2+i,fy2+i,size);
      }
      if (fx1 == fy1) 
	lux_draw_segmentsf(win,fx1+i*4,fy1+i*4,fx2+i*4,fy2+i*4,j);
      else if (fx1 == fx2)
	lux_draw_segmentsf(win,fx1+i*2,fy1+i*2,fx2+i*2,fy2+i*2,j);
      else
	lux_draw_segmentsf(win,fx1+i,fy1+i,fx2+i,fy2+i,j);
      return(1);
    }

    if ( current->win.discard_flag ) {
      if (fx1 == fy1)
	lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
      else if (fx1 == fx2) 
	lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
      else 
	for(i=0;i<n;i++) {
	  lux_convert_coordshort(win, fx1[i], fy1[i], &segs[i].x1, &segs[i].y1);
	  lux_convert_coordshort(win, fx2[i], fy2[i], &segs[i].x2, &segs[i].y2);
	}

      if (current->win.pixmap)
	XDrawSegments(current->win.display, current->win.pixmap,
		      current->win.gc, (XSegment*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawSegments(current->win.display, current->win.window,
		      current->win.gc, (XSegment*)segs, n);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      if (fx1 == fy1)
	lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
      else if (fx1 == fx2) 
	lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
      else 
	for(i=0;i<n;i++) {
	  lux_convert_coordshort(win, fx1[i], fy1[i], &segs[i].x1, &segs[i].y1);
	  lux_convert_coordshort(win, fx2[i], fy2[i], &segs[i].x2, &segs[i].y2);
	}

      if (current->win.pixmap)
	XDrawSegments(current->win.display, current->win.pixmap,
		      current->win.gc, (XSegment*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawSegments(current->win.display, current->win.window,
		      current->win.gc, (XSegment*)segs, n);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw segmentsf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(sizeof(int)+4*n*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10) {
	if (fx1 == fy1)
	  lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
	else if (fx1 == fx2) 
	  lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
	else 
	  for(i=0;i<n;i++) {
	    lux_convert_coordshort(win, fx1[i],fy1[i], &segs[i].x1,&segs[i].y1);
	    lux_convert_coordshort(win, fx2[i],fy2[i], &segs[i].x2,&segs[i].y2);
	  }

	if (current->win.pixmap)
	  XDrawSegments(current->win.display, current->win.pixmap,
			current->win.gc, (XSegment*)segs, n);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawSegments(current->win.display, current->win.window,
			current->win.gc, (XSegment*)segs, n);
	  XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in draw_segmentsf\n");
	newdata->type = NOOP;
	return(0);
      }
      else {
	j = (int)(n/2);
	if ( fx1 == fy1 ) {
	  lux_draw_segmentsf(win,fx1,fy1,fx2,fy2,j);
	  lux_draw_segmentsf(win,fx1+j*4,fy1+j*4,fx2+j*4,fy2+j*4,n-j);
	  return(1);
	}
	else if (fx1 == fx2) {
	  lux_draw_segmentsf(win,fx1,fy1,fx2,fy2,j);
	  lux_draw_segmentsf(win,fx1+j*2,fy1+j*2,fx2+j*2,fy2+j*2,n-j);
	}
	else {
	  lux_draw_segmentsf(win,fx1,fy1,fx2,fy2,j);
	  lux_draw_segmentsf(win,fx1+j,fy1+j,fx2+j,fy2+j,n-j);
	  return(1);
	}
      }
    }

    if (fx1 == fy1)
      lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
    else if (fx1 == fx2) 
      lux_convert_coordsshort(win, fx1, fy1, segs, segs, 2*n);
    else 
      for(i=0;i<n;i++) {
	lux_convert_coordshort(win, fx1[i],fy1[i], &segs[i].x1,&segs[i].y1);
	lux_convert_coordshort(win, fx2[i],fy2[i], &segs[i].x2,&segs[i].y2);
      }

    if (current->win.pixmap)
      XDrawSegments(current->win.display, current->win.pixmap,
		    current->win.gc, (XSegment*)segs, n);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawSegments(current->win.display, current->win.window,
		    current->win.gc, (XSegment*)segs, n);
      XFlush(current->win.display);
    }

    newdata->type = DRAW_SEGMENTS_F;

    newdata->data.u[0] = n;
    ptr = (float *)(&(newdata->data.i[1]));
    
    if (fx1 == fy1) 
      for(i=0;i<4*n;i++) ptr[i] = fx1[i];
    else if (fx1 == fx2)
      for(i=0;i<4*n;i+=2) {
	ptr[i]   = fx1[i/2];
	ptr[i+1] = fy1[i/2];
      }
    else
      for(i=0;i<4*n;i+=4) {
	ptr[i]   = fx1[i/4];
	ptr[i+1] = fy1[i/4];
	ptr[i+2] = fx2[i/4];
	ptr[i+3] = fy2[i/4];
      }
    return(1);
}

lux_redraw_segmentsf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    unsigned int n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    ptrf = (float *)(&((current->win.currentdata)->data.i[1]));

    lux_convert_coordsshort(win, ptrf, ptrf, (short *)segs, (short *)segs, 2*n);

    if (current->win.pixmap) {
      XDrawSegments(current->win.display, current->win.pixmap,
		    redrawgc, (XSegment*)segs, n);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawSegments(current->win.display, current->win.window,
		    redrawgc, (XSegment*)segs, n);
    }
    return(1);
}

lux_draw_rectanglesf(win,fx1,fy1,w,h,n)
Window  win;
float   *fx1, *fy1, w, h;
unsigned int n;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int    size, times;
   
    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3)/2)?(XMAXREQUESTSIZE-3)/2:MAXDRAWSIZE;


    if ((fx1 == fy1)&&(w<0.0||h<0.0)) times = 4;
    else if (fx1 == fy1) times = 2;
    else times = 1;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	lux_draw_rectanglesf(win,fx1+i*times,fy1+i*times,w,h,size);
      }
      lux_draw_rectanglesf(win,fx1+i*times,fy1+i*times,w,h,j);
      return(1);
    }

    if ( current->win.discard_flag ) {
      lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

      if (current->win.pixmap)
	XDrawRectangles(current->win.display, current->win.pixmap,
		      current->win.gc, (XRectangle*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawRectangles(current->win.display, current->win.window,
		      current->win.gc, (XRectangle*)segs, n);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

      if (current->win.pixmap)
	XDrawRectangles(current->win.display, current->win.pixmap,
		      current->win.gc, (XRectangle*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawRectangles(current->win.display, current->win.window,
		      current->win.gc, (XRectangle*)segs, n);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw rectanglesf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(sizeof(int)+
				    (n*((times==4)?4:2)+2)*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10) {
	lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);
	if (current->win.pixmap)
	  XDrawRectangles(current->win.display, current->win.pixmap,
			current->win.gc, (XRectangle*)segs, n);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawRectangles(current->win.display, current->win.window,
			current->win.gc, (XRectangle*)segs, n);
	  XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in draw_rectanglesf\n");
	newdata->type = NOOP;
	return(0);
      }
      else {
	j = (int)(n/2);
	lux_draw_rectanglesf(win,fx1,fy1,w,h,j);
	lux_draw_rectanglesf(win,fx1+j*times,fy1+j*times,w,h,n-j);
	return(1);
      }
    }

    lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

    if (current->win.pixmap)
      XDrawRectangles(current->win.display, current->win.pixmap,
		    current->win.gc, (XRectangle*)segs, n);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawRectangles(current->win.display, current->win.window,
		    current->win.gc, (XRectangle*)segs, n);
      XFlush(current->win.display);
    }

    newdata->type = DRAW_RECTANGLES_F;

    newdata->data.u[0] = n;
    ptr = (float *)(&(newdata->data.i[1]));
    ptr[0] = w;
    ptr[1] = h;
    
    if (fx1 == fy1) 
      for(i=0;i<times*n;i++) ptr[i+2] = fx1[i];
    else 
      for(i=0;i<n;i++) {
	ptr[i*2+2]   = fx1[i];
	ptr[i*2+1+2] = fy1[i];
      }
    return(1);
}

lux_redraw_rectanglesf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    float w, h;
    unsigned int  n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    ptrf = (float *)(&((current->win.currentdata)->data.i[1]));
    w = ptrf[0];
    h = ptrf[1];

    lux_convert_coordsshort_rects(win, &ptrf[2], &ptrf[2], w, h, 
				  (short *)segs, (short *)segs, n);

    if (current->win.pixmap) {
      XDrawRectangles(current->win.display, current->win.pixmap,
		    redrawgc, (XRectangle*)segs, n);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawRectangles(current->win.display, current->win.window,
		    redrawgc, (XRectangle*)segs, n);
    }
    return(1);
}

lux_draw_rectanglef(win,fx,fy,w,h)
Window  win;
float   fx, fy, w, h;
{

    register lux_wins *current;
    register lux_data *newdata;
    register float *ptr;
    int      x, y, width, height;
   
    current = get_currentwin(win);
    
    if ( current->win.discard_flag ) {
      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XDrawRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 

      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XDrawRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw rectanglef because of not enough memory\n"); 
      return(0);
    }

    newdata->data.f = (float *)malloc(4*sizeof(float));

    if (newdata->data.i == (int *)NULL) {

      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XDrawRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "\n Not enough memory in draw_rectanglef\n");
      newdata->type = NOOP;
      return(0);
    }

    lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

    if (current->win.pixmap)
      XDrawRectangle(current->win.display, current->win.pixmap,
		     current->win.gc, x,y, width, height);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawRectangle(current->win.display, current->win.window,
		     current->win.gc, x,y, width, height);
      XFlush(current->win.display);
    }

    newdata->type = DRAW_RECTANGLE_F;

    ptr = newdata->data.f;
    ptr[0] = fx;
    ptr[1] = fy;
    ptr[2] = w;
    ptr[3] = h;

    return(1);
}

lux_redraw_rectanglef(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int x, y, width, height;


    current = get_currentwin(win);

    ptrf = (current->win.currentdata)->data.f;

    lux_convert_coord_rect(win, ptrf[0], ptrf[1], ptrf[2], ptrf[3], 
				&x, &y, &width, &height);

    if (current->win.pixmap) {
      XDrawRectangle(current->win.display, current->win.pixmap,
		     redrawgc, x, y, width, height);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawRectangle(current->win.display, current->win.window,
		     redrawgc, x, y, width, height);
    }
    return(1);
}

lux_draw_arcsf(win,fx1,fy1,w,h,a1,a2,n)
Window  win;
float   *fx1, *fy1, w, h,a1,a2;
unsigned int n;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int    size, times;
   
    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3)/3)?(XMAXREQUESTSIZE-3)/3:MAXDRAWSIZE;


    if ((fx1 == fy1)&&(w<0.0||h<0.0)) times = 6;
    else if (fx1 == fy1) times = 2;
    else times = 1;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	lux_draw_arcsf(win,fx1+i*times,fy1+i*times,w,h,a1,a2,size);
      }
      lux_draw_arcsf(win,fx1+i*times,fy1+i*times,w,h,a1,a2,j);
      return(1);
    }

    if ( current->win.discard_flag ) {
      lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

      if (current->win.pixmap)
	XDrawArcs(current->win.display, current->win.pixmap,
		      current->win.gc, (XArc*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawArcs(current->win.display, current->win.window,
		      current->win.gc, (XArc*)segs, n);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

      if (current->win.pixmap)
	XDrawArcs(current->win.display, current->win.pixmap,
		      current->win.gc, (XArc*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawArcs(current->win.display, current->win.window,
		      current->win.gc, (XArc*)segs, n);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw arcsf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(sizeof(int)+
				    (n*((times==6)?6:2)+4)*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10) {
	lux_convert_coordsshort_arcs(win, fx1, fy1, w, h,a1, a2, segs, segs, n);
	if (current->win.pixmap)
	  XDrawArcs(current->win.display, current->win.pixmap,
			current->win.gc, (XArc*)segs, n);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawArcs(current->win.display, current->win.window,
			current->win.gc, (XArc*)segs, n);
	  XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in draw_arcsf\n");
	newdata->type = NOOP;
	return(0);
      }
      else {
	j = (int)(n/2);
	lux_draw_arcsf(win,fx1,fy1,w,h,a1,a2,j);
	lux_draw_arcsf(win,fx1+j*times,fy1+j*times,w,h,a1,a2,n-j);
	return(1);
      }
    }

    lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

    if (current->win.pixmap)
      XDrawArcs(current->win.display, current->win.pixmap,
		    current->win.gc, (XArc*)segs, n);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XDrawArcs(current->win.display, current->win.window,
		    current->win.gc, (XArc*)segs, n);
      XFlush(current->win.display);
    }

    newdata->type = DRAW_ARCS_F;

    newdata->data.u[0] = n;
    ptr = (float *)(&(newdata->data.i[1]));
    ptr[0] = w;
    ptr[1] = h;
    ptr[2] = a1;
    ptr[3] = a2;
    
    if (fx1 == fy1) 
      for(i=0;i<times*n;i++) ptr[i+4] = fx1[i];
    else 
      for(i=0;i<n;i++) {
	ptr[i*2+4]   = fx1[i];
	ptr[i*2+1+4] = fy1[i];
      }
    return(1);
}

lux_redraw_arcsf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    unsigned int  n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    ptrf = (float *)(&((current->win.currentdata)->data.i[1]));

    lux_convert_coordsshort_arcs(win, &ptrf[2], &ptrf[2], ptrf[0], ptrf[1], 
				ptrf[2],ptrf[3],(short *)segs,(short *)segs, n);

    if (current->win.pixmap) {
      XDrawArcs(current->win.display, current->win.pixmap,
		redrawgc, (XArc*)segs, n);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawArcs(current->win.display, current->win.window,
		redrawgc, (XArc*)segs, n);
    }
    return(1);
}

lux_draw_arcf(win,fx,fy,w,h,a1,a2)
Window  win;
float   fx, fy, w, h, a1,a2;
{

    register lux_wins *current;
    register lux_data *newdata;
    register float *ptr;
    int      x, y, width, height, ia1, ia2;
   
    current = get_currentwin(win);
    
    if ( current->win.discard_flag ) {
      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XDrawArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 

      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XDrawArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "Can't draw arcf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.f = (float *)malloc(6*sizeof(float));

    if (newdata->data.i == (int *)NULL) {

      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XDrawArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XDrawArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "\n Not enough memory in draw_arcf\n");
      newdata->type = NOOP;
      return(0);
    }


    lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

    if (current->win.pixmap)
      XDrawArc(current->win.display, current->win.pixmap,
	       current->win.gc, x,y, width, height,ia1,ia2);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XDrawArc(current->win.display, current->win.window,
		 current->win.gc, x,y, width, height,ia1,ia2);
	XFlush(current->win.display);
    }

    newdata->type = DRAW_ARC_F;

    ptr = newdata->data.f;
    ptr[0] = fx;
    ptr[1] = fy;
    ptr[2] = w;
    ptr[3] = h;
    ptr[4] = a1;
    ptr[5] = a2;

    return(1);
}

lux_redraw_arcf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int x, y, width, height, a1, a2;


    current = get_currentwin(win);

    ptrf = (current->win.currentdata)->data.f;

    lux_convert_coord_arc(win, ptrf[0], ptrf[1], ptrf[2], ptrf[3], ptrf[4],
				ptrf[5], &x, &y, &width, &height, &a1, &a2);

    if (current->win.pixmap) {
      XDrawArc(current->win.display, current->win.pixmap,
		     redrawgc, x, y, width, height, a1, a2);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XDrawArc(current->win.display, current->win.window,
		     redrawgc, x, y, width, height, a1, a2);
    }
    return(1);
}


lux_fill_rectanglesf(win,fx1,fy1,w,h,n)
Window  win;
float   *fx1, *fy1, w, h;
unsigned int n;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int    size, times;
   
    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3)/2)?(XMAXREQUESTSIZE-3)/2:MAXDRAWSIZE;


    if ((fx1 == fy1)&&(w<0.0||h<0.0)) times = 4;
    else if (fx1 == fy1) times = 2;
    else times = 1;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	lux_fill_rectanglesf(win,fx1+i*times,fy1+i*times,w,h,size);
      }
      lux_fill_rectanglesf(win,fx1+i*times,fy1+i*times,w,h,j);
      return(1);
    }

    if ( current->win.discard_flag ) { 
      lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

      if (current->win.pixmap)
	XFillRectangles(current->win.display, current->win.pixmap,
		      current->win.gc, (XRectangle*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillRectangles(current->win.display, current->win.window,
		      current->win.gc, (XRectangle*)segs, n);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

      if (current->win.pixmap)
	XFillRectangles(current->win.display, current->win.pixmap,
		      current->win.gc, (XRectangle*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillRectangles(current->win.display, current->win.window,
		      current->win.gc, (XRectangle*)segs, n);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't fill rectanglesf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(sizeof(int)+
				    (n*((times==4)?4:2)+2)*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10) {
	lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);
	if (current->win.pixmap)
	  XFillRectangles(current->win.display, current->win.pixmap,
			current->win.gc, (XRectangle*)segs, n);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillRectangles(current->win.display, current->win.window,
			current->win.gc, (XRectangle*)segs, n);
	  XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in fill_rectanglesf\n");
	newdata->type = NOOP;
	return(0);
      }
      else {
	j = (int)(n/2);
	lux_fill_rectanglesf(win,fx1,fy1,w,h,j);
	lux_fill_rectanglesf(win,fx1+j*times,fy1+j*times,w,h,n-j);
	return(1);
      }
    }

    lux_convert_coordsshort_rects(win, fx1, fy1, w, h, segs, segs, n);

    if (current->win.pixmap)
      XFillRectangles(current->win.display, current->win.pixmap,
		    current->win.gc, (XRectangle*)segs, n);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XFillRectangles(current->win.display, current->win.window,
		    current->win.gc, (XRectangle*)segs, n);
      XFlush(current->win.display);
    }

    newdata->type = FILL_RECTANGLES_F;

    newdata->data.u[0] = n;
    ptr = (float *)(&(newdata->data.i[1]));
    ptr[0] = w;
    ptr[1] = h;
    
    if (fx1 == fy1) 
      for(i=0;i<times*n;i++) ptr[i+2] = fx1[i];
    else 
      for(i=0;i<n;i++) {
	ptr[i*2+2]   = fx1[i];
	ptr[i*2+1+2] = fy1[i];
      }
    return(1);
}

lux_refill_rectanglesf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    float w, h;
    unsigned int  n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    ptrf = (float *)(&((current->win.currentdata)->data.i[1]));
    w = ptrf[0];
    h = ptrf[1];

    lux_convert_coordsshort_rects(win, &ptrf[2], &ptrf[2], w, h, 
				  (short *)segs, (short *)segs, n);

    if (current->win.pixmap) {
      XFillRectangles(current->win.display, current->win.pixmap,
		    redrawgc, (XRectangle*)segs, n);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XFillRectangles(current->win.display, current->win.window,
		    redrawgc, (XRectangle*)segs, n);
    }
    return(1);
}

lux_fill_rectanglef(win,fx,fy,w,h)
Window  win;
float   fx, fy, w, h;
{

    register lux_wins *current;
    register lux_data *newdata;
    register float *ptr;
    int      x, y, width, height;
   
    current = get_currentwin(win);
    
    if ( current->win.discard_flag ) { 
      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XFillRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 

      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XFillRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "Can't fill rectanglef because of not enough memory\n"); 
      return(0);
    }

    newdata->data.f = (float *)malloc(4*sizeof(float));

    if (newdata->data.i == (int *)NULL) {

      lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

      if (current->win.pixmap)
	XFillRectangle(current->win.display, current->win.pixmap,
		       current->win.gc, x,y, width, height);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillRectangle(current->win.display, current->win.window,
			 current->win.gc, x,y, width, height);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "\n Not enough memory in fill_rectanglef\n");
      newdata->type = NOOP;
      return(0);
    }

    lux_convert_coord_rect(win,fx,fy,w,h,&x,&y,&width,&height);

    if (current->win.pixmap)
      XFillRectangle(current->win.display, current->win.pixmap,
		     current->win.gc, x,y, width, height);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XFillRectangle(current->win.display, current->win.window,
		     current->win.gc, x,y, width, height);
      XFlush(current->win.display);
    }

    newdata->type = FILL_RECTANGLE_F;

    ptr = newdata->data.f;
    ptr[0] = fx;
    ptr[1] = fy;
    ptr[2] = w;
    ptr[3] = h;

    return(1);
}

lux_refill_rectanglef(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int x, y, width, height;


    current = get_currentwin(win);

    ptrf = (current->win.currentdata)->data.f;

    lux_convert_coord_rect(win, ptrf[0], ptrf[1], ptrf[2], ptrf[3], 
				&x, &y, &width, &height);

    if (current->win.pixmap) {
      XFillRectangle(current->win.display, current->win.pixmap,
		     redrawgc, x, y, width, height);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XFillRectangle(current->win.display, current->win.window,
		     redrawgc, x, y, width, height);
    }
    return(1);
}

lux_fill_arcsf(win,fx1,fy1,w,h,a1,a2,n)
Window  win;
float   *fx1, *fy1, w, h,a1,a2;
unsigned int n;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int    size, times;
   
    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3)/3)?(XMAXREQUESTSIZE-3)/3:MAXDRAWSIZE;


    if ((fx1 == fy1)&&(w<0.0||h<0.0)) times = 6;
    else if (fx1 == fy1) times = 2;
    else times = 1;

    if (n > size) {
      j = n;
      for(i=0;i<n-size;i+=size) {
	j = j - size;
	lux_fill_arcsf(win,fx1+i*times,fy1+i*times,w,h,a1,a2,size);
      }
      lux_fill_arcsf(win,fx1+i*times,fy1+i*times,w,h,a1,a2,j);
      return(1);
    }

    if ( current->win.discard_flag ) { 
      lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

      if (current->win.pixmap)
	XFillArcs(current->win.display, current->win.pixmap,
		      current->win.gc, (XArc*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillArcs(current->win.display, current->win.window,
		      current->win.gc, (XArc*)segs, n);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

      if (current->win.pixmap)
	XFillArcs(current->win.display, current->win.pixmap,
		      current->win.gc, (XArc*)segs, n);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillArcs(current->win.display, current->win.window,
		      current->win.gc, (XArc*)segs, n);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't fill arcsf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(sizeof(int)+
				    (n*((times==6)?6:2)+4)*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
      if (n <= 10) {
	lux_convert_coordsshort_arcs(win, fx1, fy1, w, h,a1, a2, segs, segs, n);
	if (current->win.pixmap)
	  XFillArcs(current->win.display, current->win.pixmap,
			current->win.gc, (XArc*)segs, n);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillArcs(current->win.display, current->win.window,
			current->win.gc, (XArc*)segs, n);
	  XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in fill_arcsf\n");
	newdata->type = NOOP;
	return(0);
      }
      else {
	j = (int)(n/2);
	lux_fill_arcsf(win,fx1,fy1,w,h,a1,a2,j);
	lux_fill_arcsf(win,fx1+j*times,fy1+j*times,w,h,a1,a2,n-j);
	return(1);
      }
    }

    lux_convert_coordsshort_arcs(win, fx1, fy1, w, h, a1, a2, segs, segs, n);

    if (current->win.pixmap)
      XFillArcs(current->win.display, current->win.pixmap,
		    current->win.gc, (XArc*)segs, n);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
      XFillArcs(current->win.display, current->win.window,
		    current->win.gc, (XArc*)segs, n);
      XFlush(current->win.display);
    }

    newdata->type = FILL_ARCS_F;

    newdata->data.u[0] = n;
    ptr = (float *)(&(newdata->data.i[1]));
    ptr[0] = w;
    ptr[1] = h;
    ptr[2] = a1;
    ptr[3] = a2;
    
    if (fx1 == fy1) 
      for(i=0;i<times*n;i++) ptr[i+4] = fx1[i];
    else 
      for(i=0;i<n;i++) {
	ptr[i*2+4]   = fx1[i];
	ptr[i*2+1+4] = fy1[i];
      }
    return(1);
}

lux_refill_arcsf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    unsigned int  n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    ptrf = (float *)(&((current->win.currentdata)->data.i[1]));

    lux_convert_coordsshort_arcs(win, &ptrf[2], &ptrf[2], ptrf[0], ptrf[1], 
				ptrf[2],ptrf[3],(short *)segs,(short *)segs, n);

    if (current->win.pixmap) {
      XFillArcs(current->win.display, current->win.pixmap,
		redrawgc, (XArc*)segs, n);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XFillArcs(current->win.display, current->win.window,
		redrawgc, (XArc*)segs, n);
    }
    return(1);
}

lux_fill_arcf(win,fx,fy,w,h,a1,a2)
Window  win;
float   fx, fy, w, h, a1,a2;
{

    register lux_wins *current;
    register lux_data *newdata;
    register float *ptr;
    int      x, y, width, height, ia1, ia2;
   
    current = get_currentwin(win);
    
    if ( current->win.discard_flag ) { 
      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XFillArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 

      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XFillArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "Can't fill arcf because of not enough memory\n"); 
      return(0);
    }

    newdata->data.f = (float *)malloc(6*sizeof(float));

    if (newdata->data.i == (int *)NULL) {

      lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

      if (current->win.pixmap)
	XFillArc(current->win.display, current->win.pixmap,
		 current->win.gc, x,y, width, height,ia1,ia2);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	  XFillArc(current->win.display, current->win.window,
		   current->win.gc, x,y, width, height,ia1,ia2);
	  XFlush(current->win.display);
      }
      fprintf(stderr, "\n Not enough memory in fill_arcf\n");
      newdata->type = NOOP;
      return(0);
    }


    lux_convert_coord_arc(win,fx,fy,w,h,a1,a2,&x,&y,&width,&height,&ia1,&ia2);

    if (current->win.pixmap)
      XFillArc(current->win.display, current->win.pixmap,
	       current->win.gc, x,y, width, height,ia1,ia2);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillArc(current->win.display, current->win.window,
		 current->win.gc, x,y, width, height,ia1,ia2);
	XFlush(current->win.display);
    }

    newdata->type = FILL_ARC_F;

    ptr = newdata->data.f;
    ptr[0] = fx;
    ptr[1] = fy;
    ptr[2] = w;
    ptr[3] = h;
    ptr[4] = a1;
    ptr[5] = a2;

    return(1);
}

lux_refill_arcf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int x, y, width, height, a1, a2;


    current = get_currentwin(win);

    ptrf = (current->win.currentdata)->data.f;

    lux_convert_coord_arc(win, ptrf[0], ptrf[1], ptrf[2], ptrf[3], ptrf[4],
				ptrf[5], &x, &y, &width, &height, &a1, &a2);

    if (current->win.pixmap) {
      XFillArc(current->win.display, current->win.pixmap,
		     redrawgc, x, y, width, height, a1, a2);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XFillArc(current->win.display, current->win.window,
		     redrawgc, x, y, width, height, a1, a2);
    }
    return(1);
}

lux_fill_polygonf(win,fx,fy,n,mode)  /* 2 stand for CoordModePrevious */
Window   win;
float   *fx, *fy;
unsigned int n;
int      mode;
{

    register lux_wins *current;
    register lux_data *newdata;
    register int i,j;
    register float *ptr;
    unsigned int size;
    int      m;

    current = get_currentwin(win);

    size=(MAXDRAWSIZE>(XMAXREQUESTSIZE-3))?(XMAXREQUESTSIZE-3):MAXDRAWSIZE;

    if (n>size) {
	fprintf(stderr,"ERR: Two many points in fill polygon\n     Maximum points you can draw at one time is %u \n", size);
	return(0);
    }

    if (mode == 2) m = CoordModePrevious;
    else           m = CoordModeOrigin;

    if ( current->win.discard_flag ) { 
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XFillPolygon(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, Complex, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillPolygon(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, Complex, m);
	XFlush(current->win.display);
      }
      return(1);
    }

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL ) { 
      lux_convert_coordsshort(win, fx, fy, segs, segs, n);

      if (current->win.pixmap)
	XFillPolygon(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, Complex, m);
      if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillPolygon(current->win.display, current->win.window,
		      current->win.gc, (XPoint*)segs, n, Complex, m);
	XFlush(current->win.display);
      }
      fprintf(stderr, "Can't fill polygon because of not enough memory\n"); 
      return(0);
    }

    newdata->data.i = (int *)malloc(2*sizeof(int)+2*n*sizeof(float));

    if (newdata->data.i == (int *)NULL) {
	lux_convert_coordsshort(win, fx, fy, segs, segs, n);

	if (current->win.pixmap)
	  XFillPolygon(current->win.display, current->win.pixmap,
		      current->win.gc, (XPoint*)segs, n, Complex, m);
	if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	    XFillPolygon(current->win.display, current->win.window,
			current->win.gc, (XPoint*)segs, n, Complex, m);
	    XFlush(current->win.display);
	}
	fprintf(stderr, "\n Not enough memory in fill polygon\n");
	newdata->type = NOOP;
	return(0);
    }

    lux_convert_coordsshort(win, fx, fy, segs, segs, n);

    if (current->win.pixmap)
      XFillPolygon(current->win.display, current->win.pixmap,
		 current->win.gc, (XPoint*)segs, n, Complex, m);
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL)) {
	XFillPolygon(current->win.display, current->win.window,
		   current->win.gc, (XPoint*)segs, n, Complex, m);
	XFlush(current->win.display);
    }

    newdata->type = FILL_POLYGON_F;

    newdata->data.u[0] = n;
    newdata->data.i[1] = m;
    ptr = (float *)(&(newdata->data.i[2]));
    
    if (fx == fy) 
      for(i=0;i<2*n;i++) ptr[i] = fx[i];
    else
      for(i=0;i<2*n;i+=2) {
	ptr[i]   = fx[i/2];
	ptr[i+1] = fy[i/2];
      }
    return(1);
}

lux_refill_polygonf(win)
Window  win;
{
    register lux_wins *current;
    float *ptrf;
    int  m; unsigned int n;

    current = get_currentwin(win);

    n = (current->win.currentdata)->data.u[0];
    m = (current->win.currentdata)->data.i[1];
    ptrf = (float *)(&((current->win.currentdata)->data.i[2]));

    lux_convert_coordsshort(win, ptrf, ptrf, (short *)segs, (short *)segs, n);

    if (current->win.pixmap) {
      XFillPolygon(current->win.display, current->win.pixmap,
		    redrawgc, (XPoint*)segs, n, Complex, m);
    }
    if (current->win.update_flag || (current->win.pixmap == (Pixmap)NULL))  {
      XFillPolygon(current->win.display, current->win.window,
		    redrawgc, (XPoint*)segs, n, Complex, m);
    }
    return(1);
}

lux_set_window_bgcolor(win,i)
Window win;
unsigned long i;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    XSetWindowBackground(current->win.display,current->win.window,i);

    if ( current->win.discard_flag ) return(1);

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr, "Can't set backgroud because of not enough memory\n");
	return(0);}

    newdata->data.l = (long *)malloc(sizeof(long));

    if (newdata->data.l == (long *)NULL) {
      fprintf(stderr, "\n Not enough memory in set bgcolor\n");
      newdata->type = NOOP;
      return(0);
    }

    newdata->type = SET_WINDOW_BG_COLOR;
    newdata->data.l[0] = i;  
    return(1);
}

lux_reset_window_bgcolor(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    XSetWindowBackground(current->win.display, win,
		  (current->win.currentdata)->data.l[0]);
}

lux_set_bgcolor(win,i)
Window win;
unsigned long i;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    XSetBackground(current->win.display,current->win.gc,i);

    if ( current->win.discard_flag ) return(1);

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr, "Can't set backgroud because of not enough memory\n");
	return(0);}

    newdata->data.l = (long *)malloc(sizeof(long));

    if (newdata->data.l == (long *)NULL) {
      fprintf(stderr, "\n Not enough memory in set bgcolor\n");
      newdata->type = NOOP;
      return(0);
    }

    newdata->type = SET_BG_COLOR;
    newdata->data.l[0] = i;  
    return(1);
}

lux_reset_bgcolor(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    XSetBackground(current->win.display, redrawgc,
		  (current->win.currentdata)->data.l[0]);
}

lux_set_color(win,i)
Window win;
unsigned long i;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    XSetForeground(current->win.display, current->win.gc, i);

    if ( current->win.discard_flag ) return(1);

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr, "Can't set color width because of not enough memory\n");
	return(0);}

    newdata->data.ul = (unsigned long *)malloc(sizeof(long));

    if (newdata->data.ul == (unsigned long *)NULL) {
      fprintf(stderr, "\n Not enough memory in set color\n");
      newdata->type = NOOP;
      return(0);
    }

    newdata->type = SET_COLOR;
    newdata->data.ul[0] = i;  
    return(1);
}

lux_reset_color(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    XSetForeground(current->win.display, redrawgc,
		  (current->win.currentdata)->data.ul[0]);
}

unsigned long  lux_get_fgcolor(win)
Window win;
{
    register lux_wins *current;
    XGCValues values;

    current = get_currentwin(win);

    XGetGCValues(current->win.display, current->win.gc, 
		 GCForeground, &values);
    return(values.foreground);
}

unsigned long  lux_get_bgcolor(win)
Window win;
{
    register lux_wins *current;
    XGCValues values;

    current = get_currentwin(win);

    XGetGCValues(current->win.display, current->win.gc, 
		 GCBackground, &values);
    return(values.background);
}

lux_set_update(win)  /* Automaticly update foreground */
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    current->win.update_flag = 1;

    if ( current->win.discard_flag ) return(1);

    if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr, "Can't set backgroud because of not enough memory\n");
	return(0);}

    newdata->type = SET_UPDATE;
    return(1);
}

lux_reset_update(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    current->win.update_flag = 1;
}

lux_set_discard(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    current->win.discard_flag = 1;
}

lux_set_save(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    current->win.discard_flag = 0;
}

lux_set_noupdate(win)  /* NO automaticly update foreground */
Window win;
{
    register lux_wins *current;
    register lux_data *newdata;

    current = get_currentwin(win);

    current->win.update_flag = 0;

    if ( current->win.discard_flag ) return(1);

   if ( (newdata = get_newdata(current)) == (lux_data *)NULL )
      { fprintf(stderr, "Can't set backgroud because of not enough memory\n");
	return(0);}

    newdata->type = SET_NO_UPDATE;
    return(1);
}

lux_reset_noupdate(win)
Window win;
{
    register lux_wins *current;

    current = get_currentwin(win);

    current->win.update_flag = 0;
}
