/* mcd.c                       Interface for McDraw
 * Biao Lu                     biao@eagle.drexel.edu
 *
 * Modifications by Steve McMillan, 8/96.
 *
 * Note old K&R-style arguments in functions!
 */

#include <stdio.h>
#include <stdlib.h>

#define MAXPOINTS  16000
#define MAXCOLOR   16

#define DEFAULT_CM	0
#define FILE_CM		1

typedef struct _mcdwin {
    unsigned long  win;
    int            id;
    char           name[50];
    int            linewidth;
    unsigned long  fg_color;
    unsigned long  bg_color;
    int		   colormap;
    struct _mcdwin *next;
} mcdwin;

static unsigned long   win;
static mcdwin  *window = NULL;
static int buffer_flag = 1;
static float penx = 1.0,peny = 1.0;
static float x[MAXPOINTS], y[MAXPOINTS];
static unsigned long npoints = 0;
static int   style = -1;  /* 0 for points, 1 for lines */
static float aspect_ratio =  1.0;
static unsigned int width, height;

/* Only 16 colors available now */

static unsigned long color[MAXCOLOR];
static char         *colorname[] = { 
                         "white","black","blue","purple","violet","magenta",
		         "lightblue","gray","cyan","green","limegreen",
			 "yellow","orange","brown","red","pink"};

extern unsigned long  lux_get_bgcolor();
extern unsigned long  lux_get_fgcolor();
extern unsigned long  lux_rgb_pixel();
extern unsigned long  lux_lookup_color();

#define SCREEN_SIZE_FRACTION 0.4
#define SCREEN_XPOS_FRACTION 0.3
#define SCREEN_YPOS_FRACTION 0.05

mcdxinit(aspect, nxpix, nypix, ncolor, ierr)
float *aspect;
int *nxpix, *nypix, *ncolor;
int *ierr;
{
    static int status = 1;
    static int n_color;

    lux_use_mcdraw();

    aspect_ratio = *aspect;

    if (lux_display(&width,&height)) {*ierr = 2;return 2;}

    mcdxopenwin();
    
    if (status) {status = 0; n_color = mcdxcolors();}

    *ncolor = n_color;

    *nxpix = width*SCREEN_SIZE_FRACTION;
    *nypix = (int)(*nxpix * aspect_ratio);

    if (win == 0) *ierr = 2;
    else *ierr = 0;
}

mcdxcolors()
{
    int size,im_colormap,i;

    lux_get_colorinfo(&size,&im_colormap);
/*    fprintf(stderr,"Colormap size %d ,  %d\n", size, im_colormap);*/
    if (im_colormap) {
	color[0] = lux_rgb_pixel(win,1.0,1.0,1.0); 
	color[1] = lux_rgb_pixel(win,0.0,0.0,0.0); 
	if (size == 2) return(2);
	if (size == 256) {
	    float v = 256.0;
	    for(i=2;i<MAXCOLOR;i++)
	      color[i]=(long)(v=v/1.4);
	    return(MAXCOLOR);
	}
	else {
	    for(i=2;i<size;i++) 
		color[i]=(long)(i-1);
	    return(size);
	} 
    }
    else {
	if (size >= MAXCOLOR) {
	    for(i=0;i<MAXCOLOR;i++) 
	      color[i] = lux_lookup_color(win,colorname[i]);
	    return(MAXCOLOR);
	}
	else {
	    for(i=0;i<size;i++) 
	      color[i] = lux_lookup_color(win,colorname[i]);
	    return(size);
	}
    }
}

mcdxopenwin()
{
    mcdwin *w;
    static int id_number = 0;
    static int xwin = -1, ywin = -1;
    static int inc = 25;

    int wwin, hwin;

    w = (mcdwin *)malloc(sizeof(mcdwin));
    if (window == NULL) {
	window = w;
	window->next = window;
    }
    else {
	w->next = window->next;
	window->next = w;
	window = w;
    }

    if (xwin == -1) {
	xwin = width*SCREEN_XPOS_FRACTION;
	ywin = height*SCREEN_YPOS_FRACTION;
    }

    wwin = width*SCREEN_SIZE_FRACTION;
    hwin = (int)(wwin*aspect_ratio);

    xwin += inc;
    ywin += inc;
    if (xwin + wwin > width || ywin + hwin > height) {
	xwin = width*SCREEN_XPOS_FRACTION;
	ywin = height*SCREEN_YPOS_FRACTION;
    }
    
    window->id = id_number;
    id_number++;
    sprintf(window->name,"McDraw  #%d",window->id);
    window->win = lux_openwin(xwin,ywin,wwin,hwin);
    lux_set_window_name(window->win,window->name);
    lux_setup_region(window->win, 0.0,0.0,10.0,10.0);
    lux_setup_axis(window->win, 0.0,10.0,0.0,10.0*aspect_ratio); 
    window->fg_color = lux_rgb_pixel(window->win,0.0,0.0,0.0);
    window->bg_color = lux_rgb_pixel(window->win,1.0,1.0,1.0);
    window->linewidth = 0;
    window->colormap = DEFAULT_CM;

    win = window->win;
}

mcdxsetwin(id)
int *id;
{
    mcdwin *w;
    
    w = window;
    do {
	if (window->id == *id) {
	    win = window->win; 
	    lux_uniconify_window(win);
	    lux_raise_window(win);
	    return *id;
	}
	window = window->next;
    } while(window != w);
    fprintf(stderr,"Window %d not found\n", *id);
    return -1;
}

mcdxsetwin_no_raise(id)
int *id;
{
    mcdwin *w;
    
    w = window;
    do {
	if (window->id == *id) {
	    win = window->win; 
	    return *id;
	}
	window = window->next;
    } while(window != w);
    fprintf(stderr,"Window %d not found\n", *id);
    return -1;
}

mcdxcurrwin()
{
    if (window == NULL) return -1;
    return window->id;
}

mcdxnopen()
{
    mcdwin *w;
    int n_win;

    if (window == NULL) return 0;

    w = window;
    n_win = 0;

    do {
	n_win++;
	window = window->next;
    } while(window != w);

    return n_win;
}

mcdxkillwin(id)
int *id;
{
    mcdwin *w;
    int    old_id;

    old_id = window->id;
    if (mcdxsetwin_no_raise(id) < 0) return;
    lux_freewin(win);

    if (window->next == window) {
	free(window);
	window = NULL;
    }
    else {
	w = window;
        while(window->next != w) window = window->next;
	window->next = w->next;
	win = window->win;
	free(w);
	if (old_id != *id) mcdxsetwin(&old_id);
	else {
	  mcdxsetwin(id);
	  fprintf(stderr,"Current window id = %d\n",window->id);
	}
    }

}

mcdxflush()
{
  if (style == -1) return;
  if (style == 0) {
    lux_draw_pointsf(win,x,y,npoints,0);
  }
  else if (style == 1) {
    lux_draw_segmentsf(win,x,y,x,y,npoints/2,0);
  }
  npoints = 0;
}

mcdxbuffer()
{
  buffer_flag = 1;
}

mcdxnobuffer()
{
  buffer_flag = 0;
  mcdxflush();
}

mcdxlinew(iw)
unsigned int *iw;
{
  if (npoints) mcdxflush();
  lux_set_linewidth(win,*iw);
  window->linewidth = *iw;
}

mcdxcolor(ic)
int  *ic;
{
  if (npoints) mcdxflush();

  if (window->colormap == DEFAULT_CM) {
      lux_set_color(win,color[*ic%MAXCOLOR]);
      window->fg_color = color[*ic%MAXCOLOR];
  } else {
      lux_set_color(win,*ic);
      window->fg_color = *ic;
  }
}

mcdxbackg(ic)
unsigned long *ic;
{
  if (npoints) mcdxflush();
  /* Do nothing right now */
}

mcdxmove(r,s)
float *r,*s;
{
  penx = *r;
  peny = *s;
}

mcdxdraw(r,s)
float *r,*s;
{
  if (style != 1) mcdxflush();
  style = 1;
  if (buffer_flag) {
    if (npoints < MAXPOINTS) {
      x[npoints] = penx;
      y[npoints] = peny;
      npoints++;
      x[npoints] = *r;
      y[npoints] = *s;
      npoints++;
    }
    else {
      mcdxflush();
      x[npoints] = penx;
      y[npoints] = peny;
      npoints++;
      x[npoints] = *r;
      y[npoints] = *s;
      npoints++;
    }
  }
  else
    lux_draw_linef(win,penx,peny,*r,*s);
  penx = *r;
  peny = *s;
}

mcdxpoint(r,s)
float *r,*s;
{
  if (style != 0) mcdxflush();
  style = 0;
  if (buffer_flag) {
    if (npoints < MAXPOINTS) {
      x[npoints] = *r;
      y[npoints] = *s;
      npoints++;
    }
    else {
      mcdxflush();
      x[npoints] = *r;
      y[npoints] = *s;    
      npoints++;
    }
  }
  else
    lux_draw_pointf(win,*r,*s);
  penx = *r;
  peny = *s;  
}

#define REASONABLE_WAIT_TIME 50000	/* Unit = microseconds. */

mcdxgin(r,s)
float *r, *s;
{
    int i;
    char strng[1024];

    if (npoints) mcdxflush();

    set_timeout();
    i = get_wait_time();
    set_wait_time(REASONABLE_WAIT_TIME);

    get_mouse_position(win,r,s);

    set_wait_time(i);
    reset_term(strng);
}

mcdxpolyf(r,s,n,ic)
float *r,*s;
int *n;
unsigned long *ic;
{
  unsigned long fg;

  if (npoints) mcdxflush();
  fg = lux_get_fgcolor(win);
  mcdxcolor(ic);
  lux_fill_polygonf(win,r,s,*n,0);
  lux_set_color(win, fg);
  penx = r[*n-1];
  peny = s[*n-1];
}

mcdxpolyc(r,s,n)
float *r,*s;
int *n;
{
  unsigned long fg;

  if (npoints) mcdxflush();
  fg = lux_get_fgcolor(win);
  lux_set_color(lux_get_bgcolor(win));
  lux_draw_linesf(win,r,s,r,s,*n,0);
  lux_set_color(win, fg);
  penx = r[*n-1];
  peny = s[*n-1];
}

mcdxtext(h,a,string)
float *h,*a;
char *string;
{
    if (npoints) mcdxflush();
    lux_draw_string(win,penx,peny,*h,string,0);
}

mcdxclear()
{
    if (npoints) mcdxflush();
    lux_reset_window(win);
    lux_setup_region(win, 0.0,0.0,10.0,10.0);
    lux_setup_axis(win, 0.0,10.0,0.0,10.0*aspect_ratio); 
    lux_set_color(win, window->fg_color);
    lux_set_linewidth(win, window->linewidth);
}

mcdxreset()
{
    if (npoints) mcdxflush();
    lux_reset_window(win);
    lux_setup_region(win, 0.0,0.0,10.0,10.0);
    lux_setup_axis(win, 0.0,10.0,0.0,10.0*aspect_ratio); 
}

mcdxquit()
{
  lux_quick_exit();
}

mcdxidle()
{
    float r,s;
    if (npoints) mcdxflush();
    mcdxgin(&r,&s);
}

mcdxrdline(strng)
char *strng;
{

    if (npoints) mcdxflush();

    set_timeout();

    lux_getevent();

    strng[0] = '\0';
    reset_term(strng);
}

mcdxiconifywin(id)
int *id;
{
  int old_id;

  if (window == NULL) return(0);
  old_id = window->id;
  if (mcdxsetwin_no_raise(id)<0) return 0;
  
  lux_iconify_window(win);
  mcdxsetwin_no_raise(&old_id);
}

mcdxuniconifywin(id)
int *id;
{
  int old_id;

  if (window == NULL) return(0);
  old_id = window->id;
  if (mcdxsetwin_no_raise(id)<0) return 0;
  
  lux_uniconify_window(win);
  mcdxsetwin(&old_id);
}

mcdxiconifyall()
{
  mcdwin *w;

  if (window == NULL) return(0);

  w = window;
  do {
    mcdxiconifywin(&(w->id));
    w = w->next;
  }while(w != window);
}

mcdxuniconifyall()
{
  mcdwin *w;

  if (window == NULL) return(0);

  w = window;
  do {
    mcdxuniconifywin(&(w->id));
    w = w->next;
  }while(w != window);
}

mcdxflushio()
{
    clear_buffer();
    lux_clear_keyboard_buffer();
}

mcdxsetwincolormap(filename, nc, ier)
char* filename;
int* nc;
int* ier;
{
    int n;
    if ( (n = lux_set_window_colormap(win, filename)) > 0) {
	*nc = n;
	window->colormap = FILE_CM;
	*ier = 0;
    } else
	*ier = 1;
}

#include <time.h>

lux_pause(time)
int time;
{
    /* Wait for the specified number of microseconds, if possible. */

#ifdef HAS_USLEEP
    usleep(time);
#endif

#ifdef HAS_NANOSLEEP
    struct timespec t;
    time_t i;
    i = time/1000000;
    t.tv_sec = i;
    t.tv_nsec = 1000*(time - 1000000*i);
    nanosleep(&t, &t);
#endif
}

/* Alternate names for use by name-mangling compilers: */

mcdxinit_(aspect, nxpix, nypix, ncolor, ierr)
float *aspect;
int *nxpix, *nypix, *ncolor;
int *ierr;
{
  mcdxinit(aspect, nxpix, nypix, ncolor, ierr);
}

mcdxlinew_(iw)
unsigned int *iw;
{
  mcdxlinew(iw);
}

mcdxcolor_(ic)
int *ic;
{
  mcdxcolor(ic);
}

mcdxbackg_(ic)
unsigned long *ic;
{
  mcdxbackg(ic);
}

mcdxmove_(r,s)
float *r,*s;
{
  mcdxmove(r,s);
}

mcdxdraw_(r,s)
float *r,*s;
{
  mcdxdraw(r,s);
}

mcdxpoint_(r,s)
float *r,*s;
{
  mcdxpoint(r,s);
}

mcdxgin_(r,s)
float *r, *s;
{
  mcdxgin(r,s);
}

mcdxpolyf_(r,s,n,ic)
float *r,*s;
int *n;
unsigned long *ic;
{
  mcdxpolyf(r,s,n,ic);
}

mcdxpolyc_(r,s,n)
float *r,*s;
int *n;
{
  mcdxpolyc(r,s,n);
}

mcdxtext_(h,a,string)
float *h,*a;
char *string;
{
    mcdxtext(h,a,string);
}

mcdxclear_()
{
  mcdxclear();
}

mcdxreset_()
{
  mcdxreset();
}

mcdxquit_()
{
  mcdxquit();
}

mcdxidle_()
{
  mcdxidle();
}

mcdxrdline_(strng)
char *strng;
{
  mcdxrdline(strng);
}

mcdxiconifywin_(id)
int *id;
{
  mcdxiconifywin(id);
}

mcdxuniconifywin_(id)
int *id;
{
  mcdxuniconifywin(id);
}

mcdxiconifyall_()
{
  mcdxiconifyall();
}

mcdxuniconifyall_()
{
  mcdxuniconifyall();
}

mcdxflushio_()
{
  mcdxflushio();
}

mcdxsetwincolormap_(filename, nc, ier)
char* filename;
int* nc;
int* ier;
{
  mcdxsetwincolormap(filename, nc, ier);
}

lux_pause_(time)
int time;
{
  lux_pause(time);
}

mcdxsetwin_(id)
int *id;
{
  mcdxsetwin(id);
}

mcdxsetwin_no_raise_(id)
int *id;
{
  mcdxsetwin_no_raise(id);
}

mcdxcurrwin_()
{
  mcdxcurrwin();
}

mcdxnopen_()
{
  mcdxnopen();
}

mcdxkillwin_(id)
int *id;
{
  mcdxkillwin(id);
}

mcdxflush_()
{
  mcdxflush();
}

mcdxbuffer_()
{
  mcdxbuffer();
}

mcdxnobuffer_()
{
  mcdxnobuffer();
}
