#include "win.h"

#define FABS(x)    ((x)>=0.0?(x):-((x)))
#define FINT(x)    ((x)>=0.0?(long)((x)):(long)((x)-0.999999999999999999))
#define FMOD(x,y)  ((x) - FINT((x)/(y))*(y))

#define MAXSIZE  1000     /*Maximum screen size in pixel */

extern lux_wins *get_currentwin();

lux_set_step(min, max, ln, step, start, tick)  /* Assume min < max always */
double min, max, *step, *start;
long  *tick, ln;
{
    double xlog,x,tt,eps;
    long   ilog;

    if (ln) {
	xlog = log10(max) - log10(min);
	eps  = xlog/10.0/MAXSIZE;
	ilog = (long) (xlog + 0.5);
	if (ilog>=4) {
	    *step  = pow(10.0, (double)((long)((ilog+2)/4)));
	    *tick  = (ilog+2)/4*10;
	}
	else {
	    *step  = 10.0;
	    *tick  = 10;	  
	}
	xlog = log10(min);
	*start = pow(10.0,(double)(FINT(xlog)));
/*	fprintf(stderr,"min=%lg  start=%lg step=%lg\n",min,*start,*step);*/
	while(FABS(log10(*start/min)) > eps && *start < min) *start *= 10.0;
    }
    else {
	eps = (max - min)/10.0/MAXSIZE;
	xlog = log10(max - min);
	ilog = FINT(xlog);
	x = pow(10.0, (double)ilog);
	tt = (max - min)/x;
	if (tt - (long)tt > 0.49999) tt = tt + 1.0;
	switch((long)tt) {
	    case 1: 
	    case 5: 
	    case 10:
	      *step = ((long)tt)/5.0*x; 
	      *tick = 10;
	      break;
	    case 2: *step = ((long)tt)/4.0*x; *tick = 10;  break;
	    case 4: 
	    case 8: 
	      *step = ((long)tt)/4.0*x; 
	      *tick = 10; 
	      break;
	    case 3: 
	    case 7: 
	      *step = ((long)tt+1)/4.0*x; 
	      *tick = 10;
	      break;
	    case 6: 
	      *step = ((long)tt-1)/5.0*x; 
	      *tick = 10;
	      break;
	    case 9: 
	      *step = ((long)tt+1)/5.0*x; 
	      *tick = 10;
	      break;
	    default:
	      fprintf(stderr,"Calculation error!\nCheck your axis limits!!\n");
	}


	if (FABS(min) < *step) {
	    if (FABS(*step - 2.0*FABS(min)) < eps) {
		*start = min;
	    }
	    else
	      *start = (min>=0.0)?0.0:-(*step);
	}
	else {
	    xlog = log10(FABS(min));
	    ilog = FINT(xlog);
	    x = pow(10.0, (double)ilog);
	    *start = FINT(min/x)*x;
	    if (FMOD(min-*start,*step/2.0) < eps || 
		FMOD(min-*start,*step/2.0) > *step/2.0 - eps)
	      *start = min - *step;
	}
/*	fprintf(stderr,"min=%lg  start=%lg step=%lg mod=%lg eps=%lg\n",min,*start,*step, FMOD(min-*start,*step/2.0),eps); */
	while(FABS(min-*start) > eps && min > *start) *start = *start + *step;
    }
}


lux_draw_axis(win)
Window win;
{
    register lux_wins *current;
    char  msg[20];
    double fx,fy,rfx,rfy,epsx,epsy;
    double ystep,xstep,ys,xs,x,y;
    double maxt,mayt,mixt,miyt,sxt,syt,rmaxt,rmayt,rmixt,rmiyt,rsxt,rsyt;
    long   i,xt,yt;

    current = get_currentwin(win);

    if (current->win.xmin >=  current->win.xmax ||
	current->win.ymin >=  current->win.ymax)  {
	fprintf(stderr,
		"Inproper limits setting in window %ld, can't draw axis!\n",
		win);
	fprintf(stderr, "xmin=%lg, xmax=%lg;  ymin=%lg, ymax=%lg\n",
		current->win.xmin,current->win.xmax,
		current->win.ymin,current->win.ymax);
	return(0);
    }
      

    lux_draw_linef(win, current->win.xmin, current->win.ymin, 
		        current->win.xmax, current->win.ymin);
    lux_draw_linef(win, current->win.xmax, current->win.ymin,
		        current->win.xmax, current->win.ymax);
    lux_draw_linef(win, current->win.xmax, current->win.ymax,
		        current->win.xmin, current->win.ymax);
    lux_draw_linef(win, current->win.xmin, current->win.ymax,
		        current->win.xmin, current->win.ymin);

/* 
   Setup tick-mark ponits 
   (maxt, mayt) (rmaxt, rmayt) Major tick marks(leftside, rightside) 
   (mixt, miyt) (rmixt, rmiyt) Minor tick marks(leftside, rightside) 
   (sxt,  syt)  (rsxt,  rsyt)  Small tick marks(leftside, rightside) 
*/

    if (current->win.lnx) epsx = (log10(current->win.xmax) - 
				 log10(current->win.xmin))/10.0/MAXSIZE;
    else epsx = (current->win.xmax - current->win.xmin)/10.0/MAXSIZE;
    if (current->win.lny) epsy = (log10(current->win.ymax) - 
				 log10(current->win.ymin))/10.0/MAXSIZE;
    else epsy = (current->win.ymax - current->win.ymin)/10.0/MAXSIZE;

    if (current->win.lnx) {
	double dx, minlog, maxlog;

	minlog = log10(current->win.xmin);
	maxlog = log10(current->win.xmax);
	dx = (maxlog - minlog)/60.0;
	
	maxt  = pow(10.0, dx + minlog);
	rmaxt = pow(10.0, maxlog - dx);

	mixt  = pow(10.0, dx*3.0/4.0 + minlog);
	rmixt = pow(10.0, maxlog - dx*3.0/4.0);

	sxt   = pow(10.0, dx/2.0 + minlog);
	rsxt  = pow(10.0, maxlog - dx/2.0);
    }
    else {
	double dx;

	dx = (current->win.xmax - current->win.xmin)/60.0;

	maxt  = dx + current->win.xmin;
	rmaxt = current->win.xmax - dx;

	mixt  = dx*3.0/4.0 + current->win.xmin;
	rmixt = current->win.xmax - dx*3.0/4.0;
	
	sxt   = dx/2.0 + current->win.xmin;
	rsxt  = current->win.xmax - dx/2.0;
    }

    if (current->win.lny) {
	double dy, minlog, maxlog;

	minlog = log10(current->win.ymin);
	maxlog = log10(current->win.ymax);
	dy = (maxlog - minlog)/60.0;

	mayt  = pow(10.0, dy + minlog);
	rmayt = pow(10.0, maxlog - dy);

	miyt  = pow(10.0, dy*3.0/4.0 + minlog);
	rmiyt = pow(10.0, maxlog - dy*3.0/4.0);

	syt   = pow(10.0, dy/2.0 + minlog);
	rsyt  = pow(10.0, maxlog - dy/2.0);
	
    }
    else {
	double dy;
	dy = (current->win.ymax - current->win.ymin)/60.0;

	mayt  = dy + current->win.ymin;
	rmayt = current->win.ymax - dy;

	miyt  = dy*3.0/4.0 + current->win.ymin;
	rmiyt = current->win.ymax - dy*3.0/4.0;
	
	syt   = dy/2.0 + current->win.ymin;
	rsyt  = current->win.ymax - dy/2.0;
    }


     fy = current->win.ymin;
    rfy = current->win.ymax;

    lux_set_step(current->win.xmin, current->win.xmax, current->win.lnx, 
		 &xstep, &fx, &xt);

    if (current->win.lnx) {
	while(log10(fx) <= log10(current->win.xmax) + epsx) {
	    for(i=0;i<xt;i++) {
		x = (i%10+1)*pow(10.0,(double)(i/10))*fx/xstep;
		if (x < current->win.xmin) continue;
		if (i%10 == 0) {
		    lux_draw_linef(current->win.window, x,  fy, x,  mayt);
		    lux_draw_linef(current->win.window, x, rfy, x, rmayt);
		}
		else {
		    lux_draw_linef(current->win.window, x,  fy, x,  syt);
		    lux_draw_linef(current->win.window, x, rfy, x, rsyt);
		}
	    }
	    sprintf(msg,"%lg",fx);
	    lux_draw_string(current->win.window,fx,fy,-1.1,msg,0);
	    fx = fx * xstep;
	}

	for(i=0;i<xt;i++) {
	    x = (i%10+1)*pow(10.0,(double)(i/10))*fx/xstep;
	    if (x > current->win.xmax) break;
	    if (i%10 == 0) {
		lux_draw_linef(current->win.window, x,  fy, x,  mayt);
		lux_draw_linef(current->win.window, x, rfy, x, rmayt);
	    }
	    else {
		lux_draw_linef(current->win.window, x,  fy, x,  syt);
		lux_draw_linef(current->win.window, x, rfy, x, rsyt);
	    }
	}
    }
    else {

	float axmax = (current->win.xmax < 0 ? -current->win.xmax
		       			     : current->win.xmax);
	float axmin = (current->win.xmin < 0 ? -current->win.xmin
		       			     : current->win.xmin);
	if (axmin > axmax) axmax = axmin;

	if (current->win.xmin * current->win.xmax <= 0) {

	    /* If the limits change sign, try to ensure that 0 is plotted. */

	    float newfx = 0;
	    while (newfx - xstep >= current->win.xmin) newfx -= xstep;

	    /* fprintf(stderr, "changing fx from %f to %f\n", fx, newfx); */
	    fx = newfx;
	}

	xs = xstep/xt;
	while (fx <= current->win.xmax + epsx) {

	    /* Deal with rounding error... */

	    float afx = (fx < 0 ? -fx : fx);
	    if (afx < 1.e-6*axmax) fx = 0;

	    for(i=0;i<xt;i++) {
		x = (i-xt)*xs + fx;
		if (x < current->win.xmin) continue;
		if (i%10 == 0) { 
		    lux_draw_linef(current->win.window, x,  fy, x,  mayt);
		    lux_draw_linef(current->win.window, x, rfy, x, rmayt);
		}
		else if (i%5 == 0) {
		    lux_draw_linef(current->win.window, x,  fy, x,  miyt);
		    lux_draw_linef(current->win.window, x, rfy, x, rmiyt);
		}
		else {
		    lux_draw_linef(current->win.window, x,  fy, x,  syt);
		    lux_draw_linef(current->win.window, x, rfy, x, rsyt);
		}
	    }
	    sprintf(msg,"%.4lg",fx);
	    lux_draw_string(current->win.window,fx,fy,-1.1,msg,0);
	    fx = fx + xstep;
	}

	for(i=0;i<xt;i++) {
	    x = (i-xt)*xs + fx;
	    if (x > current->win.xmax) break;
	    if (i%10 == 0) {
		lux_draw_linef(current->win.window, x,  fy, x,  mayt);
		lux_draw_linef(current->win.window, x, rfy, x, rmayt);
	    }
	    else if (i%5 == 0) {
		lux_draw_linef(current->win.window, x,  fy, x,  miyt);
		lux_draw_linef(current->win.window, x, rfy, x, rmiyt);
	    }
	    else {
		lux_draw_linef(current->win.window, x,  fy, x,  syt);
		lux_draw_linef(current->win.window, x, rfy, x, rsyt);
	    }
	}
    }


     fx = current->win.xmin;
    rfx = current->win.xmax;

    lux_set_step(current->win.ymin, current->win.ymax, current->win.lny, 
		 &ystep, &fy, &yt);

    if (current->win.lny) {
	while(log10(fy) <= log10(current->win.ymax) + epsy) {
	    for(i=0;i<yt;i++) {
		y = (i%10+1)*pow(10.0,(double)(i/10))*fy/ystep;
		if (y < current->win.ymin) continue;
		if (i%10 == 0) {
		    lux_draw_linef(current->win.window,  fx, y,  maxt, y);
		    lux_draw_linef(current->win.window, rfx, y, rmaxt, y);
		}
		else {
		    lux_draw_linef(current->win.window,  fx, y,  sxt,  y);
		    lux_draw_linef(current->win.window, rfx, y, rsxt,  y);
		}
	    }
	    sprintf(msg,"%lg",fy);
	    lux_draw_vstring(current->win.window,fx,fy,0.1,msg,0);
	    fy = fy * ystep;
	}
	for(i=0;i<yt;i++) {
	    y = (i%10+1)*pow(10.0,(double)(i/10))*fy/ystep;
	    if (y > current->win.ymax) break;
	    if (i%10 == 0) {
		lux_draw_linef(current->win.window,  fx, y,  maxt, y);
		lux_draw_linef(current->win.window, rfx, y, rmaxt, y);
	    }
	    else {
		lux_draw_linef(current->win.window,  fx, y,  sxt,  y);
		lux_draw_linef(current->win.window, rfx, y, rsxt,  y);
	    }
	}
    }
    else {

	float aymax = (current->win.ymax < 0 ? -current->win.ymax
		       			     : current->win.ymax);
	float aymin = (current->win.ymin < 0 ? -current->win.ymin
		       			     : current->win.ymin);
	if (aymin > aymax) aymax = aymin;

	if (current->win.ymin * current->win.ymax <= 0) {

	    /* If the limits change sign, try to ensure that 0 is plotted. */

	    float newfy = 0;
	    while (newfy - ystep >= current->win.ymin) newfy -= ystep;

	    /* fprintf(stderr, "changing fy from %f to %f\n", fy, newfy); */
	    fy = newfy;
	}

	ys = ystep/yt;
	while(fy <= current->win.ymax + epsy) {
	    for(i=0;i<yt;i++) {

		/* Deal with rounding error... */

		float afy = (fy < 0 ? -fy : fy);
		if (afy < 1.e-6*aymax) fy = 0;

		y = (i-yt)*ys + fy;
		if (y < current->win.ymin) continue;
		if (i%10 == 0) {
		    lux_draw_linef(current->win.window,  fx, y,  maxt, y);
		    lux_draw_linef(current->win.window, rfx, y, rmaxt, y);
		}
		else if (i%5 == 0) {
		    lux_draw_linef(current->win.window,  fx, y,  mixt, y);
		    lux_draw_linef(current->win.window, rfx, y, rmixt, y);
		}
		else {
		    lux_draw_linef(current->win.window,  fx, y,  sxt,  y);
		    lux_draw_linef(current->win.window, rfx, y, rsxt,  y);
		}
	    }
	    sprintf(msg,"%.4lg",fy);
	    lux_draw_vstring(current->win.window,fx,fy,0.1,msg,0);
	    fy = fy + ystep;
	}
	for(i=0;i<yt;i++) {
	    y = (i-yt)*ys + fy;
	    if (y > current->win.ymax) break;
	    if (i%10 == 0) {
		lux_draw_linef(current->win.window,  fx, y,  maxt, y);
		lux_draw_linef(current->win.window, rfx, y, rmaxt, y);
	    }
	    else if (i%5 == 0) {
		lux_draw_linef(current->win.window,  fx, y,  mixt, y);
		lux_draw_linef(current->win.window, rfx, y, rmixt, y);
	    }
	    else {
		lux_draw_linef(current->win.window,  fx, y,  sxt,  y);
		lux_draw_linef(current->win.window, rfx, y, rsxt,  y);
	    }
	}
    }
    return(1);
}

