
#include "xhrdplot.h"

static float linespacing;	// Set in xstarplot.  Too many references
				// to pass it around easily as an argument...
				// Keep global within this file.

			 
// Open windows for plotting and instructions, set up color schemes, etc.

void initialize_graphics(float r_factor,	     // (input) scale factor
			 bool b_flag,		     // (input) color scheme 
			 unsigned long& win,	     // ID of plot window
			 unsigned long& instr,	     // ID of status window
			 unsigned long& colwin,	     // ID of color window
			 unsigned long* c_index,     // color by index
			 unsigned long* c_energy,    // color by energy
			 int& win_size,		     // width of plot window
			 int& xorigin, int& yorigin) // origin for dialog box
{
    // This is a kludge to try to get the font size right...
    // Note that the same font is used everywhere.

    // The layout works best for 0.6 < r_factor < 1.3.

    if (r_factor >= 1.0) {
	set_default_font("9x15");
	linespacing = 1.2;
    } else if (r_factor >= 0.8) {
	set_default_font("7x13");
	if (r_factor >= 0.9)
	    linespacing = 1.15;
	else
	    linespacing = 1.1;
    } else if (r_factor >= 0.65) {
	set_default_font("6x12");
	if (r_factor >= 0.725)
	    linespacing = 1.05;
	else
	    linespacing = 1.0;
    } else if (r_factor >= 0.5) {
	set_default_font("6x10");
	if (r_factor >= 0.575)
	    linespacing = 1.05;
	else
	    linespacing = 1.0;
    } else if (r_factor >= 0.4) {
	set_default_font("6x9");
	linespacing = 1.0;
    } else {
	set_default_font("5x8");
	linespacing = 1.0;
    }

    // The twm title bar is ~27-30 pixels (*independent* of scaling),
    // but the bar doesn't count in the location of the window!
    // Mac windows need more space (frame at bottom too).
    // Assume that scale factors of 0.75 or less mean we are
    // using a Mac or PC display.  Also, for r_factor > 1.1, move
    // the  "colwin" box up and over.

    // Main window:

    int title = 27;
    int bottom = 0;
    if (r_factor <= 0.75) {
	title = 20;
	bottom = 20;
    }

    xorigin = 0;
    yorigin = 0;
    int xsize = _R_(400);
    int ysize = _R_(400);
    win = lux_openwin(xorigin, yorigin, xsize, ysize);
    lux_set_window_name(win, "Starlab");

    win_size = xsize;

    // Status window:

    if (r_factor <= 1) 
	yorigin += ysize + title + bottom;
    else
	yorigin += ysize + 2*title + bottom;

    // (Extra factor of 2 seems to improve layout...)

    ysize = _R_(260);
    if (r_factor <= 0.6) {
	xsize = (int) (xsize * 1.25);
	ysize = (int) (ysize * 1.1);
    }
    instr = lux_openwin(xorigin, yorigin, xsize, ysize);
    lux_set_window_name(instr, "Status & Commands");

    // Color window:

    if (r_factor <= 1.1 && r_factor > 0.75) {
	yorigin += ysize + title + bottom;
    } else {
	int edge = 3;
	if (r_factor <= 0.6) xsize = (int) (xsize / 1.25);
	xorigin += xsize + edge;
	yorigin = 0;
    }
    if (r_factor <= 0.6) {
	xsize = (int) (xsize * 1.25);
	ysize = (int) (ysize * 1.1);
    }
    if (r_factor <= 0.6) ysize =  (int) (ysize * 1.1);
    colwin = lux_openwin(xorigin, yorigin, xsize, ysize);
    lux_set_window_name(colwin, "Color Scheme");

    init_colors(win, c_energy, c_index, b_flag);

    lux_set_noupdate(win);
    lux_set_bgcolor(win,c_energy[background_color]);
    lux_set_window_bgcolor(win,c_energy[background_color]);
    lux_set_color(win,c_energy[default_color]);

    lux_set_noupdate(colwin);
    lux_set_bgcolor(colwin,c_energy[background_color]);
    lux_set_window_bgcolor(colwin,c_energy[background_color]);
    lux_set_color(colwin,c_energy[default_color]);

    lux_set_noupdate(instr);
    lux_set_bgcolor(instr,c_energy[background_color]);
    lux_set_window_bgcolor(instr,c_energy[background_color]);
    lux_set_color(instr,c_energy[default_color]);

    lux_setup_region(win, 1.0, 1.0, 8.0, 8.0);

    // Set up origin of dialog box.

    if (r_factor <= 1.1 && r_factor > 0.75) {
	int edge = 3;
	xorigin += xsize + edge;
	yorigin = 0;
    } else
	yorigin += ysize + title + bottom;
}

// Note overloaded function project3d:

void project3d(float x, float y, float z, float& X, float& Y,
	       float ct, float st, float cp, float sp,
	       float vx, float vy, float vz, float& dv2)
{
    float xp, yp, zp, xf, yf, zf;

    // Along z-axis rotate theta
    xp =  ct*x + st*y;
    yp = -st*x + ct*y;
    zp =  z;

    // Along y-axis rotate phi
    xf =  cp*xp + sp*zp;
    yf =  yp;
    zf = -sp*xp + cp*zp;

    // Project
    X  = yf; Y=zf;
    dv2 = (vx-xf)*(vx-xf) + (vy-yf)*(vy-yf) + (vz-zf)*(vz-zf);
}

void project3d(float x, float y, float z, float& X, float& Y,
	       float ct, float st, float cp, float sp)
{
    float xp, yp, zp, xf, yf, zf;

    // Along z-axis rotate theta
    xp =  ct*x + st*y;
    yp = -st*x + ct*y;
    zp =  z;

    // Along y-axis rotate phi
    xf =  cp*xp + sp*zp;
    yf =  yp;
    zf = -sp*xp + cp*zp;

    // Project
    X  = yf; Y=zf;
}

void draw3d_axis(unsigned long win, float lmax,
		 float ct, float st, float cp, float sp)
{
    float X[8], Y[8], Z[8], XP[8], YP[8], D[8], X0, Y0, Z0, dmax = 0.0;
    float xp,yp;
    int i, i0;
    static int od1[] = {1,2,4,3,5,6,8,7,1,2,4,3};
    static int od2[] = {2,4,3,1,6,8,7,5,5,6,8,7};

    X0 = lmax; Y0 = 0.0; Z0 = 0.0;
    for (i=0;i<8;i++) {
	X[i] = lmax*((i&0x0001)?-1:1);
	Y[i] = lmax*((i&0x0002)?-1:1);
	Z[i] = lmax*((i&0x0004)?-1:1);
	project3d(X[i], Y[i], Z[i], XP[i], YP[i],
		  ct, st, cp, sp,
		  X0, Y0, Z0, D[i]);
	if (D[i] > dmax) {dmax = D[i]; i0 = i;}
    }

    for (i=0;i<12;i++) {
	if (od1[i] == i0+1 || od2[i] == i0+1) lux_set_linestyle(win,1);
	else lux_set_linestyle(win,0);
	lux_draw_linef(win,XP[od1[i]-1],YP[od1[i]-1],XP[od2[i]-1],YP[od2[i]-1]);
    }
    lux_set_linestyle(win,0);

    // labels:
    project3d(0.0,-lmax*1.1,-lmax*1.1, xp, yp, ct,st,cp,sp);
    lux_draw_string(win, xp, yp, -0.5, "x", 0); 
    project3d(-lmax*1.1,0.0,-lmax*1.1, xp, yp, ct,st,cp,sp);
    lux_draw_string(win, xp, yp, -0.5, "y", 0); 
    project3d(-lmax*1.1,-lmax*1.1,0.0, xp, yp, ct,st,cp,sp);
    lux_draw_string(win, xp, yp, -0.5, "z", 0); 
}

void draw2d_axis(unsigned long win, float xmin, float xmax,
		 float ymin, float ymax, int k)
{
    lux_draw_axis(win);

    // X-axis label:
	    
    switch (k) {
	case 3: lux_draw_string(win, (xmin+xmax)/2, ymin, -2.2, "x", 0); break;
	case 1: lux_draw_string(win, (xmin+xmax)/2, ymin, -2.2, "y", 0); break;
	case 2: lux_draw_string(win, (xmin+xmax)/2, ymin, -2.2, "z", 0); break;
    }

    // Y-axis label:
	    
    switch (k) {
	case 3: lux_draw_vstring(win, xmin, (ymin+ymax)/2,
				 1.2, " y ", 0); break;
	case 1: lux_draw_vstring(win, xmin, (ymin+ymax)/2,
				 1.2, " z ", 0); break;
        case 2: lux_draw_vstring(win, xmin, (ymin+ymax)/2,
				 1.2, " x ", 0); break;
    }
}


void update_with_delay(unsigned long win, float t)
{
    lux_update_fg(win);
    if (t > 0.0)  lux_pause((int)(t*1000000.0));
}

// Overloaded function:

void show_instructions(unsigned long  win, float r, char* buffer,
		       int update)
{
    float statusx, statusy; 
    char s[255];
    int ptr = 0, len, line = 0;

    lux_clear_window(win);
    lux_reconvert_rcoord(win,0,0,&statusx, &statusy);
    istrstream ins(buffer,strlen(buffer));
    while(ins.get(s,255,'\n')) {
	lux_draw_image_string(win, statusx, statusy,
			      -(line+1)*linespacing, s, -1);
	ins.get(s[0]);  /* Clean "\n" */
	line = line+1;
    }
    if (update) lux_update_fg(win);
}


void show_instructions(unsigned long win, float r, char* buffer,
		       int line, int update)
{
    float statusx, statusy; 
    char s[255];
    int ptr = 0, len;

    lux_reconvert_rcoord(win,0,0,&statusx, &statusy);
    istrstream ins(buffer,strlen(buffer));
    while(ins.get(s,255,'\n')) {
	lux_draw_string(win, statusx, statusy, -(line+1)*linespacing, s, -1);
	ins.get(s[0]);  /* Clean "\n" */
	line = line+1;
    }
    if (update) lux_update_fg(win);
}


void show_main_instructions(unsigned long instr,float r,int d,int u)
{
    if (d)
	show_instructions(instr, r,
"Running in 3-D mode\n\
  d: dialog, R: rotate, i: idle, q: quit\n\
  p/P         decrease/increase point size\n\
  t           add/remove tracking\n\
  n/l         show nodes and/or links\n\
  2           enter 2-D mode\n\
  x/y/k       2-D and change viewing axis\n\
  e           toggle energy color\n\
  < > ^ V     rotate left right up down\n\
  +/-         change delay time (0: fastest)\n\
  z/Z         zoom in or out \n\
  arrow/page  shift in x, y, z\n\
  h/home      restore origin to (0,0,0)\n\
  a           resize to enclose all stars\n\
",u);

    else

	show_instructions(instr, r,
"Running in 2-D mode\n\
  d: dialog, i: idle, q: quit\n\
  p/P     decrease/increase point size\n\
  t       add/remove tracking\n\
  n/l     show nodes and/or links (r: root)\n\
  3       enter 3-D mode\n\
  x/y/k   change viewing axis\n\
  o/O     set origin using mouse\n\
  e       toggle energy color\n\
  +/-     change delay time (0: fastest)\n\
  z/Z     zoom in or out \n\
  arrow   shift horizontally/vertically\n\
  h/home  restore origin to (0,0,0)\n\
  a       resize to enclose all stars\n\
",u);
}

void format_and_show_instructions(unsigned long co, float r,
					unsigned long *c_i, int index, int tab,
					char* cstring, int line, int u)
{
    if (index > 0) lux_set_color(co, c_i[index]);

    char temp[80], temp1[80];
    sprintf(temp, " ");
    for (int i = 0; i < tab; i++) strcat(temp, TAB);
    sprintf(temp1, "%d - %s", index, cstring);
    strcat(temp, temp1);
    
    show_instructions(co, r, temp, line, u);
}

void show_color_scheme(unsigned long co, unsigned long *c_e, unsigned long *c_i,
		       float r, char ce, bool b_flag, int u)
{
    lux_clear_window(co);

    if (ce) {
	show_instructions(co, r, "Color by energy:\n",0,u);
	lux_set_color(co, c_e[default_color]);
	if (b_flag) 
	    show_instructions(co, r,
			      "  White  = default (bound to cluster)\n",1,u);
	else 
	    show_instructions(co, r,
			      "  Black  = default (bound to cluster)\n",1,u);
	lux_set_color(co, c_e[bound_single]);
	show_instructions(co, r, "  Blue   = bound to nearest neighbor\n",2,u);
	lux_set_color(co, c_e[bound_binary]);
	show_instructions(co, r, "  Green  = bound binary\n",3,u);
	lux_set_color(co, c_e[unbound_single]);
	show_instructions(co, r, "  Red    = unbound single star\n",4,u);
	lux_set_color(co, c_e[unbound_binary]);
	show_instructions(co, r, "  Gold   = unbound binary\n",5,u);
	lux_set_color(co, c_e[default_color]);

    } else {

	show_instructions(co, r, "Color by index:\n",0,u);
	lux_set_color(co, c_i[1]);

	if (b_flag) {
	    format_and_show_instructions(co, r, c_i,  1, 0, "white",   1, u);
	    format_and_show_instructions(co, r, c_i,  2, 1, RV_COLOR2, 1, u);
	    format_and_show_instructions(co, r, c_i,  3, 2, RV_COLOR3, 1, u);
	    format_and_show_instructions(co, r, c_i,  4, 0, RV_COLOR4, 2, u);
	    format_and_show_instructions(co, r, c_i,  5, 1, RV_COLOR5, 2, u);
	    format_and_show_instructions(co, r, c_i,  6, 2, RV_COLOR6, 2, u);
	    format_and_show_instructions(co, r, c_i,  7, 0, RV_COLOR7, 3, u);
	    format_and_show_instructions(co, r, c_i,  8, 1, RV_COLOR8, 3, u);
	    format_and_show_instructions(co, r, c_i,  9, 2, RV_COLOR9, 3, u);
	    format_and_show_instructions(co, r, c_i, 10, 0, RV_COLORa, 4, u);
	    format_and_show_instructions(co, r, c_i, 11, 1, RV_COLORb, 4, u);
	    format_and_show_instructions(co, r, c_i, 12, 2, RV_COLORc, 4, u);
	    format_and_show_instructions(co, r, c_i, 13, 0, RV_COLORd, 5, u);
	    format_and_show_instructions(co, r, c_i, 14, 1, RV_COLORe, 5, u);
	    format_and_show_instructions(co, r, c_i, 15, 2, RV_COLORf, 5, u);
	    format_and_show_instructions(co, r, c_i, 16, 0, RV_COLORg, 6, u);
	} else {
	    format_and_show_instructions(co, r, c_i,  1, 0, "black",   1, u);
	    format_and_show_instructions(co, r, c_i,  2, 1, NV_COLOR2, 1, u);
	    format_and_show_instructions(co, r, c_i,  3, 2, NV_COLOR3, 1, u);
	    format_and_show_instructions(co, r, c_i,  4, 0, NV_COLOR4, 2, u);
	    format_and_show_instructions(co, r, c_i,  5, 1, NV_COLOR5, 2, u);
	    format_and_show_instructions(co, r, c_i,  6, 2, NV_COLOR6, 2, u);
	    format_and_show_instructions(co, r, c_i,  7, 0, NV_COLOR7, 3, u);
	    format_and_show_instructions(co, r, c_i,  8, 1, NV_COLOR8, 3, u);
	    format_and_show_instructions(co, r, c_i,  9, 2, NV_COLOR9, 3, u);
	    format_and_show_instructions(co, r, c_i, 10, 0, NV_COLORa, 4, u);
	    format_and_show_instructions(co, r, c_i, 11, 1, NV_COLORb, 4, u);
	    format_and_show_instructions(co, r, c_i, 12, 2, NV_COLORc, 4, u);
	    format_and_show_instructions(co, r, c_i, 13, 0, NV_COLORd, 5, u);
	    format_and_show_instructions(co, r, c_i, 14, 1, NV_COLORe, 5, u);
	    format_and_show_instructions(co, r, c_i, 15, 2, NV_COLORf, 5, u);
	    format_and_show_instructions(co, r, c_i, 16, 0, NV_COLORg, 6, u);
	}
	lux_set_color(co, c_i[default_color]);
    }
}

void init_colors(unsigned long win, unsigned long *ce, unsigned long *ci,
		 bool b_flag)
{
    // Note earlier definitions:
    //		background_color = 0
    //		default_color    = 1

    if (b_flag) {

	// Reverse video color scheme			(background = black)

	ce[background_color] = ci[background_color]
	    = lux_rgb_pixel(win, 0.0,0.0,0.0);        		// black
	ce[default_color] = ci[default_color]
	    = lux_rgb_pixel(win, 1.0,1.0,1.0);           	// white
    
	// Color by energy:

	ce[bound_single]   = lux_lookup_color(win, "deepskyblue");
	ce[bound_binary]   = lux_lookup_color(win, "green");
	ce[unbound_single] = lux_lookup_color(win, "red");
	ce[unbound_binary] = lux_lookup_color(win, "goldenrod");

	// Color by index:

	ci[2]  = lux_lookup_color(win, RV_COLOR2);
	ci[3]  = lux_lookup_color(win, RV_COLOR3);
	ci[4]  = lux_lookup_color(win, RV_COLOR4);
	ci[5]  = lux_lookup_color(win, RV_COLOR5);
	ci[6]  = lux_lookup_color(win, RV_COLOR6);
	ci[7]  = lux_lookup_color(win, RV_COLOR7);
	ci[8]  = lux_lookup_color(win, RV_COLOR8);
	ci[9]  = lux_lookup_color(win, RV_COLOR9);
	ci[10] = lux_lookup_color(win, RV_COLORa);
	ci[11] = lux_lookup_color(win, RV_COLORb);
	ci[12] = lux_lookup_color(win, RV_COLORc);
	ci[13] = lux_lookup_color(win, RV_COLORd);
	ci[14] = lux_lookup_color(win, RV_COLORe);
	ci[15] = lux_lookup_color(win, RV_COLORf);
	ci[16] = lux_lookup_color(win, RV_COLORg);
    
    } else {					     // (background = white)

	ce[background_color] = ci[background_color]
	    = lux_rgb_pixel(win, 1.0,1.0,1.0);        		// black
	ce[default_color] = ci[default_color]
	    = lux_rgb_pixel(win, 0.0,0.0,0.0);           	// white
    
	// Color by energy:

	ce[bound_single]   = lux_lookup_color(win, "blue");
	ce[bound_binary]   = lux_lookup_color(win, "limegreen");
	ce[unbound_single] = lux_lookup_color(win, "red");
	ce[unbound_binary] = lux_lookup_color(win, "goldenrod");

	// Color by index:
    
	ci[2]  = lux_lookup_color(win, NV_COLOR2);
	ci[3]  = lux_lookup_color(win, NV_COLOR3);
	ci[4]  = lux_lookup_color(win, NV_COLOR4);
	ci[5]  = lux_lookup_color(win, NV_COLOR5);
	ci[6]  = lux_lookup_color(win, NV_COLOR6);
	ci[7]  = lux_lookup_color(win, NV_COLOR7);
	ci[8]  = lux_lookup_color(win, NV_COLOR8);
	ci[9]  = lux_lookup_color(win, NV_COLOR9);
	ci[10] = lux_lookup_color(win, NV_COLORa);
	ci[11] = lux_lookup_color(win, NV_COLORb);
	ci[12] = lux_lookup_color(win, NV_COLORc);
	ci[13] = lux_lookup_color(win, NV_COLORd);
	ci[14] = lux_lookup_color(win, NV_COLORe);
	ci[15] = lux_lookup_color(win, NV_COLORf);
	ci[16] = lux_lookup_color(win, NV_COLORg);
    
    }

    // Problems --> foreground color.

    for (int ke = bound_single; ke <= unbound_binary; ke++)
	if (ce[ke] == ce[background_color]) ce[ke] = ce[default_color];
    for (int ki = 2; ki < 16; ki++)
	if (ci[ki] == ci[background_color]) ci[ki] = ci[default_color];
}

void set_limits(float* origin, float lmax3d,
		int kx, float& xmin, float& xmax,
		int ky, float& ymin, float& ymax)
{
    xmin = origin[kx] - lmax3d;
    xmax = origin[kx] + lmax3d;
    ymin = origin[ky] - lmax3d;
    ymax = origin[ky] + lmax3d;
}

float interp_to_x(float r, float s, float rr, float ss, float x)
{
    return s + (ss - s) * (x - r) / (rr - r);
}


float interp_to_y(float r, float s, float rr, float ss, float y)
{
    return r + (rr - r) * (y - s) / (ss - s);
}
