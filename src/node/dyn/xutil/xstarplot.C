
//// xstarplot:  plot an N-body system in an X environment.
////
//// Options:    -a    specify viewing axis [3 = z]
////             -d    dimensionality of plot [2]
////             -D    delay time between plots, in ms [0]
////             -e    color by energy [no]
////             -f    solid-color stars [true]
////             -l    specify limits [take from initial data]
////             -p    pipe cin to cout [no]
////             -P    specify point size (relative to axis/30) [1]
////             -r    "reverse video" (black background) [true]
////             -s    window size, relative to "standard" 400 pixels [1]
////             -t    show tree structure [no]

//.............................................................................
//    version 1:  May 1989   Piet Hut            email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut      adapted to the new C++-based Starlab
//    version 3:  Dec 1992   Steve McMillan	 email: steve@zonker.drexel.edu
//			     Drexel University, Philadelphia, PA, USA
//			     Converted original ASCIIplot to X --- FIRST DRAFT!
//    version 4:  Jan 1994   Biao Lu             email: biao@eagle.drexel.edu
//                           Drexel University, Philadelphia, PA, USA
//                           Use new X interface for movie.
//	      Aug/Sep 1994   Cleaned up and expanded, SMcM
//.............................................................................
//  non-local functions: 
//    xstarplot
//.............................................................................
//  uses:
//    liblux.a (X interface)
//    The X graphics library
//.............................................................................

// NOTE: On some systems, it appears that overflow can occur, causing the
//	 win.lnx and win.lny entries to be overwritten and leading to
//	 spurious "log10" errors in draw.c (in lux_set_item, specifically).
//	 This isn't fatal, but it should be fixed someday...

#include "dyn.h"
#include <sstream>

#define HALF_ROOT2   0.707106781186547524400844362105	// Root2/2
#define ROOT2        1.414213562373095048801688724210	// sqrt(2)

#define GRAPH_WINDOW    1
#define POPUP_WINDOW    2
#define DIALOG_WINDOW   3
#define BUTTON_WINDOW   4
#define INPUT_WINDOW    5
#define TEXT_WINDOW     6    

#define NO_TYPE         0
#define OK_BUTTON       1
#define CANCEL_BUTTON   2
#define CHECK_BUTTON    3
#define OK_KEEP_BUTTON  4

// The following #defines used to be conditional on g++...
// The LUX functions are in lux/interface.c.

#define lux_openwin                 LUX_openwin
#define lux_setup_region            LUX_setup_region
#define lux_clear_current_region    LUX_clear_current_region
#define lux_setup_axis              LUX_setup_axis
#define lux_draw_linef              LUX_draw_linef
#define lux_draw_pointf             LUX_draw_pointf
#define lux_draw_rectanglef         LUX_draw_rectanglef
#define lux_draw_arcf               LUX_draw_arcf
#define lux_fill_arcf               LUX_fill_arcf
#define lux_draw_axis               LUX_draw_axis
#define lux_getevent                LUX_getevent 
#define lux_exit                    LUX_quick_exit	// Note!
#define lux_set_color               LUX_set_color
#define lux_lookup_color            LUX_lookup_color
#define lux_rgb_pixel               LUX_rgb_pixel
#define lux_draw_vstring            LUX_draw_vstring 
#define lux_draw_string             LUX_draw_string
#define lux_draw_image_string       LUX_draw_image_string
#define lux_check_keypress          LUX_check_keypress
#define lux_check_buttonpress       LUX_check_buttonpress
#define lux_open_dialog             LUX_open_dialog
#define lux_draw_palette            LUX_draw_palette 
#define lux_set_item                LUX_set_item
#define lux_get_itemvalue           LUX_get_itemvalue
#define lux_update_itemvalue        LUX_update_itemvalue
#define lux_clear_window            LUX_clear_window
#define lux_reset_window            LUX_reset_window
#define lux_update_fg               LUX_update_fg
#define lux_show_dialog             LUX_show_dialog
#define lux_reconvert_rcoord        LUX_reconvert_rcoord
#define lux_set_linestyle           LUX_set_linestyle
#define lux_set_window_name         LUX_set_window_name      
#define lux_set_window_bgcolor      LUX_set_window_bgcolor      
#define lux_set_bgcolor             LUX_set_bgcolor      
#define lux_set_noupdate            LUX_set_noupdate
#define lux_next_keypress           LUX_next_keypress

extern "C" unsigned long lux_openwin(int , int , int , int );
extern "C" int lux_set_window_name(unsigned long, char*);
extern "C" int lux_setup_region(unsigned long, float,float,float,float );
extern "C" int lux_clear_current_region(unsigned long);
extern "C" int lux_setup_axis(unsigned long, float,float,float,float );
extern "C" int lux_draw_linef(unsigned long,float,float,float,float);
extern "C" int lux_draw_pointf(unsigned long,float,float);
extern "C" int lux_draw_rectanglef(unsigned long,float,float,float,float);
extern "C" int lux_draw_arcf(unsigned long,float,float,float,float,float,float);
extern "C" int lux_fill_arcf(unsigned long,float,float,float,float,float,float);
extern "C" int lux_draw_axis(unsigned long);
extern "C" int lux_getevent();
extern "C" int lux_exit();
extern "C" int lux_set_color(unsigned long, long);
extern "C" int lux_set_window_bgcolor(unsigned long, long);
extern "C" int lux_set_bgcolor(unsigned long, long);
extern "C" unsigned long lux_rgb_pixel(unsigned long, float, float,float);
extern "C" unsigned long lux_lookup_color(unsigned long, char*);
extern "C" int lux_draw_string(unsigned long, float, float, float,
			       char*, char);
extern "C" int lux_draw_vstring(unsigned long, float, float, float,
				char*, char);
extern "C" int lux_draw_image_string(unsigned long, float, float, float,
				     char*, char);
extern "C" int lux_check_keypress(unsigned long,char);
extern "C" int lux_check_buttonpress(unsigned long);
extern "C" unsigned long lux_open_dialog(int, int, int, int);
extern "C" int lux_set_item(unsigned long, int, int, int, int,
			    int, int, char*);
extern "C" int lux_draw_palette(unsigned long);
extern "C" int lux_get_itemvalue(unsigned long, int, int, int, char*);
extern "C" int lux_update_itemvalue(unsigned long, int, int, int, char*);
extern "C" int lux_clear_window(unsigned long);
extern "C" int lux_reset_window(unsigned long);
extern "C" int lux_update_fg(unsigned long);
extern "C" int lux_show_dialog(unsigned long);
extern "C" int lux_reconvert_rcoord(unsigned long, int, int, float*, float*);
extern "C" int lux_set_linestyle(unsigned long, int);
extern "C" int lux_set_noupdate(unsigned long);
extern "C" int lux_next_keypress(unsigned long, char*, char*, char*, char*);

// These should probably be made into "interface" routines someday...

extern "C" int get_mouse_position(unsigned long, float*, float*);
extern "C" void set_default_font(char*);
extern "C" void lux_pause(int);

#define background_color 0
#define default_color    1
#define bound_single     2
#define bound_binary     3
#define unbound_single   4
#define unbound_binary   5

// Define index color scheme here (NV = normal video, RV = reverse video):

#define N_COLORS 16

// First 8 are good, next 8 need work, especially in the NV case...

#define NV_COLOR1 "black"
#define NV_COLOR2 "red"
#define NV_COLOR3 "limegreen"
#define NV_COLOR4 "blue"
#define NV_COLOR5 "gold"
#define NV_COLOR6 "magenta"
#define NV_COLOR7 "dark goldenrod"
#define NV_COLOR8 "lightpink"
#define NV_COLOR9 "aquamarine"
#define NV_COLORa "cyan"
#define NV_COLORb "lightgrey"
#define NV_COLORc "turquoise"
#define NV_COLORd "gold"
#define NV_COLORe "thistle"
#define NV_COLORf "beige"
#define NV_COLORg "plum"

#define RV_COLOR1 "white"
#define RV_COLOR2 "red"
#define RV_COLOR3 "green"
#define RV_COLOR4 "lightblue"
#define RV_COLOR5 "yellow"
#define RV_COLOR6 "magenta"
#define RV_COLOR7 "orange"
#define RV_COLOR8 "pink"
#define RV_COLOR9 "aquamarine"
#define RV_COLORa "cyan"
#define RV_COLORb "lightgrey"
#define RV_COLORc "turquoise"
#define RV_COLORd "gold"
#define RV_COLORe "thistle"
#define RV_COLORf "beige"
#define RV_COLORg "plum"

#define TAB "               "

//#define FAC3D		2.13
//#define FAC3D		1.5
#define FAC3D		1.75

#define ZOOM		ROOT2
#define PFAC		1.1892

enum   {colorenergy=1, tracking, graph3dim};
enum   {xminimum=1, xmaximum, yminimum,	ymaximum,
	pointsize, lmax3D,
	theta3D, phi3D,  
	DelayTime, dtheta3D, color,
	Origin, Xorigin, Yorigin, Zorigin,
	View2D, view2D, originstar, OriginStar};
enum   {ok=1, ok_keep, cancel};

char   temp_buffer[255];	// Convenient to share this temporary space.

void accumulate_potential_energy(dyn* bj, dyn*bi,
				 real& epot, real& rmin, dynptr& bmin)

// Determine the potential energy of bi with respect to bj, along
// with bi's nearest neighbor and minimum neighbor distance.

{
    if (bj->get_oldest_daughter() != (dyn*)NULL)
	for (dyn * bb = bj->get_oldest_daughter(); bb != (dyn*)NULL;
	    bb = bb->get_younger_sister()) {
	    accumulate_potential_energy(bb, bi, epot, rmin, bmin);
	}
    else
	if (bi != bj) {
	    vec d_pos = bi->get_pos() - bj->get_pos();
	    real mi = bi->get_mass();
	    real mj = bj->get_mass();
	    real r = sqrt(d_pos * d_pos);
	    epot += -mi*mj/r;
	    if (r < rmin) {rmin = r; bmin = bj;}
	}
}

// compute_energies: calculate the appropriate color code for particle bi
//		     relative to "particle" bj (usually the root node).

void compute_energies(dyn* bj, dyn* bi, char& c)
{
    c = default_color;

    real   epot = 0, rmin = VERY_LARGE_NUMBER;
    dynptr bmin = bi;

    accumulate_potential_energy(bj, bi, epot, rmin, bmin);

    real   mi = bi->get_mass();
    real   mj = bmin->get_mass();
    vec d_vel = bi->get_vel() - bmin->get_vel();
    vec d_pos = bi->get_pos() - bmin->get_pos();
    real   r = sqrt(d_pos * d_pos);
    real   e0 = (0.5 * d_vel * d_vel - (mi + mj)/r);

    if (e0 < 0) {
	real   epot1 = 0, rmin1 = VERY_LARGE_NUMBER;
	dynptr bmin1 = bi;

	accumulate_potential_energy(bj, bmin, epot1, rmin1, bmin1);

	if (bi == bmin1) {
	    real e  = - mi*mj / r;
	    vec R_vel = (mi*bi->get_vel()+mj*bmin->get_vel())/(mi+mj);
	    real ekin = 0.5*(mi+mj)*R_vel*R_vel;

	    if (epot + epot1 - 2*e + ekin < 0) c = bound_binary;
	    else c = unbound_binary;

	} else c = bound_single;

    } else {
	vec vel = bi->get_vel();
	real ekin = 0.5*bi->get_mass()*vel*vel;
	
	if (ekin + epot > 0.0) c = unbound_single;
    }
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


static float linespacing = 1.0;

void show_instructions(unsigned long  win, float r, char* buffer, int update)
{
    float statusx, statusy; 
    char s[255];
    int ptr = 0, len, line = 0;

    lux_clear_window(win);
    lux_reconvert_rcoord(win,0,0,&statusx, &statusy);
    istringstream ins(buffer);
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
    istringstream ins(buffer);
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

local void format_and_show_instructions(unsigned long co, float r,
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
    if (b_flag) {

	// Reverse video color scheme

	// Color by energy:

	ce[background_color] = ci[background_color]
	    = lux_rgb_pixel(win, 0.0,0.0,0.0);        		// white
	ce[default_color] = ci[default_color]
	    = lux_rgb_pixel(win, 1.0,1.0,1.0);           	// black
    
	ce[bound_single] = lux_lookup_color(win, "deepskyblue");
	ce[bound_binary] = lux_lookup_color(win, "green");
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
    
    } else {

	// Color by energy:

	ce[background_color] = ci[background_color]
	    = lux_rgb_pixel(win, 1.0,1.0,1.0);        		// white
	ce[default_color] = ci[default_color]
	    = lux_rgb_pixel(win, 0.0,0.0,0.0);           	// black
    
	ce[bound_single] = lux_lookup_color(win, "blue");
	ce[bound_binary] = lux_lookup_color(win, "limegreen");
	ce[unbound_single] = lux_lookup_color(win, "red");
	ce[unbound_binary] = lux_lookup_color(win, "goldenrod");

	// Color by index:
    
	ci[2] = lux_lookup_color(win, NV_COLOR2);
	ci[3] = lux_lookup_color(win, NV_COLOR3);
	ci[4] = lux_lookup_color(win, NV_COLOR4);
	ci[5] = lux_lookup_color(win, NV_COLOR5);
	ci[6] = lux_lookup_color(win, NV_COLOR6);
	ci[7] = lux_lookup_color(win, NV_COLOR7);
	ci[8] = lux_lookup_color(win, NV_COLOR8);
	ci[9]  = lux_lookup_color(win, NV_COLOR9);
	ci[10] = lux_lookup_color(win, NV_COLORa);
	ci[11] = lux_lookup_color(win, NV_COLORb);
	ci[12] = lux_lookup_color(win, NV_COLORc);
	ci[13] = lux_lookup_color(win, NV_COLORd);
	ci[14] = lux_lookup_color(win, NV_COLORe);
	ci[15] = lux_lookup_color(win, NV_COLORf);
	ci[16] = lux_lookup_color(win, NV_COLORg);
    
    }
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

local int nearest_index(dyn* b, float r, float s, int kx, int ky)
{
    real r2_min = VERY_LARGE_NUMBER;
    int index = -1;
    for_all_leaves(dyn, b, bb) {
	real dx = bb->get_pos()[kx] - r;
	real dy = bb->get_pos()[ky] - s;
	real r2 = dx*dx + dy*dy;
	if (r2 < r2_min) {
	    r2_min = r2;
	    index = bb->get_index();
	}
    }
    return index;
}

local int count_nodes(dyn* b)
// Return the total number of nodes BELOW the node b.
{
    if(b->get_oldest_daughter() == NULL)
	return 0;
    else {
	int  n = 0;
	for (dyn* daughter = b->get_oldest_daughter();
	     daughter != NULL; daughter = daughter->get_younger_sister())
	    n += 1 + count_nodes(daughter);
	return n;
    }
}

//-------------------------------------------------------------------------

// These shouldn't really be global, but keep them this way for now...

static int kx = 1, ky = 2; 

static unsigned long  win, dia, instr, colwin;
static int    status = 0;

static float  xmin, xmax, ymin, ymax, point_size, lmax3d, origin[3];

// The plotting box is determined completely by origin and lmax3d.
// Its projection onto the viewing plane is xmin, xmax, ymin, ymax

static float  theta = 0.33, costheta = cos(0.33), sintheta = sin(0.33);
static float  dtheta = 0.03;
static float  phi = 0.33, cosphi = cos(0.33), sinphi = sin(0.33);
static int    graph3d = 1;
static int    nodes = 0, links = 0, root = 0;
static char   track = 0, cenergy = 0;
static char   c = 0;
static char   myrotate = 0;
static float  delay_time;
static unsigned long c_energy[10], c_index[N_COLORS+1];

static real local_offset[3];
static int  origin_star;

//-------------------------------------------------------------------------

local void convert_relative_to_absolute(dyn* b)
{
    if (b->get_parent()) b->inc_pos(b->get_parent()->get_pos());
    for_all_daughters(dyn, b, bb) convert_relative_to_absolute(bb);
}

local void draw_star_point(unsigned long win, float r, float s,
			   float point_size, bool f_flag)

// Plot a point of size point_size at (r,s), filled or unfilled
// according to the value of f_flag.

{
    if (f_flag)
	lux_fill_arcf(win, r-point_size/2, s-point_size/2,
		      point_size, point_size, 0.0, 360.0);
    else
	lux_draw_arcf(win, r-point_size/2, s-point_size/2, 
		      point_size,point_size, 0.0, 360.0);
}

local float interp_to_x(float r, float s, float rr, float ss, float x)
{
    return s + (ss - s) * (x - r) / (rr - r);
}

local float interp_to_y(float r, float s, float rr, float ss, float y)
{
    return r + (rr - r) * (y - s) / (ss - s);
}

local void draw_links_2d(dyn* b, float r, float s)
{
    for_all_daughters(dyn, b, bb) {

	float rr = bb->get_pos()[kx]-local_offset[kx];
	float ss = bb->get_pos()[ky]-local_offset[ky];
	float stemp;

	// Deal with the 9 possible locations of (rr,ss) with respect
	// to the box defined by xmin, xmax, ymin, and ymax.

	if (rr < xmin) {

	    if (ss < ymin) {
		float stemp = interp_to_x(r, s, rr, ss, xmin);
		if (stemp >= ymin)
		    rr = xmin, ss = stemp;
		else
		    rr = interp_to_y(r, s, rr, ss, ymin), ss = ymin;
	    } else if (ss > ymax) {
		if ((stemp = interp_to_x(r, s, rr, ss, xmin)) <= ymax)
		    rr = xmin, ss = stemp;
		else
		    rr = interp_to_y(r, s, rr, ss, ymax), ss = ymax;
	    } else
		ss = interp_to_x(r, s, rr, ss, xmin), rr = xmin;

	} else if (rr > xmax) {

	    if (ss < ymin) {
		if ((stemp = interp_to_x(r, s, rr, ss, xmax)) >= ymin)
		    rr = xmax, ss = stemp;
		else
		    rr = interp_to_y(r, s, rr, ss, ymin), ss = ymin;
	    } else if (ss > ymax) {
		if ((stemp = interp_to_x(r, s, rr, ss, xmax)) <= ymax)
		    rr = xmax, ss = stemp;
		else
		    rr = interp_to_y(r, s, rr, ss, ymax), ss = ymax;
	    } else
		ss = interp_to_x(r, s, rr, ss, xmax), rr = xmax;

	} else {

	    if (ss < ymin)
		rr = interp_to_y(r, s, rr, ss, ymin), ss = ymin;
	    else if (ss > ymax)
		rr = interp_to_y(r, s, rr, ss, ymax), ss = ymax;

	}
	lux_draw_linef(win, r, s, rr, ss);
    }
}					   
local void draw_links_3d(dyn* b, float r, float s)
{
    for_all_daughters(dyn, b, bb) {

	float X = (float)bb->get_pos()[0] - local_offset[0] - origin[0];
	float Y = (float)bb->get_pos()[1] - local_offset[1] - origin[1];
	float Z = (float)bb->get_pos()[2] - local_offset[2] - origin[2];

	// Don't attempt to deal with the 27 possible locations
	// of (X, Y, Z) with respect to the box!

	float rr, ss;
	project3d(X, Y, Z, rr, ss, costheta, sintheta, cosphi, sinphi);
	lux_draw_linef(win, r, s, rr, ss);
    }
}					   

local int xplot_stars(dyn* b, int k, int f_flag)
{
    int  n_stars = b->n_leaves();
    int  n_nodes = count_nodes(b);
    float r, s;

    if (graph3d) {

	for_all_nodes(dyn, b, bi) if (root || bi->get_parent()) {
	    if (nodes || (bi->get_oldest_daughter() == NULL)) {
		
		float X = (float)bi->get_pos()[0] - local_offset[0] - origin[0];
		float Y = (float)bi->get_pos()[1] - local_offset[1] - origin[1];
		float Z = (float)bi->get_pos()[2] - local_offset[2] - origin[2];

		if ( (X > (-lmax3d+point_size)) && (X < (lmax3d-point_size)) 
		    && (Y > (-lmax3d+point_size)) && (Y < (lmax3d-point_size)) 
		    && (Z > (-lmax3d+point_size)) && (Z < (lmax3d-point_size))){

		    // Should really sort by depth here...

		    project3d(X, Y, Z, r, s, costheta, sintheta,
			      cosphi, sinphi);

		    // Determine the color to use.

		    bool temp_flag = f_flag;
		    if (bi->get_oldest_daughter()) {
			lux_set_color(win, lux_lookup_color(win, "grey"));
			temp_flag =  1 - f_flag;
		    } else {
			if (cenergy) { 
			    compute_energies(b, bi, c);
			    lux_set_color(win,c_energy[c]);
			} else if (bi->get_index() > 0) {

			    // Wrap the color map.

			    int ii = bi->get_index();
			    while (ii > N_COLORS) ii -= N_COLORS;

			    lux_set_color(win,c_index[ii]);
			}
		    }
		    
		    if (track) 
			lux_draw_pointf(win, r, s);
		    else
			draw_star_point(win, r, s, point_size, f_flag);

		    if (links && bi->get_oldest_daughter())
			draw_links_3d(bi, r, s);
		}
	    }
	}

    } else {

	// Make a list of nodes, *including* the root node if root = 1.
	// Root will be at the start of the list if it is being displayed.

	dynptr* p = new dynptr[n_nodes+root];

	int ip = 0;
	for_all_nodes(dyn, b, bi) if (root || bi != b) p[ip++] = bi;
	if (ip != n_nodes+root) {
	    cerr << "xplot_stars: n_nodes = " << n_nodes+root
		<< " counted " << ip << endl;
	    exit(1);
	}

	// Sort by k (note that k = 1, 2, or 3 for x, y, or z):

	for (ip = 0; ip < n_nodes+root; ip++)
	    for (int jp = ip+1; jp < n_nodes+root; jp++)
		if (p[jp]->get_pos()[k-1] < p[ip]->get_pos()[k-1]) {
		    dynptr bb = p[jp];
		    p[jp] = p[ip];
		    p[ip] = bb;
		}

	// Plot ordered by depth.

	for (ip = 0; ip < n_nodes+root; ip++) {
	    dyn* bi = p[ip];
	    if ( (root || bi != b)
		&& (nodes || bi->get_oldest_daughter() == NULL)) {

		r = (float)bi->get_pos()[kx] - local_offset[kx];
		s = (float)bi->get_pos()[ky] - local_offset[ky];

		if ( (r > (xmin+point_size)) && (r < (xmax-point_size)) 
		    && (s > (ymin+point_size)) && (s < (ymax-point_size)) ) {

		    // Determine the color to use.

		    bool temp_flag = f_flag;
		    if (bi->get_oldest_daughter()) {
			lux_set_color(win, lux_lookup_color(win, "grey"));
			temp_flag =  1 - f_flag;
		    } else {
			if (cenergy) { 
			    compute_energies(b, bi, c);
			    lux_set_color(win, c_energy[c]);
			} else if (bi->get_index() > 0) {

			    // Wrap the color map.

			    int ii = bi->get_index();
			    while (ii > N_COLORS) ii -= N_COLORS;

			    lux_set_color(win,c_index[ii]);
			}
		    }

		    if (track) 
			lux_draw_pointf(win, r, s);
		    else
			draw_star_point(win, r, s, point_size, temp_flag);

		    if (links && bi->get_oldest_daughter())
			draw_links_2d(bi, r, s);
		}
	    }
	}
	delete p;
    }
    return n_stars;
}

//-------------------------------------------------------------------------

// Convenient scaling macro:

#define _R_(i)  ((int)( ((float)i)*r_factor + 0.5 ))

/*-----------------------------------------------------------------------------
 *  xstarplot.C  --  project the positions of all particles onto the screen
 *                   input:   pn: a pointer to a nbody system,
 *		                k: the number of the coordinate axis along which
 *		                   the projection is directed, with {k,x,y}
 *		                   having a right-handed orientation,
 *                           xmax: displayed x-axis spans [-xmax, xmax]
 *                           ymax: displayed y-axis spans [-ymax, ymax]
 *-----------------------------------------------------------------------------
 */

void  xstarplot(dyn* b, float r_factor, int& k, int d, float lmax,
		float size, float D, int ce,
		bool b_flag, bool f_flag, bool t_flag,
		int init_flag)
{
    int xorigin, yorigin;
    static int win_size;
    float r, s;

    if (init_flag == 0) {

	// Open windows for plot and instructions.

	if (t_flag) nodes = links = 1;

	if (!status) {

	    // Sort of a kludge to try to get the font size right...
	    // Note that the same font is used everywhere.

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

	    // twm title bar is ~27-30 pixels (independent of scaling), but
	    // the bar doesn't count in the location of the window!
	    // Mac windows need more space (frame at bottom too).
	    // Assume that scale factors less than 0.75 mean we are
	    // using a Mac display.  For r_factor > 1.1, move the 
	    // "colwin" box up and over.

	    int title = 27;
	    int bottom = 0;
	    if (r_factor < 0.75) {
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

	    // ???

	    if (r_factor <= 1) 
		yorigin += ysize + title + bottom;
	    else
		yorigin += ysize + 2*title + bottom;

	    ysize = _R_(260);
	    instr = lux_openwin(xorigin, yorigin, xsize, ysize);
	    lux_set_window_name(instr, "Status & Commands");

	    if (r_factor <= 1.1) {
		yorigin += ysize + title + bottom;
	    } else {
		int edge = 3;
		xorigin += xsize + edge;
		yorigin = 0;
	    }
	    ysize = _R_(140);
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

	    if (r_factor <= 1.1) {
		int edge = 3;
		xorigin += xsize + edge;
		yorigin = 0;
	    } else
		yorigin += ysize + title + bottom;
	}

	lux_clear_window(win);
	lux_update_fg(win);
	lux_clear_window(colwin);
	lux_update_fg(colwin);
	lux_clear_window(instr);
	lux_update_fg(instr);

	// Determine axes to plot:

	switch (k) {
	    case 1:	kx = 1; ky = 2; graph3d = 0; break;
	    case 2:	kx = 2; ky = 0; graph3d = 0; break;
	    case 3:	kx = 0; ky = 1; graph3d = 0; break;
	    default: 	cerr << "xstarplot: k = " << k
		             << ": illegal value; choose from {1, 2, 3}"
			     << endl;
			exit(0);
	}

	switch(d) {
	    case 2:     graph3d = 0;  break;
	    case 3:     graph3d = 1;  break;
	    default:    cerr << "xstarplot: d = " << d
		             << " illegal value; choose from {2, 3)"
			     << endl;
			exit(0);
	}

	cenergy = ce?1:0;
	show_color_scheme(colwin, c_energy, c_index,
			  r_factor, cenergy, b_flag,1);

	delay_time = D;

	if (lmax > 0.0)

	    lmax3d = lmax;

	else {

	    lmax3d = 0;

	    // Use all dimensions in determining the initial scale!

	    for_all_leaves(dyn, b, bi)
		for (int kk = 0; kk < 3; kk++)
		    lmax3d = Starlab::max(lmax3d, abs(bi->get_pos()[kk]));

	    // Round lmax3d up to something reasonable:

	    lmax3d *= 1.2;
	    if (lmax3d <= 0) lmax3d = 1;

	    real m_log = log10(lmax3d);
	    int i_log = (int) m_log;
	    if (m_log - i_log > 0.699) i_log++;
	    real scale = pow(10.0, i_log);
	    lmax3d = ((int) (lmax3d/scale) + 1) * scale;
	}

	for (int ko = 0; ko < 3; ko++) origin[ko] = 0;

	set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

	point_size = lmax3d/30.0;
	if (size > 0.0) point_size *= size;
	
	if (graph3d) {
	    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
			   -FAC3D*lmax3d, FAC3D*lmax3d);      
	    draw3d_axis(win, lmax3d, costheta, sintheta, cosphi, sinphi);
	} else {
	    lux_setup_axis(win, xmin, xmax, ymin, ymax);
	    draw2d_axis(win, xmin, xmax, ymin, ymax, k);
	}

	if (!status) {

	    int i = 0;
// 	    cout << "setting dialog...\n";

	    // Initialize dialog stuff (note: y is measured up from bottom!)

	    real save = r_factor;
	    if (r_factor > 1) r_factor = 1;
	    dia = lux_open_dialog(xorigin, yorigin, _R_(512), _R_(550));

//	    cout << i++ << endl;
	    // 3-D origin:

	    lux_set_item(dia, Origin, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(500), 6, "origin");

//	    cout << i++ << endl;
	    sprintf(temp_buffer, "%f", origin[0]);
	    lux_set_item(dia, Xorigin, INPUT_WINDOW, NO_TYPE,
			 _R_(160), _R_(500), 10, temp_buffer);
//	    cout << i++ << endl;
	    sprintf(temp_buffer, "%f", origin[1]);
	    lux_set_item(dia, Yorigin, INPUT_WINDOW, NO_TYPE,
			 _R_(280), _R_(500), 10, temp_buffer);
//	    cout << i++ << endl;
	    sprintf(temp_buffer, "%f", origin[2]);
	    lux_set_item(dia, Zorigin, INPUT_WINDOW, NO_TYPE,
			 _R_(400), _R_(500), 10, temp_buffer);

//	    exit(1);
//	    cout << i++ << endl;
	    // Minima and maxima of 2-D display:

	    lux_set_item(dia, xminimum, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(460), 4, "xmin");
	    lux_set_item(dia, xmaximum, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(420), 4, "xmax");
	    lux_set_item(dia, yminimum, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(380), 4, "ymin");
	    lux_set_item(dia, yminimum, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(340), 4, "ymax");

	    sprintf(temp_buffer, "%f", xmin);
	    lux_set_item(dia, xminimum, INPUT_WINDOW, NO_TYPE,
			 _R_(160), _R_(460), 10, temp_buffer);
	    sprintf(temp_buffer, "%f", xmax);
	    lux_set_item(dia, xmaximum, INPUT_WINDOW, NO_TYPE,
			 _R_(160), _R_(420), 10, temp_buffer);
	    sprintf(temp_buffer, "%f", ymin);
	    lux_set_item(dia, yminimum, INPUT_WINDOW, NO_TYPE,
			 _R_(160), _R_(380), 10, temp_buffer);
	    sprintf(temp_buffer, "%f", ymax);
	    lux_set_item(dia, ymaximum, INPUT_WINDOW, NO_TYPE,
			 _R_(160), _R_(340), 10, temp_buffer);

//	    cout << i++ << endl;
	    // 3-D cube size:

	    lux_set_item(dia, lmax3D, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(300), 10, "cube size");
	    sprintf(temp_buffer, "%f", 2*lmax3d);
	    lux_set_item(dia, lmax3D, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(300), 10, temp_buffer);

//	    cout << i++ << endl;
	    // 3-D projection direction:

	    lux_set_item(dia, theta3D, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(260), 10, "theta");
	    sprintf(temp_buffer, "%f", theta);
	    lux_set_item(dia, theta3D, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(260), 10, temp_buffer);
	    lux_set_item(dia, phi3D, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(220), 10, "phi");
	    sprintf(temp_buffer, "%f", phi);
	    lux_set_item(dia, phi3D, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(220), 10, temp_buffer);
	    lux_set_item(dia, dtheta3D, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(180), 10, "dtheta");
	    sprintf(temp_buffer, "%f", dtheta);
	    lux_set_item(dia, dtheta3D, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(180), 10, temp_buffer);

//	    cout << i++ << endl;
	    // Point size:

	    lux_set_item(dia, pointsize, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(140), 10, "point size");
	    sprintf(temp_buffer, "%f", point_size);
	    lux_set_item(dia, pointsize, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(140), 10, temp_buffer);

//	    cout << i++ << endl;
	    // Delay time:

	    lux_set_item(dia, DelayTime, TEXT_WINDOW, NO_TYPE,
			 _R_(70), _R_(100), 10, "delay time");
	    sprintf(temp_buffer, "%f", delay_time);
	    lux_set_item(dia, DelayTime, INPUT_WINDOW, NO_TYPE,
			 _R_(180), _R_(100), 10, temp_buffer);

//	    cout << i++ << endl;
	    // 2-D view axis:

	    sprintf(temp_buffer, " %d ", k);
	    lux_set_item(dia, view2D, INPUT_WINDOW, NO_TYPE,
			 _R_(320), _R_(460), 3, temp_buffer);
	    lux_set_item(dia, View2D, TEXT_WINDOW, NO_TYPE,
			 _R_(360), _R_(460), 8, "2-D view axis");

//	    cout << i++ << endl;
	    // Energy/tracking/3D selection boxes:

	    temp_buffer[0] = cenergy;
	    lux_set_item(dia, colorenergy, BUTTON_WINDOW, CHECK_BUTTON,
			 _R_(320), _R_(420), 1, temp_buffer);
	    temp_buffer[0] = track;
	    lux_set_item(dia, tracking, BUTTON_WINDOW, CHECK_BUTTON,
			 _R_(320), _R_(380), 1, temp_buffer);
	    temp_buffer[0] = graph3d;
	    lux_set_item(dia, graph3dim, BUTTON_WINDOW, CHECK_BUTTON,
			 _R_(320), _R_(340), 1, temp_buffer);

	    lux_set_item(dia, colorenergy, TEXT_WINDOW, CHECK_BUTTON,
			 _R_(350), _R_(420), 11, "color by energy");
	    lux_set_item(dia, tracking, TEXT_WINDOW, CHECK_BUTTON,
			 _R_(350), _R_(380), 14, "track particle");
	    lux_set_item(dia, graph3dim, TEXT_WINDOW, CHECK_BUTTON,
			 _R_(350), _R_(340), 8, "3-D graph");

//	    cout << i++ << endl;
	    // Origin star:

	    sprintf(temp_buffer, " %d ", 0);
	    lux_set_item(dia, originstar, INPUT_WINDOW, NO_TYPE,
			 _R_(320), _R_(300), 3, temp_buffer);
	    lux_set_item(dia, OriginStar, TEXT_WINDOW, NO_TYPE,
			 _R_(360), _R_(300), 8, "Origin Star");

//	    lux_draw_palette(dia);

	    status = 1;

//	    cout << i++ << endl;
	    // OK/cancel:

	    lux_set_item(dia, ok, BUTTON_WINDOW, OK_BUTTON,
			 _R_(100), _R_(25), 9, "OK, CLEAR");
	    lux_set_item(dia, ok_keep, BUTTON_WINDOW, OK_KEEP_BUTTON,
			 _R_(250), _R_(25), 8, "OK, KEEP");
	    lux_set_item(dia, cancel, BUTTON_WINDOW, CANCEL_BUTTON,
			 _R_(400), _R_(25), 6, "CANCEL");

	    r_factor = save;

	} else ;  // Here we should update the all the values for dialog
	    
    } else if (init_flag < 0) {
	show_instructions(instr, r_factor,
		  "End of Data.  Idle now\n  r: replay, (c,q): quit", 1);
	lux_getevent();
	exit(0);
    }

//    lux_getevent();
//    exit(1);
    

//----------------------------------------------------------------------------//

    // End of initialization.  Start by redrawing axes and instructions(?)

    if (!track) {
	lux_clear_current_region(win);
	if (graph3d) draw3d_axis(win, lmax3d, costheta, sintheta,
				 cosphi, sinphi);
    }

    show_main_instructions(instr, r_factor, graph3d, 1);

//----------------------------------------------------------------------------//

    // Always express all positions and velocities relative to the root node.
    // (There is no requirement that the root node be at rest at the origin...)

    {
	vec root_pos = b->get_pos();
	vec root_vel = b->get_vel();
	for_all_nodes(dyn, b, bi) {
	    bi->inc_pos(-root_pos);
	    bi->inc_vel(-root_vel);
	}
    }

    // Determine local offset, if any.

    if (origin_star > 0) {
	for_all_leaves(dyn, b, bi)
	    if (bi->get_index() == origin_star) {
		for (int kk = 0; kk < 3; kk++)
		    local_offset[kk] = bi->get_pos()[kk];
		break;
	    }
    } else
	for (int kk = 0; kk < 3; kk++) local_offset[kk] = 0;

    // Plot the data points (sorted by k-component in the 2D case):

    int n_stars = xplot_stars(b, k, f_flag);

//----------------------------------------------------------------------------//

    // Headers, etc.

    lux_set_color(win,c_energy[default_color]);

    // Other labels:

    if (graph3d) {
	sprintf(temp_buffer, "N = %d (snapshot #%5d) max=%5.3f",
		n_stars, init_flag + 1, lmax3d);
	lux_draw_image_string(win, 0.0, lmax3d*FAC3D, 0.5, temp_buffer, 0);
    } else {
	sprintf(temp_buffer, "N = %d  (snapshot #%5d)", n_stars, init_flag + 1);
	lux_draw_image_string(win, (xmin+xmax)/2.0, ymax, 0.5, temp_buffer, 0);
    }
    update_with_delay(win, delay_time);

//----------------------------------------------------------------------------//

    // Loop through the input temp_buffer until no more events remain.
    // lux_next_keypress gets and discards the next key in the input stream.

    char key, string[20], shift, control;
    while (lux_next_keypress(win, &key, string, &shift, &control)) {

        // A "defined" key has been pressed.  See if it is one we want.

        if (key == 0			// key = 0 for non-ASCII character
    			 		// e.g. Up, Down, Home, etc.--see win.c
	    || key == 'h'
	    || key == 'a') {

	    // "Shift" functions:

	    int change = 0;

	    if (key == 'h')		// Keyboard lookalikes for
		change = 4;		// function keys.
	    else if (key == 'a')
		change = 5;

	    else if (strcmp(string, "Up") == 0) {
		change = ky + 1;
		origin[ky] += lmax3d/2;
	    } else if (strcmp(string, "Down") == 0) {
		change = ky + 1;
		origin[ky] -= lmax3d/2;
	    } else if (strcmp(string, "Right") == 0) {
		change = kx + 1;
		origin[kx] += lmax3d/2;
	    } else if (strcmp(string, "Left") == 0) {
		change = kx + 1;
		origin[kx] -= lmax3d/2;
	    } else if (strcmp(string, "PgUp") == 0) {
		change = k;
		origin[k-1] += lmax3d/2;
	    } else if (strcmp(string, "PgDn") == 0) {
		change = k;
		origin[k-1] -= lmax3d/2;
	    } else if (strcmp(string, "Home") == 0)
		change = 4;
	    else if (strcmp(string, "R11") == 0)
		change = 5;

	    if (change) {

		if (change == 4) {
		    origin[kx] = origin[ky] = origin[k-1] = 0;
		} else if (change == 5) {
		    lmax3d = 0;

		    // Use all dimensions in redetermining the scale!

		    for_all_leaves(dyn, b, bi)
			for (int kk = 0; kk < 3; kk++)
			    lmax3d = Starlab::max(lmax3d, abs(bi->get_pos()[kk]
						     - local_offset[kk]));
		
		    // Round lmax3d up to something reasonable:
		
		    lmax3d *= 1.2;
		    if (lmax3d <= 0) lmax3d = 1;
		
		    real m_log = log10(lmax3d);
		    int i_log = (int) m_log;
		    if (m_log - i_log > 0.699) i_log++;
		    real scale = pow(10.0, i_log);
		    lmax3d = ((int) (lmax3d/scale) + 1) * scale;

		    origin[kx] = origin[ky] = origin[k-1] = 0;		
		    point_size  = lmax3d/30.0;
		    if (size > 0.0) point_size *= size;
		}

		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

		// Update all entries in the dialog box:

		if (change == 1 || change >= 4) {
		    sprintf(temp_buffer,"%f",origin[0]);
		    lux_update_itemvalue(dia, Xorigin,
					 INPUT_WINDOW, NO_TYPE, temp_buffer);
		}
		if (change == 2 || change >= 4) {
		    sprintf(temp_buffer,"%f",origin[1]);
		    lux_update_itemvalue(dia, Yorigin,
					 INPUT_WINDOW, NO_TYPE, temp_buffer);
		}
		if (change == 3 || change >= 4) {
		    sprintf(temp_buffer,"%f",origin[2]);
		    lux_update_itemvalue(dia, Zorigin,
					 INPUT_WINDOW, NO_TYPE, temp_buffer);
		}

		if (change == 5) {
		    sprintf(temp_buffer,"%f",point_size);
		    lux_update_itemvalue(dia, pointsize, INPUT_WINDOW,
					 NO_TYPE, temp_buffer);
		    sprintf(temp_buffer,"%f",2*lmax3d);
		    lux_update_itemvalue(dia, lmax3D,
					 INPUT_WINDOW,NO_TYPE,temp_buffer);
		}

		sprintf(temp_buffer,"%f",xmin);
		lux_update_itemvalue(dia,xminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmax);
		lux_update_itemvalue(dia,xmaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymin);
		lux_update_itemvalue(dia,yminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymax);
		lux_update_itemvalue(dia,ymaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

		// No redrawing of axes is needed for 3D graphs for change != 5.
		
		if (!graph3d) {

		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		} else if (change == 5) {
		    lux_clear_window(win);
		    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		}
	    }

	} else {			// Key is the ASCII code.

	    // Check for "quit" first.

	    if (key == 'q') {
		show_instructions(instr, r_factor,
		    "\n   Quitting...                                        ",
				  1);
		lux_exit();	// Note that this is now quick_exit
	    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//

	    if (graph3d && !track) {
		int button;

		// 3-D operations:

		if (key== 'R') {

		    // Handle static rotation using lux_check_keypress:

		    show_instructions(instr, r_factor,
		      "Pause...\n  r: rotate, (p,c): continue                ",
				      1);

		    do {
			if (lux_check_keypress(win,'r')) {
			    show_instructions(instr, r_factor,
		      "Rotate... p: pause, c: continue, c: exit rotation     ",
					      0);
			    myrotate = 1;

			    do {
				theta = theta + dtheta;
				if (theta < 0.0) theta += 2*PI;
				else if (theta > PI*2) theta -= 2*PI;
				sprintf(temp_buffer, "%f", theta);
				lux_update_itemvalue(dia, theta3D, INPUT_WINDOW,
						     NO_TYPE, temp_buffer);
				costheta = cos(theta);
				sintheta = sin(theta);
				lux_clear_current_region(win);
				draw3d_axis(win, lmax3d, costheta, sintheta,
					    cosphi, sinphi);
				for_all_leaves(dyn, b, bi) {
				    float X, Y, Z;
				    X = (float)bi->get_pos()[0] - origin[0];
				    Y = (float)bi->get_pos()[1] - origin[1];
				    Z = (float)bi->get_pos()[2] - origin[2];
				    if ( (X > (-lmax3d+point_size))
					&& (X < (lmax3d-point_size)) 
					&& (Y > (-lmax3d+point_size))
					&& (Y < (lmax3d-point_size)) 
					&& (Z > (-lmax3d+point_size))
					&& (Z < (lmax3d-point_size))) {
					project3d( X, Y, Z, r, s,
						  costheta, sintheta,
						  cosphi, sinphi);

					if (cenergy) { 
					    compute_energies(b, bi, c);
					    lux_set_color(win,c_energy[c]);
					} else if (bi->get_index()>0) {
					    if (bi->get_index() <= N_COLORS)
						lux_set_color(win,
						    c_index[bi->get_index()]);
					    else lux_set_color(win,c_index[1]);
					}

					if (track) 
					    lux_draw_pointf(win,r,s);
					else
					    if (f_flag)
						lux_fill_arcf(win,
							      r-point_size/2,
							      s-point_size/2, 
							      point_size,
							      point_size,
							      0.0, 360.0);
					    else
						lux_draw_arcf(win,
							      r-point_size/2,
							      s-point_size/2, 
							      point_size,
							      point_size,
							      0.0, 360.0);
				    }
				}
				update_with_delay(win, delay_time);

				if (lux_check_keypress(win,'p')) 
				    while(!lux_check_keypress(win,'c'));

			    } while(!lux_check_keypress(win,'c'));

			    // Flush "c"s from input stream:

			    while(lux_check_keypress(win,'c'));
			}
		    } while(!lux_check_keypress(win,'p')
			    && !lux_check_keypress(win,'c') && !myrotate);

		    myrotate = 0;
		    lux_set_color(win,c_energy[default_color]);

		    // *** End of static rotation loop. ***

		} else {

		    button = lux_check_buttonpress(win);

		    if (key == '^') {
			phi = phi - dtheta;
			if (phi < -PI) phi = phi + 2.0*PI;
			cosphi = cos(phi); sinphi = sin(phi);
			sprintf(temp_buffer, "%f", phi);
			lux_update_itemvalue(dia, phi3D, INPUT_WINDOW,
					     NO_TYPE, temp_buffer);
		    } else if (key == 'V') {
			phi = phi + dtheta;
			if (phi > PI) phi = phi - 2.0*PI;
			cosphi = cos(phi); sinphi = sin(phi);
			sprintf(temp_buffer, "%f", phi);
			lux_update_itemvalue(dia, phi3D, INPUT_WINDOW,
					     NO_TYPE, temp_buffer);
		    } else if (key == '<') {
			theta = theta + dtheta;
			if (theta > 2*PI) theta -= 2*PI;
			sprintf(temp_buffer, "%f", theta);
			lux_update_itemvalue(dia, theta3D, INPUT_WINDOW,
					     NO_TYPE, temp_buffer);
			costheta = cos(theta);
			sintheta = sin(theta);
		    } else if (key == '>') {
			theta = theta - dtheta;
			if (theta < 0.0) theta += 2*PI;
			sprintf(temp_buffer, "%f", theta);
			lux_update_itemvalue(dia, theta3D, INPUT_WINDOW,
					     NO_TYPE, temp_buffer);
			costheta = cos(theta);sintheta = sin(theta);
		    }
		}
	    }		// End of 3-D.

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//

	    // Modify point size:

	    if (key == 'p' || key == 'P') {
		if (key == 'p') {
		    if (point_size > 4*lmax3d/win_size) point_size /= PFAC;
		} else
		    point_size *= PFAC;
		sprintf(temp_buffer,"%f",point_size);
		lux_update_itemvalue(dia, pointsize, INPUT_WINDOW,
				     NO_TYPE, temp_buffer);
	    }

	    // Toggle tree display:

	    if (key == 'n') {
		nodes = 1 - nodes;
		if (nodes == 0) links = 0;
	    }

	    if (key == 'l') {
		links = 1 - links;
		if (links == 1) nodes = 1;
	    }

	    if (key == 'r') root = 1 - root;

	    // Toggle tracking:

	    if (key == 't') {
		track = 1 - track;
		temp_buffer[0] = track;
		lux_update_itemvalue(dia, tracking, BUTTON_WINDOW, CHECK_BUTTON,
				     temp_buffer);
	    }

	    // Toggle 2D/3D mode:

	    if (key == '3' && (!graph3d)) {
		graph3d = 1;
		lux_update_itemvalue(dia, graph3dim, BUTTON_WINDOW,
				     CHECK_BUTTON,"\1");
		lux_clear_window(win);
		lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
			       -FAC3D*lmax3d, FAC3D*lmax3d);
		draw3d_axis(win, lmax3d, costheta, sintheta, cosphi, sinphi);
	    } else if (graph3d &&
		       (key == '2' || key == 'x' || key == 'y' || key == 'k')) {
		graph3d = 0;
		lux_update_itemvalue(dia, graph3dim, BUTTON_WINDOW,
				     CHECK_BUTTON,"\0");
		lux_clear_window(win);
		lux_setup_axis(win, xmin, xmax, ymin, ymax);
		draw2d_axis(win, xmin, xmax, ymin, ymax, k);
	    }

	    if (key == 'e') {
		temp_buffer[0] = cenergy = !cenergy;
		lux_update_itemvalue(dia, colorenergy, BUTTON_WINDOW,
				     CHECK_BUTTON,temp_buffer);
		show_color_scheme(colwin, c_energy, c_index,
				  r_factor, cenergy, b_flag, 1);
	    }

	    // Update speed:

	    if (key == '+' || key == '=') {
		delay_time = delay_time + 0.05;
		sprintf(temp_buffer, "%f",delay_time);
		lux_update_itemvalue(dia, DelayTime, INPUT_WINDOW,
				     NO_TYPE,temp_buffer);
	    }else if (key == '-') {
		delay_time = delay_time - 0.05;
		if (delay_time < 0.0) delay_time = 0.0;
		sprintf(temp_buffer, "%f",delay_time);
		lux_update_itemvalue(dia, DelayTime, INPUT_WINDOW,
				     NO_TYPE,temp_buffer);
	    } if (key == '0') {
		delay_time = 0.0;
		sprintf(temp_buffer, "%f",delay_time);
		lux_update_itemvalue(dia, DelayTime, INPUT_WINDOW,
				     NO_TYPE,temp_buffer);

		// Flush the queue:

		while(lux_check_keypress(win,'+')
		      || lux_check_keypress(win,'-'));
	    }

	    // Zoom in/out (fixed axes and origin):

	    if (key == 'z') {

		lmax3d /= ZOOM;
		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		point_size /= ZOOM;

		sprintf(temp_buffer,"%f",point_size);
		lux_update_itemvalue(dia, pointsize, INPUT_WINDOW,
				     NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",2*lmax3d);
		lux_update_itemvalue(dia, lmax3D,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmin);
		lux_update_itemvalue(dia,xminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmax);
		lux_update_itemvalue(dia,xmaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymin);
		lux_update_itemvalue(dia,yminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymax);
		lux_update_itemvalue(dia,ymaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

		if (graph3d) {
		    lux_clear_window(win);
		    lux_setup_axis(win,-FAC3D*lmax3d,FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		} else {
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);
		}

	    } else if (key == 'Z') {

		lmax3d *= ZOOM;
		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		point_size *= ZOOM;

		sprintf(temp_buffer,"%f",point_size);
		lux_update_itemvalue(dia, pointsize, INPUT_WINDOW,
				     NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",2*lmax3d);
		lux_update_itemvalue(dia, lmax3D,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmin);
		lux_update_itemvalue(dia,xminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmax);
		lux_update_itemvalue(dia,xmaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymin);
		lux_update_itemvalue(dia,yminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymax);
		lux_update_itemvalue(dia,ymaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

		if (graph3d) {
		    lux_clear_window(win);
		    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		} else {
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);
		}
	    }

	    // Change projection axis (fixed origin, scale; x=1, y=2, z=3):

	    if (key == 'k') {
		if (graph3d) {
		    theta = -PI/2; phi = PI/2.0;
		    costheta = cos(theta); sintheta = sin(theta);
		    cosphi = cos(phi); sinphi = sin(phi);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );
		} else {
		    k = 3; kx = 0; ky = 1;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		    sprintf(temp_buffer, " %d ", k);
		    lux_update_itemvalue(dia, view2D, INPUT_WINDOW, NO_TYPE,
					 temp_buffer);
		}
	    } else if (key == 'y') {
		if (graph3d) {
		    theta = -PI/2.0; phi = 0.0;
		    costheta = cos(theta); sintheta = sin(theta);
		    cosphi = cos(phi); sinphi = sin(phi);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );
		} else {
		    k = 2; kx = 2; ky = 0;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		    sprintf(temp_buffer, " %d ", k);
		    lux_update_itemvalue(dia, view2D, INPUT_WINDOW, NO_TYPE,
					 temp_buffer);
		}
	    } else if (key == 'x') {
		if (graph3d) {
		    theta = 0.0; phi = 0.0;
		    costheta = cos(theta); sintheta = sin(theta);
		    cosphi = cos(phi); sinphi = sin(phi);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );
		} else {
		    k = 1; kx = 1; ky = 2;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		    sprintf(temp_buffer, " %d ", k);
		    lux_update_itemvalue(dia, view2D, INPUT_WINDOW, NO_TYPE,
					 temp_buffer);
		}
	    }

	    // Select new origin: 'o' ==> point in space, 'O' ==> star

	    if (!graph3d && key == 'o') {
		show_instructions(instr, r_factor,
	        " Use mouse-right to select new origin...\n",1);
		float r, s;
		get_mouse_position(win, &r, &s);

		// Looks like r and s come back in the correct units!

		origin[kx] = r;
		origin[ky] = s;
		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

		lux_clear_window(win);
		lux_setup_axis(win, xmin, xmax, ymin, ymax);
		draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		sprintf(temp_buffer,"%f",origin[0]);
		lux_update_itemvalue(dia, Xorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",origin[1]);
		lux_update_itemvalue(dia, Yorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",origin[2]);
		lux_update_itemvalue(dia, Zorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",xmin);
		lux_update_itemvalue(dia,xminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmax);
		lux_update_itemvalue(dia,xmaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymin);
		lux_update_itemvalue(dia,yminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymax);
		lux_update_itemvalue(dia,ymaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

		origin_star = -1;
		sprintf(temp_buffer," %d",origin_star);
		lux_update_itemvalue(dia,originstar,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

	    } else if (!graph3d && key == 'O') {

		show_instructions(instr, r_factor,
	        " Use mouse-right to select new origin star...\n",1);
		float r, s;
		get_mouse_position(win, &r, &s);

		// Looks like r and s come back in the correct units!

		// Find the star closest to the (2-D) mouse position.

		origin_star = nearest_index(b, r, s, kx, ky);

		if (origin_star > 0) {
		    for (int kk = 0; kk < 3; kk++) origin[kk] = 0;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		}

		lux_clear_window(win);
		lux_setup_axis(win, xmin, xmax, ymin, ymax);
		draw2d_axis(win, xmin, xmax, ymin, ymax, k);

		sprintf(temp_buffer,"%f",origin[0]);
		lux_update_itemvalue(dia, Xorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",origin[1]);
		lux_update_itemvalue(dia, Yorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",origin[2]);
		lux_update_itemvalue(dia, Zorigin,
				     INPUT_WINDOW, NO_TYPE, temp_buffer);
		sprintf(temp_buffer,"%f",xmin);
		lux_update_itemvalue(dia,xminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",xmax);
		lux_update_itemvalue(dia,xmaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymin);
		lux_update_itemvalue(dia,yminimum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
		sprintf(temp_buffer,"%f",ymax);
		lux_update_itemvalue(dia,ymaximum,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);

		sprintf(temp_buffer," %d",origin_star);
		lux_update_itemvalue(dia,originstar,
				     INPUT_WINDOW,NO_TYPE,temp_buffer);
	    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//

	    // Handle dialog box input:

	    {
		int return_value;
		return_value = 0;
		while(lux_check_keypress(win,'r')); // Clean garbage

		if (key == 'i') {
		    show_instructions(instr, r_factor,
	        " Idle...\n  d: dialog, r: replay, c: continue: q-quit",1);
		    return_value = lux_getevent();
		} else if (key == 'd') {
		    show_instructions(instr, r_factor,
"Dialog mode keyboard options\n\n\
  Dialog window:\n\
    o   OK, keep dialog window\n\
    O   OK, close dialog window\n\
    c   CANCEL, keep dialog window\n\
    C   CANCEL, close dialog window\n\n\
  Other windows:\n\
    r   replay\n\
    c   continue (keep dialog window,\n\
                  ignore dialog changes)\n\
    q   quit\n\
",1);
		    lux_show_dialog(dia); return_value = lux_getevent();
		}
		
		// Get all the values if OK button being clicked (clicking
		// OK forces a return value of 3).
		
		if (return_value == 3) {

		    // Read all data from the dialog box.  Note that it is
		    // important that the box be maintained correctly.
		    // Any errors in the data displayed in the box will be
		    // propogated to the "real" data by this procedure.

		    lux_get_itemvalue(dia, view2D, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer, "%d", &k);
		    if (k == 1)
			kx = 1, ky = 2;
		    else if (k == 2)
			kx = 2, ky = 0;
		    else
			kx = 0, ky = 1;

		    lux_get_itemvalue(dia, originstar, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer, "%d", &origin_star);

		    lux_get_itemvalue(dia, Xorigin, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",origin);

		    lux_get_itemvalue(dia, Yorigin, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&origin[1]);

		    lux_get_itemvalue(dia, Zorigin, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&origin[2]);

		    lux_get_itemvalue(dia, lmax3D, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&lmax3d);
		    lmax3d /= 2;			// Display twice lmax3d

		    // ------------------------------------------------------

		    // NOTE: xmin, xmax, ymin, ymax are IRRELEVANT, since they
		    // will be forced to be consistent with origin and lmax3d.

		    lux_get_itemvalue(dia, xminimum, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&xmin);

		    lux_get_itemvalue(dia, xmaximum, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&xmax);

		    lux_get_itemvalue(dia, yminimum, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&ymin);

		    lux_get_itemvalue(dia, ymaximum, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&ymax);

		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

		    // Make sure entries in the dialog box are consistent:

		    sprintf(temp_buffer,"%f",xmin);
		    lux_update_itemvalue(dia,xminimum,
					 INPUT_WINDOW,NO_TYPE,temp_buffer);
		    sprintf(temp_buffer,"%f",xmax);
		    lux_update_itemvalue(dia,xmaximum,
					 INPUT_WINDOW,NO_TYPE,temp_buffer);
		    sprintf(temp_buffer,"%f",ymin);
		    lux_update_itemvalue(dia,yminimum,
					 INPUT_WINDOW,NO_TYPE,temp_buffer);
		    sprintf(temp_buffer,"%f",ymax);
		    lux_update_itemvalue(dia,ymaximum,
					 INPUT_WINDOW,NO_TYPE,temp_buffer);

		    // ------------------------------------------------------

		    lux_get_itemvalue(dia, pointsize, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&point_size);

		    lux_get_itemvalue(dia, theta3D, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer, "%f", &theta);
		    costheta = cos(theta); sintheta = sin(theta);

		    lux_get_itemvalue(dia, phi3D, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer, "%f", &phi);
		    cosphi = cos(phi); sinphi = sin(phi);

		    lux_get_itemvalue(dia, dtheta3D, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer, "%f", &dtheta);

		    lux_get_itemvalue(dia, DelayTime, INPUT_WINDOW, NO_TYPE,
				      temp_buffer);
		    sscanf(temp_buffer,"%f",&delay_time);

		    lux_get_itemvalue(dia, colorenergy, BUTTON_WINDOW,
				      CHECK_BUTTON, temp_buffer);
		    cenergy = temp_buffer[0];
		    show_color_scheme(colwin, c_energy, c_index,
				      r_factor, cenergy, b_flag, 1);

		    lux_get_itemvalue(dia, tracking, BUTTON_WINDOW,
				      CHECK_BUTTON,  temp_buffer);
		    track = temp_buffer[0];

		    lux_get_itemvalue(dia, graph3dim, BUTTON_WINDOW,
				      CHECK_BUTTON, temp_buffer);
		    graph3d = temp_buffer[0];

		    lux_clear_window(win);

		    if (graph3d) {
			lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
				       -FAC3D*lmax3d, FAC3D*lmax3d);
			draw3d_axis(win, lmax3d, costheta, sintheta,
				    cosphi, sinphi);
		    } else {
			lux_setup_axis(win, xmin, xmax, ymin, ymax);
			draw2d_axis(win, xmin, xmax, ymin, ymax, k);
		    }
		}
	    }
	}
    }
}

/*-----------------------------------------------------------------------------
 *  main  --  driver to use  xstarplot()  as a tool. 
 *               The argument -a is interpreted as the axis along which to
 *            view the N-body system: the number of the coordinate axis along
 *            which the projection is directed, with {k,x,y} having a
 *            right-handend orientation.
 *               If an argument -l is provided, it defines the maximum
 *	      lengths along the remaining axes.
 *-----------------------------------------------------------------------------
 */
main(int argc, char** argv)
{
    dyn *b;

    // Establish some defaults:

    int  k = 3;		// Note x = 1, y = 2, z = 3 here!
    float lmax = -1.0;
    int cenergy = 0;
    int d = 2;
    float D = 0.0;
    float rel_point_size = -1.0;
    float scale = 1.0;

    bool  b_flag = TRUE;        // if TRUE, the background is black
    bool  f_flag = TRUE;        // if TRUE, the star is solid
    bool  p_flag = FALSE;	// if TRUE, output stdin to stdout
    bool  t_flag = FALSE;	// if TRUE, show tree structure

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:d:D:efl:pP:rs:t";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
      switch(c) {

	    case 'a': k = atoi(poptarg);
		      break;
	    case 'd': d = atoi(poptarg);
		      break;
	    case 'D': D = atof(poptarg);
		      break;
	    case 'e': cenergy = 1;
		      break;
	    case 'f': f_flag = FALSE;
		      break;
	    case 'l': lmax = atof(poptarg);
		      break;
	    case 'p': p_flag = TRUE;
		      break;
	    case 'P': rel_point_size = atof(poptarg);
		      break;
	    case 'r': b_flag = FALSE;
		      break;
	    case 's': scale = atof(poptarg);
		      break;
	    case 't': t_flag = TRUE;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}
        
    int gfx_counter = 0;
    while (b = get_dyn()) {
	convert_relative_to_absolute(b);
        xstarplot(b, scale, k, d, lmax, rel_point_size, D, cenergy,
		  b_flag, f_flag, t_flag, gfx_counter++);
	if (p_flag) put_node(b);
	rmtree(b);
    }

    if (p_flag) cout << "End of data\n" << flush;

    // idle, wait to leave

    xstarplot(b, scale, k, d, lmax, rel_point_size, D, cenergy,
	      b_flag, f_flag, t_flag, -1);
}
