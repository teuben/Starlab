
#include <strstream.h>
#include "stdinc.h"

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

enum   {colorenergy=1, tracking, graph3dim,
	xminimum, xmaximum, yminimum, ymaximum,
	basepointsize, pointscalemode, lmax3D,
	theta3D, phi3D,  
	DelayTime, dtheta3D, color,
	Origin, Xorigin, Yorigin, Zorigin,
	View2D, view2D, originstar, OriginStar};
enum   {ok=1, ok_keep, cancel};

// From gfx_util.C:

void initialize_graphics(float, bool,
			 unsigned long&, unsigned long&, unsigned long&,
			 unsigned long*, unsigned long*,
			 int&, int&, int&);
void project3d(float, float, float, float&, float&,
	       float, float, float, float,
	       float, float, float, float&);
void project3d(float, float, float, float&, float&,
	       float, float, float, float);

void draw3d_axis(unsigned long, float, float, float, float, float);
void draw2d_axis(unsigned long, float, float, float, float, int);
void update_with_delay(unsigned long, float);

void show_instructions(unsigned long, float, char*, int);
void show_instructions(unsigned long, float, char*, int, int);

void show_main_instructions(unsigned long, float, int, int);
void format_and_show_instructions(unsigned long, float,
					unsigned long*, int, int,
					char*, int, int);
void show_color_scheme(unsigned long, unsigned long*, unsigned long*,
		       float, char, bool, int);
void init_colors(unsigned long, unsigned long*, unsigned long*, bool);
void set_limits(float*, float, int, float&, float&, int, float&, float&);
float interp_to_x(float, float, float, float, float);
float interp_to_y(float, float, float, float, float);

// From hrd_util.h
void initialize_hrd_graphics(float, bool,
                         unsigned long&, unsigned long&, unsigned long&,
                         unsigned long*, unsigned long*,
                         int&, int&, int&);
void show_hrd_color_scheme(unsigned long, unsigned long*, unsigned long*,
                       float, char, bool, int);

//-------------------------------------------------------------------------

// Smallest reasonable dot, assumes lmax3d defined in context.

#define SMALL_DOT_SIZE (lmax3d/60.0)

// Convenient scaling macro (assumes r_factor defined in context):

#define _R_(i)  ((int)( ((float)i)*r_factor + 0.5 ))

//-------------------------------------------------------------------------
