#include "sdyn3.h"
#include "xstarplot.h"
#include "dyn_util.h"

//-------------------------------------------------------------------------

// These shouldn't really be global, but keep them this way for now...

unsigned long  win, dia, instr, colwin;
unsigned long c_energy[10], c_index[N_COLORS+1];
int      init_status = 0;

// The plotting box is determined completely by origin and lmax3d.
// Its projection onto the viewing plane is xmin, xmax, ymin, ymax

float  xmin, xmax, ymin, ymax, base_point_size, lmax3d, origin[3];

int   win_size;
int   point_scale_mode = 0;
int   kx = 1, ky = 2, kproj = 3;
int   origin_star = -1;
int   nodes = 0, links = 0, root = 0;
float theta = 0.33, costheta = cos(0.33), sintheta = sin(0.33), dtheta = 0.03;
float phi = 0.33, cosphi = cos(0.33), sinphi = sin(0.33);
real  local_offset[3];
float delay_time;
float r_factor = 1.0;

char  graph3d = 1, track = 0, cenergy = 0; // Convenient to distinguish these.

char   temp_buffer[255];	// Convenient to share this temporary space.

//-------------------------------------------------------------------------

local void check_for_input(unsigned long win, sdyn3* b,
			   bool b_flag, bool f_flag)
{
    char key, string[20], shift, control;

    while (lux_next_keypress(win, &key, string, &shift, &control))
	cerr << "key = " << key << endl;
}

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

void xstarplot(sdyn3* b, float scale, int k, int d, float lmax,
	       int point_mode, float rel_point_size,
	       float D, int ce,
	       bool b_flag, bool f_flag, bool t_flag,
	       int init_flag)
{
    r_factor = scale;
    int xorigin, yorigin;

    initialize_graphics(r_factor, b_flag,
			win, instr, colwin,
			c_index, c_energy,
			win_size, xorigin, yorigin);

    lux_clear_window(win);
    lux_update_fg(win);

    while (1) check_for_input(win, b, b_flag, f_flag);
}

main(int argc, char** argv)
{
    sdyn3 *b;

    int   k = 3;
    int   d = 2;
    float D = 0.0;
    int   cenergy = 0;
    float lmax = -1.0;
    int   point_mode = 0;
    float rel_point_size = -1.0;
    float scale = 1.0;

    bool  b_flag = TRUE;
    bool  f_flag = TRUE;
    bool  o_flag = FALSE;
    bool  t_flag = FALSE;

    xstarplot(b, scale, k, d, lmax,
	      point_mode, rel_point_size, D, cenergy,
	      b_flag, f_flag, t_flag, 0);
}
