
//// xstarplot2:  plot an N-body system in an X environment.
////              This version uses the NEW interpolation routines.
////
//// Options:    -a    specify viewing axis [1/2/3 = x/y/z]
////             -b    allow backward steps [no]
////             -d    dimensionality of plot [2]
////             -D    time interval between frames [0.015625 = 1/64]
////             -e    color by energy [no]
////             -E    toggle display of stars flagged as escapers [true]
////             -f    solid-color stars [true]
////             -F    input file [run.out]
////             -k    maximum size of CM indicator is k/N [2]
////             -l    show links and nodes [no]
////             -L    specify limits [take from initial data]
////             -m    highlight multiple systems [no]
////             -o    pipe cin to cout [no]
////             -p    point scale mode [0]
////             -P    specify point size (relative to axis/30) [1, min = 0.5]
////             -r    "reverse video" (black background) [true]
////             -s    window size, relative to "standard" 400 pixels [1]
////             -t    show tree structure [no]
////             -u    highlight unperturbed systems [no]
////             -v    verbose mode (worldbundle details) [0]
////
////
//// Note:  There appears to be a problem with the X interface on some
////        systems that causes xstarplot to crash if too many run-time
////        user commands are received in a short space of time (e.g.
////        try holding down the ">" key in 3D mode).  Presently, the
////        only fix is to avoid excessive keystrokes in the graphics
////        window at run time...

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
//	    	  Sep 1995   Extended to sdyn3 and nonzero radii, SMcM
//	    	  Aug 1996   Added forward/backward step modes, SMcM
//	    	  Aug 1996   Reworked for hdyn and nonzero radii, SMcM
//                Jul 1998   New end-of-data behavior; more robust
//                           X interface, SMcM
//.............................................................................
//  non-local functions: 
//    xstarplot22
//.............................................................................
//  uses:
//    liblux.a (X interface)
//    The X graphics library
//.............................................................................

// NOTE: On some systems, it appears that overflow can occur, causing the
//	 win.lnx and win.lny entries to be overwritten and leading to
//	 spurious "log10" errors in draw.c (in lux_set_item, specifically).
//	 This isn't fatal, but it should be fixed someday...

#include "worldline.h"
#include "xstarplot.h"
#include "dyn_util.h"

// Use dyn if we only want to deal with dynamical quantities.
// Use pdyn if we want to see physical quantities too.

#define DYN pdyn
#define DYNPTR pdynptr

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
int   nodes = 0, links = 0, root = 0, multiples = 0, unperturbed = 0;
float theta = 0.33, costheta = cos(0.33), sintheta = sin(0.33), dtheta = 0.03;
float phi = 0.33, cosphi = cos(0.33), sinphi = sin(0.33);
real  local_offset[3];

#define MIN_DELAY  1.0
float delay_time = MIN_DELAY;

float r_factor = 1.0;

char  graph3d = 1, track = 0, cenergy = 0; // convenient to distinguish these

char   temp_buffer[255];	// convenient to share this temporary space

//-------------------------------------------------------------------------

// base_point_size defines the size of the default "dot" size on the display.
// It is used (with scaling) for stars and (without scaling) for nodes.

local void set_base_point_size(DYN* b, float rel_point_size)
{
    base_point_size = 2 * SMALL_DOT_SIZE;	// default point size

    if (point_scale_mode == 1) {

	// Scale so the average initial mass has a reasonable size.

	real m_tot = 0;
	real n_tot = 0;
	for_all_leaves(DYN, b, bi) {
	    n_tot++;
	    m_tot += bi->get_mass();
	}

	if (n_tot > 0) base_point_size /= sqrt(m_tot/n_tot);	// NB sqrt below

    } else if (point_scale_mode == 2 || point_scale_mode == 3) {

	// In the case point_scale_mode = 2, use unscaled (i.e. true)
	// radii if rel_point_size < 0.

	if (point_scale_mode == 2 && rel_point_size < 0) {

	    base_point_size = 1;

	} else {

	    // Scale so the average initial radius has a reasonable size.

	    real r_tot = 0;
	    real n_tot = 0;
	    for_all_leaves(DYN, b, bi) {
		n_tot++;
		r_tot += bi->get_radius();
	    }

	    if (n_tot > 0) {
		if (point_scale_mode == 2)
		    base_point_size /= r_tot/n_tot;
		else
		    base_point_size /= sqrt(r_tot/n_tot);
	    }
	}
    }

    // Optional additional scaling:

    if (rel_point_size > 0.0)  base_point_size *= rel_point_size;
}

local float get_point_size(DYN* bi)
{
    // Scaling by mass is too extreme.  Use sqrt(mass).
    // Also, impose a lower limit on the point size.

    // For pdyns, radius is in solar units.  Better place a limit so
    // supergiants don't just vanish!

    if (point_scale_mode == 1)
	return Starlab::max(SMALL_DOT_SIZE,
		   base_point_size * sqrt(bi->get_mass()));
    else if (point_scale_mode == 2)
	return Starlab::max(SMALL_DOT_SIZE,
		   base_point_size * Starlab::min(20., bi->get_radius()));
    else if (point_scale_mode == 3)
	return Starlab::max(SMALL_DOT_SIZE,
		   base_point_size * Starlab::min(20., sqrt(bi->get_radius())));
    else
	return base_point_size;
}

local void draw_star_point(unsigned long win, float r, float s,
			   float actual_point_size, bool f_flag)

// Plot a point of size actual_point_size at (r,s), filled or unfilled
// according to the value of f_flag.

{
    // May be able to improve appearance of larger filled points by
    // drawing a circle around them...

    if (f_flag) {
	lux_fill_arcf(win, r - actual_point_size/2, s - actual_point_size/2,
		      actual_point_size, actual_point_size, 0.0, 360.0);
	if (actual_point_size < 12*lmax3d/win_size) return;
	actual_point_size *= 0.8;
    }
    lux_draw_arcf(win, r - actual_point_size/2, s - actual_point_size/2, 
		  actual_point_size, actual_point_size, 0.0, 360.0);
}

local void draw_links_2d(DYN* b, float r, float s)
{
    for_all_daughters(DYN, b, bb) {

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

local void draw_links_3d(DYN* b, float r, float s)
{
    for_all_daughters(DYN, b, bb) {

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

local int clean_index(DYN* b)
{
    int clean_name = 0;

    if (b->get_index() > 0) 
	clean_name = b->get_index();
    else if (b->get_name()) {
	int i = 0;
	if (sscanf(b->get_name(), "%d", &i)) clean_name = i;
    }

    return clean_name;
}

static real max_cm = 2;

local void plot_star(DYN *bi, float r, float s,
		     float actual_point_size, int f_flag)
{
    // Plot a point at (r, s) representing star bi.
    // Logic regarding colors and resetting link colors needs
    // reconsideration...

    // Determine the color to use.

    bool temp_flag = f_flag;
    if (bi->get_oldest_daughter()) {

	lux_set_color(win, lux_lookup_color(win, "grey"));
	temp_flag =  1 - f_flag;

    } else {

	if (cenergy) {

	    char c;

	    compute_energies(bi->get_root(), bi, c);
	    lux_set_color(win, c_energy[c]);

	} else if (clean_index(bi) > 0) {

	    // Wrap the color map.

	    int ii = clean_index(bi);
	    while (ii > N_COLORS) ii -= N_COLORS;

	    lux_set_color(win,c_index[ii]);
	}
    }

    if (track)

	lux_draw_pointf(win, r, s);

    else {

	real scale = 1.0;
	real reset_grey = false;

	// Deal with possible highlighted systems.
	//
	// CM color schemes:
	//
	//	green		any unbound binary
	//	light blue	bound perturbed binary
	//	blue		bound unperturbed binary
	//	yellow		triple system
	//	orange		quadruple system
	//	red		quintuple or higher-order system

	DYN *od = bi->get_oldest_daughter();

	if (multiples && od && bi->is_top_level_node()) {

	    if (bi->n_leaves() == 2) {

		// PRL(od->get_kepler());

		if (od->get_kepler() == (kepler*)2) {

		    // Highlight lightly perturbed top-level binaries.

		    scale = 2.0;
		    lux_set_color(win, lux_lookup_color(win, "purple"));
		    reset_grey = true;

		} else if (!od->get_kepler()) {

		    // Highlight perturbed top-level binaries.

		    scale = 2.0;
		    lux_set_color(win, lux_lookup_color(win, "lightblue"));
		    reset_grey = true;
		}

	    } else {

		// Highlight multiples.

		scale = 4.0;
		if (bi->n_leaves() == 3)
		    lux_set_color(win, lux_lookup_color(win, "yellow"));
		else if (bi->n_leaves() == 4)
		    lux_set_color(win, lux_lookup_color(win, "orange"));
		else
		    lux_set_color(win, lux_lookup_color(win, "red"));

		reset_grey = true;
	    }
	}

	if (unperturbed && od
	    && od->get_kepler() && od->get_kepler() != (kepler*)2) {

	    // Highlight all unperturbed binaries.

	    scale = 2.0;
	    lux_set_color(win, lux_lookup_color(win, "blue"));
	    reset_grey = true;
	}

	// Plot the point.

	real psize = scale*actual_point_size;

	if (scale > 1.0) {

	    // Compare binary size to actual_point_size and use the
	    // larger of the two.  Best to use semi-major axis, so
	    // we need *velocity* information in this case.  Also,
	    // don't bother with special treatment of unbound pairs.

	    real dr = abs(od->get_pos()
			   - od->get_younger_sister()->get_pos());
	    real dv2 = square(od->get_vel()
			       - od->get_younger_sister()->get_vel());

	    real sma2 = 1 / (1/dr - 0.5*dv2/bi->get_mass());	// >0 if bound

	    // Change color for unbound binary.

	    if (sma2 <= 0) {
		sma2 = -sma2;
		if (scale == 2.0)
		    lux_set_color(win, lux_lookup_color(win, "green"));
	    }

	    // Place limits on point size.

	    if (sma2 > max_cm) sma2 = Starlab::max(dr, max_cm);

	    psize = Starlab::max(psize, 1.2*sma2);
	    //PRC(bi->format_label()); PRC(dr); PRC(dv2); PRL(sma2);
	}

	draw_star_point(win, r, s, psize, temp_flag);

	if (bi->get_oldest_daughter()
	    && psize > 3*actual_point_size) {
	    lux_set_color(win, lux_lookup_color(win, "grey"));
	    draw_star_point(win, r, s, actual_point_size, temp_flag);
	}

	if (links && reset_grey)
	    lux_set_color(win, lux_lookup_color(win, "grey"));
    }
}

local int plot_all_stars(DYN* b, int f_flag)
{
    int  n_stars = b->n_leaves();
    if (b->get_oldest_daughter() == NULL) n_stars = 0;

    int  n_nodes = count_nodes(b);
    float r, s;

    if (graph3d) {

	for_all_nodes(DYN, b, bi) if (root || bi->get_parent()) {
	    if (nodes
		|| (multiples && bi->is_top_level_node())
		|| (unperturbed && bi->get_oldest_daughter())
		|| (bi->get_oldest_daughter() == NULL)) {
		
		float X = (float)bi->get_pos()[0] - local_offset[0] - origin[0];
		float Y = (float)bi->get_pos()[1] - local_offset[1] - origin[1];
		float Z = (float)bi->get_pos()[2] - local_offset[2] - origin[2];

		float actual_point_size = get_point_size(bi);

		if (   (X > (-lmax3d + actual_point_size))
		    && (X < ( lmax3d - actual_point_size)) 
		    && (Y > (-lmax3d + actual_point_size))
		    && (Y < ( lmax3d - actual_point_size)) 
		    && (Z > (-lmax3d + actual_point_size))
		    && (Z < ( lmax3d - actual_point_size))) {

		    // Should really sort by depth here...

		    project3d(X, Y, Z, r, s,
			      costheta, sintheta,
			      cosphi, sinphi);

		    plot_star(bi, r, s, actual_point_size, f_flag);

		    if (links && bi->get_oldest_daughter())
			draw_links_3d(bi, r, s);
		}
	    }
	}

    } else {

	// Make a list of all nodes, *including* the root node if root = 1.
	// Root will be at the start of the list if it is being displayed.

	DYNPTR* p = new DYNPTR[n_nodes+root];

	int ip = 0;
	for_all_nodes(DYN, b, bi) if (root || bi != b) p[ip++] = bi;
	if (ip != n_nodes+root) {
	    cerr << "plot_all_stars: n_nodes = " << n_nodes+root
		<< " counted " << ip << endl;
	    exit(1);
	}

	// Sort by kproj (note that kproj = 1, 2, or 3 for x, y, or z):

	for (ip = 0; ip < n_nodes+root; ip++)
	    for (int jp = ip+1; jp < n_nodes+root; jp++)
		if (p[jp]->get_pos()[kproj-1] < p[ip]->get_pos()[kproj-1]) {
		    DYNPTR bb = p[jp];
		    p[jp] = p[ip];
		    p[ip] = bb;
		}

	// Plot ordered by depth.

	for (ip = 0; ip < n_nodes+root; ip++) {
	    DYN * bi = p[ip];
	    if ( (root || bi != b)
		&& (nodes
		    || (multiples && bi->is_top_level_node())
		    || (unperturbed && bi->get_oldest_daughter())
		    || (bi->get_oldest_daughter() == NULL))) {

		r = (float)bi->get_pos()[kx] - local_offset[kx];
		s = (float)bi->get_pos()[ky] - local_offset[ky];

		float actual_point_size = get_point_size(bi);

		if (   (r > (xmin + actual_point_size))
		    && (r < (xmax - actual_point_size)) 
		    && (s > (ymin + actual_point_size))
		    && (s < (ymax - actual_point_size)) ) {

		    if (bi->get_index() != -42)
			plot_star(bi, r, s, actual_point_size, f_flag);

		    if (links && bi->get_oldest_daughter())
			draw_links_2d(bi, r, s);
		}
	    }
	}
	delete [] p;
    }
    return n_stars;
}

local int type(int which)
{
    if (which == tracking || which == graph3dim || which == colorenergy)
	return BUTTON_WINDOW;
    else
	return INPUT_WINDOW;
}

local int subtype(int which)
{
    if (which == tracking || which == graph3dim || which == colorenergy)
	return CHECK_BUTTON;
    else
	return NO_TYPE;
}

local void set_temp_buffer(int which)
{
    // Translate ID into a string representing the value.

    char c = -1;
    int i = -1;
    float f = 1.e30;	// Don't use VERY_LARGE_NUMBER here (DEC may object!)

    switch (which) {
	case colorenergy:	c = cenergy; break;
	case tracking:		c = track;   break;
	case graph3dim:		c = graph3d; break;
	case xminimum:		f = xmin; break;
	case xmaximum:		f = xmax; break;
	case yminimum:		f = ymin; break;
	case ymaximum:		f = ymax; break;
	case basepointsize:	f = base_point_size; break;
	case pointscalemode:	i = point_scale_mode; break;
	case lmax3D:		f = 2*lmax3d; break;
	case theta3D:		f = theta; break;
	case phi3D:		f = phi; break;
	case DelayTime:		f = delay_time; break;
	case dtheta3D:		f = dtheta; break;
	case Xorigin:		f = origin[0]; break;
	case Yorigin:		f = origin[0]; break;
	case Zorigin:		f = origin[0]; break;
	case view2D:		i = kproj; break;
	case originstar:	i = origin_star; break;
	default:		cerr << "xstarplot22: diag error...\n";
    }

    if (f < 1.e30)
	sprintf(temp_buffer, "%f", f);
    else if (c == -1)
	sprintf(temp_buffer, " %d ", i);
    else
	temp_buffer[0] = c;
}

local void set_diag_item(int which, char* id, int label_on_left,
			 int x1, int x2, int y)
{
    set_temp_buffer(which);

    if (label_on_left) {	// label to left of (real) input box

	lux_set_item(dia, which, TEXT_WINDOW, NO_TYPE,
		     _R_(x1), _R_(y), strlen(id), id);
	lux_set_item(dia, which, INPUT_WINDOW, NO_TYPE,
		     _R_(x2), _R_(y), 10, temp_buffer);

    } else {			// label to right of (integer) input box/button

	if (type(which) == INPUT_WINDOW) {

	    lux_set_item(dia, which, INPUT_WINDOW, NO_TYPE,
			 _R_(x1), _R_(y), 3, temp_buffer);
	    lux_set_item(dia, which, TEXT_WINDOW, NO_TYPE,
			 _R_(x2), _R_(y), strlen(id), id);

	} else {

	    lux_set_item(dia, which, BUTTON_WINDOW, CHECK_BUTTON,
			 _R_(x1), _R_(y), 1, temp_buffer);
	    lux_set_item(dia, which, TEXT_WINDOW, CHECK_BUTTON,
			 _R_(x2), _R_(y), strlen(id), id);
	}
    }
}

local void initialize_dialog(int xorigin, int yorigin)
{
    // Initialize dialog box material (note: y is measured up from bottom!)

    float save = r_factor;
    if (r_factor > 1) r_factor = 1;

    int xsize = 550;
    int ysize = 550;
    if (r_factor <= 0.6) xsize = (int) (xsize / (r_factor/0.6));
    dia = lux_open_dialog(xorigin, yorigin, _R_(xsize), _R_(ysize));

    // ---------- 3-D origin info across top of box (special case): ----------

     lux_set_item(dia, Origin, TEXT_WINDOW, NO_TYPE,
		 _R_(70), _R_(500), 6, "origin");

    sprintf(temp_buffer, "%f", origin[0]);
    lux_set_item(dia, Xorigin, INPUT_WINDOW, NO_TYPE,
		 _R_(160), _R_(500), 10, temp_buffer);

    sprintf(temp_buffer, "%f", origin[1]);
    lux_set_item(dia, Yorigin, INPUT_WINDOW, NO_TYPE,
		 _R_(280), _R_(500), 10, temp_buffer);

    sprintf(temp_buffer, "%f", origin[2]);
    lux_set_item(dia, Zorigin, INPUT_WINDOW, NO_TYPE,
		 _R_(400), _R_(500), 10, temp_buffer);


    // ---------- Left-hand column of dialog box: ----------

    //		xmin
    //		xmax
    //		ymin
    //		ymax
    //		3-D cube size
    //		projection theta
    //		projection phi
    //		projection dtheta
    //		point scale
    //		delay time

    // Minima and maxima of 2-D display:

    set_diag_item(xminimum, "xmin", 1, 70, 160, 460);
    set_diag_item(xmaximum, "xmax", 1, 70, 160, 420);
    set_diag_item(yminimum, "ymin", 1, 70, 160, 380);
    set_diag_item(ymaximum, "ymax", 1, 70, 160, 340);

    // 3-D cube size:

    set_diag_item(lmax3D, "cube size", 1, 70, 180, 300);

    // 3-D projection direction:

    set_diag_item(theta3D,  "theta",  1, 70, 180, 260);
    set_diag_item(phi3D,    "phi",    1, 70, 180, 220);
    set_diag_item(dtheta3D, "dtheta", 1, 70, 180, 180);

    // Point size:

    set_diag_item(basepointsize, "point scale", 1, 70, 180, 140);

    // Delay time:

    set_diag_item(DelayTime, "delay time", 1, 70, 180, 100);


    // ---------- Right-hand column of dialog box: ----------

    //		origin star
    //		2-D view axis
    //		3-D plot?

    //		point display mode
    //		color by energy?
    //		track particles?

    // Origin star:

    set_diag_item(originstar, "Origin Star", 0, 320, 360, 460);

    // 2-D view axis:

    set_diag_item(view2D, "2-D view axis", 0, 320, 360, 420);

    // 3D graphics:

    set_diag_item(graph3dim, "3-D graph", 0, 320, 350, 380);

    // Point display mode

    set_diag_item(pointscalemode, "point display mode", 0, 320, 360, 300);

    // Color by energy:

    set_diag_item(colorenergy, "color by energy", 0, 320, 350, 260);

    // Track particles:

    set_diag_item(tracking, "track particles", 0, 320, 350, 220);

    // lux_draw_palette(dia);

    // OK/cancel:

    lux_set_item(dia, ok, BUTTON_WINDOW, OK_BUTTON,
		 _R_(100), _R_(25), 9, "OK, CLEAR");
    lux_set_item(dia, ok_keep, BUTTON_WINDOW, OK_KEEP_BUTTON,
		 _R_(250), _R_(25), 8, "OK, KEEP");
    lux_set_item(dia, cancel, BUTTON_WINDOW, CANCEL_BUTTON,
		 _R_(400), _R_(25), 6, "CANCEL");

    r_factor = save;
}

local void make_relative_to_root(DYN* b)
{
    vec root_pos = b->get_pos();
    vec root_vel = b->get_vel();
    for_all_nodes(DYN, b, bi) {
	bi->inc_pos(-root_pos);
	bi->inc_vel(-root_vel);
    }
}

local void update_diag_item(int which)
{
    set_temp_buffer(which);
    lux_update_itemvalue(dia, which, type(which), subtype(which),
			 temp_buffer);
}

local void get_diag_string(int which) {
}

// Overloaded function:

local void read_diag_item(int which, float& f)
{
    lux_get_itemvalue(dia, which, type(which), subtype(which),
		      temp_buffer);
    sscanf(temp_buffer, "%f", &f);
}

local void read_diag_item(int which, int& i)
{
    lux_get_itemvalue(dia, which, type(which), subtype(which),
		      temp_buffer);
    sscanf(temp_buffer, "%d", &i);
}

local void read_diag_item(int which, char& c)
{
    lux_get_itemvalue(dia, which, type(which), subtype(which),
		      temp_buffer);
    c = temp_buffer[0];
}

local void update_from_dialog(bool r_flag)
{
    // Read all data from the dialog box.  Note that it is
    // important that the box be maintained correctly.
    // Any errors in the data displayed in the box will be
    // propogated to the "real" data by this procedure.

    read_diag_item(view2D, kproj);
    if (kproj == 1)
	kx = 1, ky = 2;
    else if (kproj == 2)
	kx = 2, ky = 0;
    else
	kx = 0, ky = 1;

    read_diag_item(originstar, origin_star);

    read_diag_item(Xorigin, origin[0]);
    read_diag_item(Yorigin, origin[1]);
    read_diag_item(Zorigin, origin[2]);

    read_diag_item(lmax3D, lmax3d);
    lmax3d /= 2;			// display twice lmax3d

    // ------------------------------------------------------

    // NOTE: xmin, xmax, ymin, ymax are IRRELEVANT, since they
    // will be forced to be consistent with origin and lmax3d.

    read_diag_item(xminimum, xmin);
    read_diag_item(xmaximum, xmax);
    read_diag_item(yminimum, ymin);
    read_diag_item(ymaximum, ymax);

    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

    // Make sure entries in the dialog box are consistent:

    update_diag_item(xminimum);
    update_diag_item(xmaximum);
    update_diag_item(yminimum);
    update_diag_item(ymaximum);

    // ------------------------------------------------------

    read_diag_item(basepointsize,base_point_size);

    read_diag_item(theta3D, theta);
    costheta = cos(theta);
    sintheta = sin(theta);

    read_diag_item(phi3D, phi);
    cosphi = cos(phi);
    sinphi = sin(phi);

    read_diag_item(dtheta3D, dtheta);

    read_diag_item(DelayTime, delay_time);
    if (delay_time < MIN_DELAY) delay_time = MIN_DELAY;

    read_diag_item(colorenergy, cenergy);
    show_color_scheme(colwin, c_energy, c_index,
		      r_factor, cenergy, r_flag, 1);

    read_diag_item(tracking, track);
    read_diag_item(graph3dim, graph3d);

    lux_clear_window(win);

    if (graph3d) {
	lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
		            -FAC3D*lmax3d, FAC3D*lmax3d);
	draw3d_axis(win, lmax3d, costheta, sintheta,
		    cosphi, sinphi);
    } else {
	lux_setup_axis(win, xmin, xmax, ymin, ymax);
	draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);
    }

    read_diag_item(pointscalemode, point_scale_mode);
}

local void show_static_rotation(DYN* b, bool f_flag)
{
    // Handle static rotation using lux_check_keypress:

    show_instructions(instr, r_factor,
		      "Pause...\n  r: rotate, (p,c): continue                ",
		      1);

    char rotate = 0;
    do {
	if (lux_check_keypress(win,'r')) {
	    show_instructions(instr, r_factor,
		      "Rotate... p: pause, c: continue, c: exit rotation     ",
			      0);
	    rotate = 1;

	    do {
		theta += dtheta;
		if (theta < 0.0)
		    theta += 2*PI;
		else if (theta > PI*2)
		    theta -= 2*PI;
		costheta = cos(theta);
		sintheta = sin(theta);
		update_diag_item(theta3D);

		lux_clear_current_region(win);
		draw3d_axis(win, lmax3d, costheta, sintheta,
			    cosphi, sinphi);

		for_all_leaves(DYN, b, bi) {

		    float X = (float)bi->get_pos()[0] - origin[0];
		    float Y = (float)bi->get_pos()[1] - origin[1];
		    float Z = (float)bi->get_pos()[2] - origin[2];

		    float actual_point_size = get_point_size(bi);

		    if (   (X > (-lmax3d + actual_point_size))
			&& (X < ( lmax3d - actual_point_size)) 
			&& (Y > (-lmax3d + actual_point_size))
			&& (Y < ( lmax3d - actual_point_size)) 
			&& (Z > (-lmax3d + actual_point_size))
			&& (Z < ( lmax3d - actual_point_size))){

			float r, s;
			project3d(X, Y, Z, r, s,
				  costheta, sintheta,
				  cosphi, sinphi);

			if (cenergy) {
			    char c;
			    compute_energies(b, bi, c);
			    lux_set_color(win,c_energy[c]);
			} else if (clean_index(bi)>0) {
			    if (clean_index(bi) <= N_COLORS)
				lux_set_color(win,
					      c_index[clean_index(bi)]);
			    else
				lux_set_color(win,c_index[1]);
			}

			if (track) 
			    lux_draw_pointf(win,r,s);
			else
			    if (f_flag)
				lux_fill_arcf(win,
					      r - actual_point_size/2,
					      s - actual_point_size/2, 
					      actual_point_size,
					      actual_point_size,
					      0.0, 360.0);
			    else
				lux_draw_arcf(win,
					      r - actual_point_size/2,
					      s - actual_point_size/2, 
					      actual_point_size,
					      actual_point_size,
					      0.0, 360.0);
		    }
		}
		update_with_delay(win, Starlab::max(MIN_DELAY, delay_time));

		if (lux_check_keypress(win,'p')) 
		    while(!lux_check_keypress(win,'c'));

	    } while(!lux_check_keypress(win,'c'));

	    // Flush "c"s from input stream:

	    while(lux_check_keypress(win,'c'));
	}

    } while(!lux_check_keypress(win,'p')
	    && !lux_check_keypress(win,'c') && !rotate);

    rotate = 0;
    lux_set_color(win,c_energy[default_color]);
}

local char check_for_input(unsigned long win, DYN* b,
			   real &t, real &dt, bool &E_flag,
			   bool r_flag, bool f_flag, bool eod)
{
    // Check for (and immediately act upon) interactive input.
    // Loop through the input buffer until no more events remain.

    // lux_next_keypress gets and discards the next key in the input stream.

    char key, string[20], shift, control;
    bool key_pressed = false;

    while (lux_next_keypress(win, &key, string, &shift, &control)) {

        // A "defined" key has been pressed.  See if it is one we want.

	key_pressed = true;

        if (key == 0			// key = 0 for non-ASCII character
    			 		// e.g. Up, Down, Home, etc.--see win.c
	    || key == 'h'
	    || key == 'a') {

	    // "Shift" functions:

	    int change = 0;

	    if (key == 'h')		// keyboard lookalikes for
		change = 4;		// function keys
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
		change = kproj;
		origin[kproj-1] += lmax3d/2;
	    } else if (strcmp(string, "PgDn") == 0) {
		change = kproj;
		origin[kproj-1] -= lmax3d/2;
	    } else if (strcmp(string, "Home") == 0)
		change = 4;
	    else if (strcmp(string, "R11") == 0)
		change = 5;

	    if (change) {

		if (change == 4) {
		    origin[kx] = origin[ky] = origin[kproj-1] = 0;
		} else if (change == 5) {

		    real save = lmax3d;
		    lmax3d = 0;

		    // Use all dimensions in redetermining the scale!

		    for_all_leaves(DYN, b, bi)
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

		    origin[kx] = origin[ky] = origin[kproj-1] = 0;		
		    base_point_size *= lmax3d/save;
		}

		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);

		// Update all entries in the dialog box, as needed:

		update_diag_item(xminimum);
		update_diag_item(xmaximum);
		update_diag_item(yminimum);
		update_diag_item(ymaximum);

		if (change == 1 || change >= 4) update_diag_item(Xorigin);
		if (change == 2 || change >= 4) update_diag_item(Yorigin);
		if (change == 3 || change >= 4) update_diag_item(Zorigin);

		if (change == 5) {
		    update_diag_item(basepointsize);
		    update_diag_item(lmax3D);
		}

		// No redrawing of axes is needed for 3D graphs for change != 5.
		
		if (!graph3d) {
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);
		} else if (change == 5) {
		    lux_clear_window(win);
		    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		}
	    }

	} else {			// key is the ASCII code

	    // Check for "quit" first.

	    if (key == 'q') {
		show_instructions(instr, r_factor,
		    "\n   Quitting...                                        ",
				  1);
		lux_exit();		// note that this is now quick_exit
	    }

	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	    if (graph3d && !track) {
		int button;

		// 3-D operations:

		if (key== 'R') {

		    show_static_rotation(b, f_flag);

		} else {

		    button = lux_check_buttonpress(win);

		    if (key == '^') {
			phi = phi - dtheta;
			if (phi < -PI) phi = phi + 2.0*PI;
			cosphi = cos(phi);
			sinphi = sin(phi);
			update_diag_item(phi3D);
		    } else if (key == 'V') {
			phi = phi + dtheta;
			if (phi > PI) phi = phi - 2.0*PI;
			cosphi = cos(phi);
			sinphi = sin(phi);
			update_diag_item(phi3D);
		    } else if (key == '<') {
			theta = theta + dtheta;
			if (theta > 2*PI) theta -= 2*PI;
			costheta = cos(theta);
			sintheta = sin(theta);
			update_diag_item(theta3D);
		    } else if (key == '>') {
			theta = theta - dtheta;
			if (theta < 0.0) theta += 2*PI;
			costheta = cos(theta);
			sintheta = sin(theta);
			update_diag_item(theta3D);
		    }
		}
	    }

	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	    if (key == 'p' || key == 'P') {	// modify base point size

		if (key == 'p') {
		    if (base_point_size > 4*lmax3d/win_size)
			base_point_size /= PFAC;
		} else
		    base_point_size *= PFAC;

		update_diag_item(basepointsize);

	    } else if (key == 'm') {		// toggle multiple display

		multiples = 1 - multiples;
		if (multiples == 1) nodes = 1;

	    } else if (key == 'n') {		// toggle tree display

		nodes = 1 - nodes;
		if (nodes == 0) links = 0;

	    } else if (key == 'l') {		// toggle links/nodes

		links = 1 - links;
		if (links == 1) nodes = 1;

	    } else if (key == 'r') {		// show root node

		root = 1 - root;

	    } else if (key == 't') {		// toggle tracking

		track = 1 - track;
		update_diag_item(tracking);

	    } else if (key == 'u') {		// toggle unperturbed display

		unperturbed = 1 - unperturbed;

	    } else if (key == '3' && !graph3d) {    // enter 3D mode

		graph3d = 1;
		update_diag_item(graph3dim);
		lux_clear_window(win);
		lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
			       -FAC3D*lmax3d, FAC3D*lmax3d);
		draw3d_axis(win, lmax3d, costheta, sintheta, cosphi, sinphi);

		// Refresh the color window.

//		show_color_scheme(colwin, c_energy, c_index,
//				  r_factor, cenergy, r_flag, 1);

	    } else if (graph3d &&
		       (key == '2' || key == 'x'
			  || key == 'y' || key == 'k')) {   // enter 2D mode

		graph3d = 0;
		update_diag_item(graph3dim);
		lux_clear_window(win);
		lux_setup_axis(win, xmin, xmax, ymin, ymax);
		draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		// Refresh the color window.

//		show_color_scheme(colwin, c_energy, c_index,
//				  r_factor, cenergy, r_flag, 1);

	    } else if (key == 'E') {		// toggle escaper display

		E_flag = !E_flag;

	    } else if (key == 'e') {		// color by energy

		cenergy = !cenergy;
		update_diag_item(colorenergy);
		show_color_scheme(colwin, c_energy, c_index,
				  r_factor, cenergy, r_flag, 1);

	    } else if (key == '*') {    	// increase dt

//#define DTFAC	1.77827941003892275873		// fourth root of 10
#define DTFAC	1.41421356237309514547		// square root of 2

		dt *= DTFAC;

	    } else if (key == '/') {		// decrease dt

		dt /= DTFAC;

	    } else if (key == '\\') {		// reverse dt

		dt = -dt;

	    } else if (key == 'j') {

		cerr << "jump to time: " << flush;
		cin >> t;
		t -= dt;
		return 'j';

	    } else if (key == 'J') {

		cerr << "time, dt: " << flush;
		cin >> t >> dt;
		t -= dt;
		return 'J';

	    } else if (key == '>') {

		max_cm *= 2;

	    } else if (key == '<') {

		max_cm /= 2;

	    } else if (key == '+' || key == '=') {    // increase delay time

		delay_time *= 2;
		update_diag_item(DelayTime);

	    } else if (key == '-') {		// decrease delay time

		delay_time /= 2;
		if (delay_time < MIN_DELAY) delay_time = MIN_DELAY;
		update_diag_item(DelayTime);

	    } else if (key == '0') {		// set delay time = "0"

		delay_time = MIN_DELAY;
		update_diag_item(DelayTime);

		// Flush the queue:

		while(lux_check_keypress(win,'+')
		      || lux_check_keypress(win,'-'));

	    } else if (key == 'R') {		// color by energy

		show_color_scheme(colwin, c_energy, c_index,
				  r_factor, cenergy, r_flag, 1);

	    } else if (key == 'z') {		// zoom in
						// (fixed axes and origin)

		lmax3d /= ZOOM;
		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		base_point_size /= ZOOM;

		if (graph3d) {
		    lux_clear_window(win);
		    lux_setup_axis(win,-FAC3D*lmax3d,FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		} else {
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);
		}

		update_diag_item(basepointsize);
		update_diag_item(lmax3D);
		update_diag_item(xminimum);
		update_diag_item(xmaximum);
		update_diag_item(yminimum);
		update_diag_item(ymaximum);

	    } else if (key == 'Z') {		// zoom out

		lmax3d *= ZOOM;
		set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		base_point_size *= ZOOM;

		if (graph3d) {
		    lux_clear_window(win);
		    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
				   -FAC3D*lmax3d, FAC3D*lmax3d);
		    draw3d_axis(win, lmax3d, costheta, sintheta,
				cosphi, sinphi);
		} else {
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);
		}

		update_diag_item(basepointsize);
		update_diag_item(lmax3D);
		update_diag_item(xminimum);
		update_diag_item(xmaximum);
		update_diag_item(yminimum);
		update_diag_item(ymaximum);

	    } else if (key == 'k') {		// change projection axis to z
						// (fixed origin, scale)
		if (graph3d) {

		    theta = -PI/2;
		    phi = PI/2.0;
		    costheta = cos(theta);
		    sintheta = sin(theta);
		    cosphi = cos(phi);
		    sinphi = sin(phi);

		    update_diag_item(theta3D);
		    update_diag_item(phi3D);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );

		} else {

		    kproj = 3;
		    kx = 0;
		    ky = 1;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		    update_diag_item(view2D);

		}

	    } else if (key == 'y') {		// change projection axis to y

		if (graph3d) {

		    theta = -PI/2.0;
		    phi = 0.0;
		    costheta = cos(theta);
		    sintheta = sin(theta);
		    cosphi = cos(phi);
		    sinphi = sin(phi);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );

		    update_diag_item(theta3D);
		    update_diag_item(phi3D);

		} else {

		    kproj = 2;
		    kx = 2;
		    ky = 0;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		    update_diag_item(view2D);

		}

	    } else if (key == 'x') {		// change projection axis to x

		if (graph3d) {
		    theta = 0.0;
		    phi = 0.0;
		    costheta = cos(theta);
		    sintheta = sin(theta);
		    cosphi = cos(phi);
		    sinphi = sin(phi);

		    while (lux_check_keypress(win,'<')
			   || lux_check_keypress(win,'>')
			   || lux_check_keypress(win,'^')
			   || lux_check_keypress(win,'V') );

		    update_diag_item(theta3D);
		    update_diag_item(phi3D);

		} else {

		    kproj = 1;
		    kx = 1;
		    ky = 2;
		    set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
		    lux_clear_window(win);
		    lux_setup_axis(win, xmin, xmax, ymin, ymax);
		    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		    update_diag_item(view2D);

		}
	    }

	    // Select new origin: 'o' ==> point in space, 'O' ==> star

	    if (key == 'o' && !graph3d) {

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
		draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		origin_star = -1;

		update_diag_item(Xorigin);
		update_diag_item(Yorigin);
		update_diag_item(Zorigin);
		update_diag_item(xminimum);
		update_diag_item(xmaximum);
		update_diag_item(yminimum);
		update_diag_item(ymaximum);
		update_diag_item(originstar);

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
		draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);

		update_diag_item(Xorigin);
		update_diag_item(Yorigin);
		update_diag_item(Zorigin);
		update_diag_item(xminimum);
		update_diag_item(xmaximum);
		update_diag_item(yminimum);
		update_diag_item(ymaximum);
		update_diag_item(originstar);

	    } else if (key == 'b' || key == 'B') {

		return 'b';

	    } else if (key == 'c' || key == 'C') {

		return 'c';

	    } else if (key == 'f' || key == 'F' || key == 's' || key == 'S') {

		return 'f';

	    }

	    // -------------------------------------------------------------

	    while(lux_check_keypress(win,'r'));  // clean up any replay garbage

	    // Idle mode and dialog box input:

	    int return_value = 0;

	    if (key == 'i') {

		show_instructions(instr, r_factor,
		  " Idle...\n  d: dialog, r: replay, c: continue, q: quit",1);
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
		lux_show_dialog(dia);
		return_value = lux_getevent();
	    }
	    
	    // Update all values if OK button was clicked (clicking
	    // OK forces a return value of 3).
	    
	    if (return_value == 3) update_from_dialog(r_flag);

	}
    }

    if (key_pressed) {

	// Hitting any key will cause the instruction and color
	// windows to be updated immediately.

	show_main_instructions(instr, r_factor, graph3d, 1, eod,
			       nodes, links, multiples, unperturbed, root);
	show_color_scheme(colwin, c_energy, c_index,
			  r_factor, cenergy, r_flag, 1);

	return ' ';
    }

    return 0;
}

/*-----------------------------------------------------------------------------
 *  xstarplot22.C  --  project the positions of all particles onto the screen
 *                    input:   pn: a pointer to a nbody system,
 *		                k: the number of the coordinate axis along which
 *		                   the projection is directed, with {k,x,y}
 *		                   having a right-handed orientation,
 *                            xmax: displayed x-axis spans [-xmax, xmax]
 *                            ymax: displayed y-axis spans [-ymax, ymax]
 *-----------------------------------------------------------------------------
 */

#define INSTR_FREQ_MAX 32
#define INSTR_FREQ_FAC 64
static int instr_freq = INSTR_FREQ_MAX;
static int count_instr = 0;

void xstarplot22(DYN* b, float scale, int k, int d, float lmax,
		int point_mode, float rel_point_size,
		real &t, real &dt, float D, int ce,
		bool &E_flag, bool r_flag, bool f_flag, bool t_flag,
		int init_flag, int& step_mode, bool eod)
{
    if (!b) return;

    r_factor = scale;

    if (!eod && init_flag == 0) {

	int xorigin, yorigin;

	if (t_flag) nodes = links = 1;

	// Open windows for plotting and instructions, set up color
	// schemes, etc.

	if (init_status == 0) initialize_graphics(r_factor, r_flag,
						  win, instr, colwin,
						  c_index, c_energy,
						  win_size, xorigin, yorigin);
	lux_clear_window(win);
	lux_update_fg(win);
	lux_clear_window(colwin);
	lux_update_fg(colwin);
	lux_clear_window(instr);
	lux_update_fg(instr);

	// Determine which axes to plot:

	kproj = k;
	switch (kproj) {
	    case 1:	kx = 1; ky = 2; graph3d = 0; break;
	    case 2:	kx = 2; ky = 0; graph3d = 0; break;
	    case 3:	kx = 0; ky = 1; graph3d = 0; break;
	    default: 	cerr << "xstarplot22: k = " << k
		             << ": illegal value; choose from {1, 2, 3}"
			     << endl;
			exit(0);
	}

	switch(d) {
	    case 2:     graph3d = 0;  break;
	    case 3:     graph3d = 1;  break;
	    default:    cerr << "xstarplot22: d = " << d
		             << " illegal value; choose from {2, 3)"
			     << endl;
			exit(0);
	}

	point_scale_mode = point_mode;

	cenergy = ce?1:0;
	show_color_scheme(colwin, c_energy, c_index,
			  r_factor, cenergy, r_flag, 1);

	delay_time = D;

	if (lmax > 0.0)

	    lmax3d = lmax;

	else {

	    lmax3d = 0;

	    // Use all dimensions in determining the initial scale!

	    for_all_leaves(DYN, b, bi)
		for (int kk = 0; kk < 3; kk++)
		    lmax3d = Starlab::max(lmax3d, abs(bi->get_pos()[kk]));

	    // Round lmax3d up to something reasonable:

	    lmax3d *= 1.2;

	    // Hmmm...  Maybe we should exclude some outliers.

	    lmax3d /= 2;

	    if (lmax3d <= 0) lmax3d = 1;

	    real m_log = log10(lmax3d);
	    int i_log = (int) m_log;
	    if (m_log - i_log > 0.699) i_log++;
	    real scale = pow(10.0, i_log);
	    lmax3d = ((int) (lmax3d/scale) + 1) * scale;
	}

	for (int ko = 0; ko < 3; ko++)
	    origin[ko] = 0;

	set_limits(origin, lmax3d, kx, xmin, xmax, ky, ymin, ymax);
	set_base_point_size(b, rel_point_size);
	
	if (graph3d) {
	    lux_setup_axis(win, -FAC3D*lmax3d, FAC3D*lmax3d, 
			   -FAC3D*lmax3d, FAC3D*lmax3d);      
	    draw3d_axis(win, lmax3d, costheta, sintheta, cosphi, sinphi);
	} else {
	    lux_setup_axis(win, xmin, xmax, ymin, ymax);
	    draw2d_axis(win, xmin, xmax, ymin, ymax, kproj);
	}

	// Complete the initialization with the dialog box.

	if (init_status == 0)
	    initialize_dialog(xorigin, yorigin);

	init_status = 1;
	    
    }

    //------------------------------------------------------------------------

    // End of initialization or start of next plot.
    // Delete current points, and possibly redraw axes and instructions.

    if (track == 0) {
	lux_clear_current_region(win);
	if (graph3d) draw3d_axis(win, lmax3d, costheta, sintheta,
				 cosphi, sinphi);
    }

    if (!eod) {

	// Showing these instructions too frequently causes malloc problems
	// in lux_draw_image_string (for unknown reasons; possibly this is
	// a timing problem of some sort with the X display).
	//
	// Work around the problem by reducing the display frequency.
	// See also gfx_util for further limitations.

	instr_freq = INSTR_FREQ_FAC / b->n_leaves();
	if (instr_freq > INSTR_FREQ_MAX) instr_freq = INSTR_FREQ_MAX;
	if (instr_freq < 1) instr_freq = 1;

	if (count_instr++%instr_freq == 0)
	    show_main_instructions(instr, r_factor, graph3d, 1, eod,
				   nodes, links, multiples, unperturbed, root);
    }

    //------------------------------------------------------------------------

    // Always express all positions and velocities relative to the root node.
    // (There is no requirement that the root node be at rest at the origin...)

    make_relative_to_root(b);

    // Determine local offset, if any.

    if (origin_star > 0) {
	for_all_leaves(DYN, b, bi)
	    if (clean_index(bi) == origin_star) {
		for (int kk = 0; kk < 3; kk++)
		    local_offset[kk] = bi->get_pos()[kk];
		break;
	    }
    } else
	for (int kk = 0; kk < 3; kk++)
	    local_offset[kk] = 0;

    // Plot the data points (sorted by kproj-component in the 2D case):

    int n_stars = plot_all_stars(b, f_flag);

    // Update headers.

    lux_set_color(win, c_energy[default_color]);

    if (graph3d) {

	//sprintf(temp_buffer, "   N = %d (snapshot #%5d) max=%5.3f   ",
	//	n_stars, init_flag + 1, lmax3d);
	sprintf(temp_buffer, "   N = %d (t = %.6f) max=%5.3f   ",
		n_stars, b->get_real_system_time(), lmax3d);
	lux_draw_image_string(win, 0.0, lmax3d*FAC3D, 0.5, temp_buffer, 0);

    } else {

	//sprintf(temp_buffer, "   N = %d  (snapshot #%5d)   ",
	//	n_stars, init_flag + 1);
	sprintf(temp_buffer, "   N = %d  (t = %.6f)   ",
		n_stars, b->get_real_system_time());
	lux_draw_image_string(win, (xmin+xmax)/2.0, ymax, 0.5, temp_buffer, 0);

    }

    // Update the display.

    update_with_delay(win, Starlab::max(MIN_DELAY, delay_time));

    // Deal with user input.  In step mode, wait for
    // 'b', 'c', or 's' before continuing.

    while (true) {

	char c = check_for_input(win, b, t, dt, E_flag, r_flag, f_flag, eod);

 	// Special treatment required for forward/backward step modes.

	if (c == 'b' || c == 'f') {		// Enter/continue step mode

	    step_mode = (c == 'b' ? -1 : 1);
	    return;

	} else if (c == 'c') {			// Cancel step mode

	    if (eod)
		step_mode = 1;
	    else
		step_mode = 0;

	    return;

	} else if (c == 'j' || c == 'J') {

	    step_mode = 1;
	    return;
	}

	// In eod mode, hitting 'b' or 'f' will go to a new frame.
	// Hitting any other key will cause the current frame to be
	// redisplayed (with any specified changes made).

	if (eod && c != 0) {
	    step_mode = 0;
	    return;
	}
	if (!eod && step_mode == 0) return;

	// Wait for 10 milliseconds before retrying.

	lux_pause(10000);
    }

}

//======================================================================

local void print_worldlines(worldbundleptr wh[], int nh)
{
    for (int ih = 0; ih < nh; ih++) {
	cerr << "worldbundle " << ih << ":" << endl;
	wh[ih]->dump(4);
    }
}

local void print_escapes(worldbundleptr wh[], int nh)
{
    DYN *b = create_interpolated_tree2(wh[0], 0);
    for_all_daughters(DYN, b, bb) {

	PRL(bb->format_label());

	// Follow the complete worldline set (all worldbundles) of bb.

	for (int ih = 0; ih < nh; ih++) {

	    worldline *w = wh[ih]->find_worldline(bb);
	    if (w->get_start_esc_flag() || w->get_end_esc_flag()) {
		PRC(ih); PRL(w->get_t_esc());
	    }
	}
    }
}

/*-----------------------------------------------------------------------------
 *  main  --  driver to use  xstarplot22()  as a tool. 
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
    // Establish some defaults:

    int   k = 3;		// note x = 1, y = 2, z = 3 here!
    int   d = 2;
    float D = 16*MIN_DELAY;
    real dt = 0.015625;		// powers of 2 are preferred, but not essential
    int   cenergy = 0;
    float lmax = -1.0;
    int   point_mode = 0;
    float rel_point_size = 0.5;
    float scale = 1.0;

    char infile[128];
    strcpy(infile, "run.out");

    bool f_flag = TRUE;		// if TRUE, the star is solid
    bool o_flag = FALSE;	// if TRUE, output stdin to stdout
    bool r_flag = TRUE;		// if TRUE, the background is black
    bool t_flag = FALSE;	// if TRUE, show tree structure
    bool E_flag = true;		// if true, show escapers

    int verbose = 0;

    check_help();

    extern char *poptarg;
    char* params = "a:bd:D:eEfF:lL:mop:P:rs:tuv.";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'a': k = atoi(poptarg);	// projection axis [z]
		      break;
	    case 'd': d = atoi(poptarg);	// number of dimensions [2]
		      break;
	    case 'D': dt = atof(poptarg);	// delay between frames [0.01]
		      if (dt < 0)		// usual convention
			  dt = pow(2.0, dt);
		      break;
	    case 'e': cenergy = 1;		// color by energy [no]
		      break;
	    case 'E': E_flag = !E_flag;
	    	      break;
	    case 'f': f_flag = FALSE;		// fill star points [yes]
		      break;
	    case 'F': strcpy(infile, poptarg);
	    	      break;
	    case 'k': max_cm = atoi(poptarg);	// max CM size (/N) [2]
		      break;
	    case 'L': lmax = atof(poptarg);	// axes +/- lmax [use data]
		      break;
	    case 'l': links = 1 - links;
	    	      nodes = links;
		      break;
	    case 'm': multiples = 1 - multiples;
	    	      break;
	    case 'o': o_flag = TRUE;		// output stdin to stdout [no]
		      break;
	    case 'p': point_mode = atoi(poptarg);	// point scale mode [0]
		      break;
	    case 'P': rel_point_size = atof(poptarg);   // point scale factor
		      if (rel_point_size < 0.5)
			  rel_point_size = 0.5;
		      break;
	    case 'r': r_flag = FALSE;		// black background [yes]
		      break;
	    case 's': scale = atof(poptarg);	// rescale display [1.0]
		      break;
	    case 't': t_flag = TRUE;		// show tree links [no]
		      break;
	    case 'u': unperturbed = 1 - unperturbed;
	    	      break;
	    case 'v': if (poptarg)
			  verbose = atoi(poptarg);
	    	      else
			  verbose = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], params);
	              get_help();
		      exit(0);
	}

    ifstream s(infile);
    if (!s) {
	cerr << "Data file " << infile << " not found." << endl;
	exit(1);
    }

    // Read in an array of worldbundles (a "worldhistory", when this
    // becomes the next level of structure in the world hierarchy).
    // Eventually, the display can start as soon as the first
    // worldbundle is ready, and the remainder can be read in
    // asynchronously.

    worldbundleptr wb, wh[1024];

    int nh = 0;
    while (nh < 1024 && (wb = read_bundle(s, verbose))) wh[nh++] = wb;

    // PRL(mass_scale_factor());

    cerr << endl << "statistics on " << nh << " worldbundle";
    if (nh != 1) cerr << "s";
    cerr << ":" << endl;

    int nwtot = 0, nstot = 0, netot = 0;
    for (int ih = 0; ih < nh; ih++) {
	wb = wh[ih];
	real t = wb->get_t_min();
	int nw = wb->get_nw(), ns = count_segments(wb), ne = count_events(wb);
	cerr << "worldbundle " << ih << ": "
	     << nw << " worldlines, "
	     << ns << " segments, "
	     << ne << " events, t = "
	     << wb->get_t_min() << " to " << wb->get_t_max()
	     << endl;
	nwtot += nw;
	nstot += ns;
	netot += ne;
    }
    cerr << "totals: " << nwtot << " worldlines, "
	 << nstot << " segments, " << netot << " events"
	 << endl << endl;

    // print_worldlines(wh, nh);
    // print_escapes(wh, nh);

    set_center(wh, nh, 1, true);
    PRC(get_center()); PRL(get_center_id());

    // Now display the data.

    bool eod = false;
    int step_mode = 0;
    int init = 0;

    int ih = 0;
    wb = wh[ih];
    real t = wb->get_t_min();

    bool max_cm_set = false;

    while (1) {

	DYN *root = create_interpolated_tree2(wb, t);

	if (!E_flag) {

	    // Flag non-members (not switchable at runtime...).

	    int n_mem = 0;
	    for_all_nodes(DYN, root, bb)
		if (bb != root) {
		    if (is_member(wb, bb))
			n_mem++;
		    else
			bb->set_index(-42);		// or whatever
		}
	    // PRC(t); PRL(n_mem);
	}

	if (!max_cm_set) {
	    max_cm /= root->n_daughters();
	    max_cm_set = true;
	}

	if (root) {
	    convert_relative_to_absolute(root);

	    real dts = dt;
	    xstarplot22(root, scale, k, d, lmax,
		       point_mode, rel_point_size,
		       t, dt, D, cenergy,
		       E_flag, r_flag, f_flag, t_flag, init++,
		       step_mode, eod);
	    if (dt != dts) PRL(dt);
	    //rmtree(root);
	}

	t += dt;

	// Move to the next worldbundle, or loop, as necessary.

#define EPS 1.e-12

	if (dt > 0 && t > wb->get_t_max() + EPS) {
	    if (++ih >= nh) {
		ih = 0;
		wb = wh[ih];
		t = wb->get_t_min();
	    } else
		wb = wh[ih];
	}

	// May be moving backwards, so check that too.

	if (dt < 0 && t < wb->get_t_min() - EPS) {
	    if (--ih < 0) {
		ih = nh-1;
		wb = wh[ih];
		t = wb->get_t_max();
	    } else
		wb = wh[ih];
	}

    }
}
