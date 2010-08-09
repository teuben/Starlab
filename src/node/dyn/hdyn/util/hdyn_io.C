
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// hdyn_io:  Starlab hdyn-specific I/O functions.

#include "hdyn.h"
#include "util_io.h"
#include "../evolve/kira_defaults.h"

#ifndef TOOLBOX

// Initialize all static hdyn data here.

kira_counters *hdyn::kc				= NULL;
kira_options *hdyn::options			= NULL;
kira_diag *hdyn::diag				= NULL;

bool hdyn::use_dstar				= false;

real hdyn::stellar_encounter_criterion_sq	= 0;
real hdyn::stellar_merger_criterion_sq		= 0;
real hdyn::stellar_capture_criterion_sq		= 0;

hdyn** hdyn::perturbed_list			= NULL;
int  hdyn::n_perturbed				= 0;

real hdyn::eta					= DEFAULT_ETA;
real hdyn::eps					= DEFAULT_EPS;
real hdyn::eps2					= (DEFAULT_EPS*DEFAULT_EPS);

real hdyn::d_min_fac				= DEFAULT_D_MIN_FAC;
real hdyn::d_min_sq				= 0;	// compute at runtime
real hdyn::lag_factor				= DEFAULT_LAG_FACTOR;
real hdyn::mbar					= 0;	// compute at runtime

real hdyn::gamma2				= DEFAULT_GAMMA;
real hdyn::gamma23				= VERY_LARGE_NUMBER;

real hdyn::initial_step_limit			= 0;
real hdyn::step_limit				= 0;
real hdyn::unpert_step_limit			= 1;

real hdyn::scaled_stripping_radius		= -1;	// no value set

int hdyn::max_slow_factor			= DEFAULT_MAX_SLOW_FACTOR;
real hdyn::max_slow_perturbation		= DEFAULT_MAX_SLOW_PERTURBATION;
real hdyn::max_slow_perturbation_sq	     = DEFAULT_MAX_SLOW_PERTURBATION_SQ;

unsigned int hdyn::config			= 0;
bool hdyn::restart_grape_flag			= false;

//acc_function_ptr hdyn::kira_calculate_top_level_acc_and_jerk = NULL;

int hdyn::n_threads				= 0;
real hdyn::thread_cpu				= 0;


void check_sanity_of_timestep(hdyn *b, xreal & time, real & timestep)
{
    if (timestep > 0) {

	// Sanity check for timestep...

	if (fmod(time,  timestep) != 0.0) {

	    cerr << endl << "unsynchronized timestep error for "
		 << b->format_label() << endl;

	    cerr.precision(HIGH_PRECISION);
	    PRL(time);
	    PRL(timestep);
	    PRL(1.0/timestep);
	    PRL(fmod(time,  timestep));
	    cerr << " correcting..." << endl;

	    real corrected_timestep = 1.0;
	    while (corrected_timestep > timestep*1.1)
		corrected_timestep *= 0.5;

	    timestep = corrected_timestep;
	    real  err;
	    if ((err = fmod(time,  timestep)) != 0.0) {
		time = time - err;
		if (err > 0.5*timestep)
		    time += timestep;
	    }
	    PRL(time);
	    PRL(timestep);
	}
    }
}

static bool read_xreal = false;

void hdyn::print_static(ostream& s)		// default = cerr
{
    _dyn_::print_static(s);

    s << "kc = " << kc << endl;
    s << "options = " << options << endl;
    s << "diag = " << diag << endl;

    s << "use_dstar = " << use_dstar << endl;

    s << "stellar_encounter_criterion_sq = "
      << stellar_encounter_criterion_sq << endl;
    s << "stellar_merger_criterion_sq = "
      << stellar_merger_criterion_sq << endl;
    s << "stellar_capture_criterion_sq = "
      << stellar_capture_criterion_sq << endl;

    s << "perturbed_list = " << perturbed_list << endl;
    s << "n_perturbed = " << n_perturbed << endl;

    s << "eta = " << eta << endl;
    s << "eps = " << eps << endl;
    s << "eps2 = " << eps2 << endl;

    s << "d_min_fac = " << d_min_fac << endl;
    s << "d_min_sq = " << d_min_sq << endl;
    s << "lag_factor = " << lag_factor << endl;
    s << "mbar = " << mbar << endl;

    s << "gamma2 = " << gamma2 << endl;
    s << "gamma23 = " << gamma23 << endl;

    s << "initial_step_limit = " << initial_step_limit << endl;
    s << "step_limit = " << step_limit << endl;
    s << "unpert_step_limit = " << unpert_step_limit << endl;

    s << "scaled_stripping_radius = " << scaled_stripping_radius << endl;
  
    s << "max_slow_factor = " << max_slow_factor << endl;
    s << "max_slow_perturbation = " << max_slow_perturbation << endl;
    s << "max_slow_perturbation_sq = " << max_slow_perturbation_sq << endl;
}

// *** NOTE: check_and_correct_node() is no longer inherited from the
// *** dyn class.  For restart purposes, we don't want to change the
// *** input data in any way, other than restoring the center of mass
// *** offset (discussed below).  Thus, we can check for consistency
// *** in the snapshot file, but we do not alter it.

bool hdyn::check_and_correct_node(bool verbose)	// default = false
{
    // cerr << "hdyn::check_and_correct_node: "; PRL(this);

    bool ok = true;

    if (oldest_daughter) {		// for hdyn, check that masses,
					// pos and vel are all consistent
	real m = 0;
	vec p = 0, v = 0;
	bool low = false;

	for_all_daughters(hdyn, this, bb) {
	    if (bb->oldest_daughter)
		ok &= bb->check_and_correct_node(verbose);
	    real mm = bb->get_mass();
	    m += mm;
	    p += mm*bb->get_pos();
	    v += mm*bb->get_vel();
	}
	if (!ok) low = true;

	if (!twiddles(m, mass)) ok = false;
	if (!twiddles(abs(p), 0)) ok = false;
	if (!twiddles(abs(v), 0)) ok = false;

	mass = m;
	if (m > 0) {
	    p /= m;
	    v /= m;
	}

	if (!ok && verbose && parent == NULL) {
	    cerr << "check_and_correct_node: found ";
	    if (low)
		cerr << "low-level";
	    else
		cerr << "top-level";
	    cerr << " mass/pos/vel inconsistency" << endl;
	    PRC(mass); PRL(m);
	    PRL(p);
	    PRL(v);
	}
    }

    // Checking over.  As of 7/04, we always place the root node at rest
    // at the origin for hdyn internal data (SMcM).

    if (is_root()) offset_com();
    return ok;
}

// Allow user to turn off timestep checking...

static bool check_timestep = true;

void set_hdyn_check_timestep(bool check)	// default = true
{
    check_timestep = check;
}

istream & hdyn::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s, input_line),
//	   strcmp(END_DYNAMICS, input_line)) {
	   !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (val) {

	    // See xreal notes in dyn_io.C...

	    if (!strcmp("real_system_time", keyword)) {

		read_xreal = true;
		last_real = true;

	    } else if (!strcmp("system_time", keyword)) {

		// Check input format before reading.

		if (!last_real) read_xreal = false;

		if (read_xreal) {

		    //cerr << "hdyn::scan_dyn_story: input "
		    //     << "time data type is xreal"
		    //     << endl;

		    set_system_time(get_xreal_from_input_line(input_line));

		    // Note dangerous to separate root time and system time
		    //   -- should probably never be different...

		    // PRL(input_line);
		    // PRL(system_time);
		    // xprint(system_time);

		} else {

		    // if (sizeof(xreal) != sizeof(real))	// crude test...
		    //     cerr << "hdyn::scan_dyn_story: input "
		    //		<< "time data type is real"
		    //		<< endl;

		    real_system_time = system_time = strtod(val, NULL);

		}
		// PRC(system_time); xprint(system_time);

	    } else {

		last_real = false;

		if (!strcmp("t", keyword)) {

		    if (read_xreal)
			time = get_xreal_from_input_line(input_line);
		    else
			time = strtod(val, NULL);

		} else if (!strcmp("dt", keyword))
		    timestep = strtod(val, NULL);
		else if (!strcmp("m", keyword))
		    mass = strtod(val, NULL);
		else if (!strcmp("r", keyword))
		    set_vector_from_input_line(pos, input_line);
		else if (!strcmp("v", keyword)) {
		    set_vector_from_input_line(vel, input_line);
		    posvel = pos*vel;
		} else if (!strcmp("a", keyword))
		    set_vector_from_input_line(acc, input_line);
		else if (!strcmp("pot", keyword))
		    pot = strtod(val, NULL);
		else if (!strcmp("R_eff", keyword)) {
		    real r = strtod(val, NULL);
		    set_radius(r);
		} else if (!strcmp("steps", keyword))
		    steps = strtod(val, NULL);
		else if (!strcmp("dir_f", keyword))
		    direct_force = strtod(val, NULL);
		else if (!strcmp("indir_f", keyword))
		    indirect_force = strtod(val, NULL);

		// NOTE:  Complete initialization of unperturbed and slow
		// binary structures requires knowledge of the tree structure
		// that is available only after the entire tree has been
		// read in.  The get_hdyn() macro completes the setup of
		// these parameters after the tree is known.

		// Unperturbed motion:

		else if (!strcmp("dt_u", keyword))
		    unperturbed_timestep = strtod(val, NULL);
		else if (!strcmp("full_u", keyword))
		    fully_unperturbed = strtol(val, NULL, 10);

		// Slow binary motion:

		else if (!strcmp("slow_kappa", keyword)) {

		    // Component of a slow binary.  The information needed to
		    // recognize and reconstruct the slow structure is attached
		    // to the ELDER sister only.

		    int k = strtol(val, NULL, 10);

		    if (k > 1) {

			// Create the slow structure.  Note that we won't need
			// to apply modifications to acc or the sister data.

			slow = new slow_binary(k);
			slow->set_dtau(timestep/k);

			// t_init, t_apo, and tau will be set in due course...

		    }

		} else if (!strcmp("slow_t_init", keyword)) {

		    if (slow) {
			slow->set_t_init( strtod(val,NULL) );
		    }

		} else if (!strcmp("slow_t_apo", keyword)) {

		    if (slow) {
			slow->set_t_apo( strtod(val,NULL) );
		    }

		} else if (!strcmp("slow_tau", keyword)) {

		    if (slow) {
			slow->set_tau( strtod(val,NULL) );
			slow->init_tau_pred();
		    }

		} else

		    add_story_line(dyn_story, input_line);
	    }
	}
    }

    if (check_timestep && !fully_unperturbed)
        check_sanity_of_timestep(this, time, timestep);

    return s;
}

// Temporary flag to force unformatted output:

static bool write_unformatted = false;

// Accessors:

void set_write_unformatted(bool u)	// default = true;
{
    write_unformatted = u;
}

bool get_write_unformatted()
{
    return write_unformatted;
}

// Another kludge -- we need to know if this print_dyn_story() is
// part of a complete snapshot, or just a single node or subtree.
// Set complete_system_dump before put_node, unset it afterward.

static bool complete_system_dump = false;

void set_complete_system_dump(bool d)			// default = true
{
    complete_system_dump = d;
}

// Another (sort of) kludge...  We want the output to have the root node be
// the center of mass of the top-level nodes.  This is standard elsewhere in
// Starlab, but is not the case in hdyn/kira, for operational reasons.  Thus
// we do the offsetting at output time.  An extra complication is that tdyn
// output relies on the full dump performed earlier.  We save the system
// center of mass each time we do a complete system dump from the root node.
// The saved data are used for every pos and vel output, defining the root
// data and offsetting the daughters.

static vec fulldump_com_pos = 0, fulldump_com_vel = 0;

ostream & hdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    int precision = 0;
    int use_floats = 0;

    // Two basic branches here:
    //
    //		1. short_output (full_dump mode)
    //		2. write_unformatted
    //
    // Note that write_unformatted is only relevant when short_output > 0.
    //
    // Short_output options (see also write_unformatted):
    //
    //		0:	normal (long) output
    // 		1:	short output using current time, pos, ...
    //		2:	short output using current time, predicted pos, ...
    //		3:	as for 2, but mark as defunct
    //		4:	new (Steve, 7/01) to allow restart; combines hdyn
    //			and tdyn...  Don't allow unformatted output.
    //
    // Short output is:  time, mass, pos, vel (options 1-3).
    //
    // Particle name is printed in node/util/tree_io.C

    bool short_short = (short_output && short_output != 4);

    // See if we need to recompute the com pos and vel.

    if (complete_system_dump && is_root()) {
	compute_com(this, fulldump_com_pos, fulldump_com_vel);
	fulldump_com_pos -= pos;				// relative
	fulldump_com_vel -= vel;				// to root
    }

    if (write_unformatted && short_short) {

	// In this case, we have to do everything ourselves, and we
	// can't rely on the base class output functions to help us.

	// Note that short_output is implicit.

	if (is_root()) put_real_number(s, "  system_time  =  ",
				       real_system_time);

	// Write time, mass, pos, and vel as unformatted data.
	// Allows significantly faster I/O, and converting doubles
	// to floats saves additional space.

	// *** Must coordinate with tdyn_io.C. ***

	vec putpos = pos;
	vec putvel = vel;

	if (short_output > 1) {
	    putpos = pred_pos;
	    putvel = pred_vel;
	}

	if (is_root()) {
	    putpos += fulldump_com_pos;
	    putvel += fulldump_com_vel;
	} else if (is_top_level_node()) {
	    putpos -= fulldump_com_pos;
	    putvel -= fulldump_com_vel;
	}

	if (use_floats) {

	    // Write floats.

	    s << "t64mpv32 =" << endl;
	    write_unformatted_real(s, system_time);
	    write_unformatted32_real(s, mass);
	    write_unformatted32_vector(s, putpos);
	    write_unformatted32_vector(s, putvel);

	} else {

	    // Write doubles.

	    s << "tmpv =" << endl;
	    write_unformatted_real(s, system_time);
	    write_unformatted_real(s, mass);
	    write_unformatted_vector(s, putpos);
	    write_unformatted_vector(s, putvel);
	}

    } else {

	// For formatted output, we can rely on the _dyn_ output
	// functions to start us off as usual...

	vec save_pos = pos, save_vel = vel;
	if (is_root()) {
	    pos += fulldump_com_pos;
	    vel += fulldump_com_vel;
	} else if (is_top_level_node()) {
	    pos -= fulldump_com_pos;
	    vel -= fulldump_com_vel;
	}

	_dyn_::print_dyn_story(s, print_xreal, short_output);

	pos = save_pos;
	vel = save_vel;

	if (!short_short) {

	    put_real_number(s, "  steps  =  ", steps);
	    put_real_number(s, "  dir_f  =  ", direct_force);
	    put_real_number(s, "  indir_f  =  ", indirect_force);

	    // Special cases for which only partial information is saved:

	    if (kep) {

		// Component of an unperturbed binary -- may be fully or
		// partially unperturbed.  Store sufficient information for
		// restart.  Other kepler data can be recomputed from hdyn.
		//
		// BOTH components carry full_u and dt_u information (no
		// particular reason for this, but retain for consistency).

		put_integer(s, "  full_u  =  ", fully_unperturbed);
		put_real_number(s, "  dt_u  =  ", unperturbed_timestep);
	    }

	    if (slow) {

		// Component of a slow binary.  The value of kappa serves
		// as a flag for slow motion.  The information needed for
		// restart is attached to the ELDER component only.

		if (!elder_sister) {
		    put_integer(s, "  slow_kappa  =  ", get_kappa());
		    put_real_number(s, "  slow_t_init  =  ",
				    slow->get_t_init());
		    put_real_number(s, "  slow_t_apo  =  ",
				    slow->get_t_apo());
		    put_real_number(s, "  slow_tau  =  ", slow->get_tau());
		}
	    }
	}
    }

    if (short_output) {

	// Convenient to output Star quantities here too.
	// (Star output is now suppressed in tree_io in the
	// short_output case, and the tdyn input is simplest
	// when all relevant data are in Dyn.)

	if (use_sstar && sbase) {

	    // Print out basic stellar information.
	    // See inc/star/star_support.h for element_type.

	    // put_string(s,      "  S  =  ",
	    //		  type_short_string(sbase->get_element_type()));
	    put_integer(s,     "  S  =  ", (int)sbase->get_element_type());

	    if (write_unformatted && short_short) {

		// Always write floats for T and L.

		s << "TL =" << endl;
		write_unformatted32_real(s, sbase->temperature());
		write_unformatted32_real(s, sbase->get_luminosity());

	    } else {
		put_real_number(s, "  T  =  ", sbase->temperature());
		put_real_number(s, "  L  =  ", sbase->get_luminosity());
	    }
	}

	// Use kep as a general-purpose flag here (careful!).

	if (kep)
	    put_integer(s, "  kep  =  ", 1);
#if 1
	else {

	    // *** Hook for future treatment of lightly perturbed binaries.
	    // *** Indicator is kep = 2; currently not used by tdyn.

	    if (this != root			// must check this...
		&& is_low_level_node()
		&& get_perturbation_squared() < SMALL_TDYN_PERT_SQ)
		put_integer(s, "  kep  =  ", 2);
	}
#endif

	// Defunct node:
	
	if (short_output == 3)
	    put_integer(s, "  defunct  =  ", 1);

	// Add information on the system center to the root node.
	// Must do this explicitly here because story output is
	// suppressed in short_output mode.

	if (is_root()) {

	    vec center_pos, center_vel;
	    int which = get_std_center(this, center_pos, center_vel);
	    put_real_vector(s, "  center_pos  =  ", center_pos);
	    put_real_vector(s, "  center_vel  =  ", center_vel);
	    put_integer(s, "  center_type  =  ", which);

	    // Add extra information on other centers:

	    if (get_external_field() > 0) {
		refine_cluster_mass2(this);
		if (find_qmatch(get_dyn_story(), "bound_center_pos")) {
		    vec bound_pos = getvq(get_dyn_story(),
					     "bound_center_pos");
		    vec bound_vel = getvq(get_dyn_story(),
					     "bound_center_vel");
		    put_real_vector(s, "  bound_center_pos  =  ", bound_pos);
		    put_real_vector(s, "  bound_center_vel  =  ", bound_vel);
		}
	    }
	}

	// Extra output for complete system dumps only.

	if (complete_system_dump) {

	    if (external_field) {

		// Not needed if short_output = 4, as the data are already
		// contained in the dyn story...

		if (short_output != 4) {

		    if (!is_root()) {
			hdyn *top = get_top_level_node();

			// The esc flag is probably redundant...

			bool esc = false;
			if (find_qmatch(top->get_dyn_story(), "esc"))
			    esc = getiq(top->get_dyn_story(), "esc");
			put_integer(s, "  esc  =  ", esc);

			if (find_qmatch(top->get_dyn_story(), "t_esc"))
			    put_real_number(s, "t_esc  =  ",
					    getrq(top->get_dyn_story(),
						  "t_esc"));
		    }
		}
	    }
	}
    }

    return s;
}

typedef struct {
    hdyn *b;
    real  d;
} bd_pair, *bd_pair_ptr;

local int compare(const void *pi, const void *pj)
{
    if (((bd_pair_ptr) pi)->d > ((bd_pair_ptr) pj)->d)
        return 1;
    else if (((bd_pair_ptr) pi)->d < ((bd_pair_ptr) pj)->d)
        return -1;
    else
        return 0;
}

//#define SORT_PERTURBERS

void hdyn::print_perturber_list(ostream & s, const char *pre)
{
    // NOTE: 'this' is a binary center of mass.

    s << pre << "perturber list of " << format_label()
      << " [" << perturber_list << "]" << flush;

    if (!valid_perturbers || perturber_list == NULL) {
	s << " is invalid" << endl;
	return;
    }

    int np = n_perturbers;

    if (np <= 0) {
	s << " is empty (no perturbers)" << endl;
	return;
    }

    s << "  (" << np << " perturbers):";
    if (np > MAX_PERTURBERS) np = MAX_PERTURBERS;	// shouldn't happen

    s << endl << pre << "    ";

#ifndef SORT_PERTURBERS

    // Just write out the list without sorting it.

    for (int i = 0; i < np; i++) {
        hdyn *p = perturber_list[i];
	if (p && p->is_valid())
	    p->print_label(s);
	else
	    s << "(null)";
	if ((i+1)%10 == 0) {
	    if (i < np-1)
		s << endl << pre << "    ";
	} else
	    s << " ";
    }

#else

    // Sort the perturbers by distance (between top-level nodes) before
    // writing out the list.

    // Note from Steve (8/03): for unknown reasons, the following code
    // seems to be unstable under g++ 3.2.2 with optimization -O2 under
    // RedHat 9.  The symptom is that subtle changes are apparently made
    // in some members of the list...  The result is that invocation of
    // this function affects the outcome of the simulation.  The dynamical
    // declaration (with dimension MAX_PERTURBERS or np) causes the
    // problem by itself.  The static declaration seems to work (note the
    // commented-out delete statement at the end), but the code still
    // seems to break when the qsort below is restored.  Without qsort,
    // we may as well use the simpler code above (which is OK...).

    // bd_pair_ptr bd_list = new bd_pair[MAX_PERTURBERS];
    bd_pair bd_list[MAX_PERTURBERS];

    hdyn *top = get_top_level_node();

    for (int i = 0; i < np; i++) {
	bd_list[i].b = perturber_list[i];
	if (perturber_list[i] && perturber_list[i]->is_valid()) {
	    bd_list[i].d = square(perturber_list[i]
					  ->get_top_level_node()->pos
				  - top->pos);
	} else
	    bd_list[i].d = VERY_LARGE_NUMBER;
    }

    // qsort((void *)bd_list, (size_t)np, sizeof(bd_pair), compare);

    for (int i = 0; i < np; i++) {
        hdyn *p = bd_list[i].b;
	if (p && p->is_valid())
	    p->print_label(s);
	else
	    s << "(null)";
	if ((i+1)%10 == 0) {
	    if (i < np-1)
		s << endl << pre << "    ";
	} else
	    s << " ";
    }

    // delete [] bd_list;

#endif

    s << endl;

    // s << "leaving print_perturber_list" << endl << flush;
}

void hdyn::find_print_perturber_list(ostream & s, const char* pre)
{
    s << pre << "perturber list node for " << format_label() << " is ";

    hdyn* pnode = find_perturber_node();

    if (!pnode) {
	s << "NULL" << endl;
	return;
    } else
	s << pnode->format_label() << endl;

    pnode->print_perturber_list(s, pre);
}

#else

main(int argc, char** argv)
{
    hdyn  * b;
    check_help();

    while (b = get_hdyn()) {
	cout << "TESTING put_hdyn:" << endl;
        put_node(b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	delete b;
    }
    cerr << "Normal exit" << endl;
}

#endif
