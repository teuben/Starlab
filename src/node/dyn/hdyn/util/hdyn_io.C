
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// hdyn_io:  Starlab hdyn-specific I/O functions.

#include "hdyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize all static hdyn data here (?)

kira_counters *hdyn::kc		= NULL;    
kira_options *hdyn::options	= NULL;    
kira_diag *hdyn::diag		= NULL;    

int  hdyn::tidal_type	= 0;
real hdyn::omega	= 0;
real hdyn::omega_sq	= 0;
real hdyn::alpha1	= 0;
real hdyn::alpha3	= 0;

bool hdyn::use_dstar	= false;

real hdyn::stellar_encounter_criterion_sq = 1;
real hdyn::stellar_merger_criterion_sq = 1;
real hdyn::stellar_capture_criterion_sq = 1;

hdyn** hdyn::perturbed_list = NULL;
int  hdyn::n_perturbed	= 0;

real hdyn::eta		= 0;
real hdyn::eps		= 0;
real hdyn::eps2		= 0;

real hdyn::d_min_sq	= 0;
real hdyn::lag_factor	= 0;
real hdyn::mbar		= 0;

real hdyn::gamma2	= 0;    
real hdyn::gamma23	= VERY_LARGE_NUMBER;   

real hdyn::initial_step_limit	= 0;
real hdyn::step_limit		= 0;
real hdyn::unpert_step_limit	= 1;

real hdyn::scaled_stripping_radius = -1;

int hdyn::max_slow_factor = 1;			// default is to suppress slow
						// motion, for now
real hdyn::max_slow_perturbation = 1.e-4;	// should be tied to parameters
						// for unperturbed motion
real hdyn::max_slow_perturbation_sq = 1.e-8;


void check_sanity_of_timestep(xreal & time, real & timestep)
{
    if (timestep > 0) {

	// Sanity check for timestep...

	if (fmod(time,  timestep) != 0.0) {

	    cerr << " impossible timestep error \n";

	    cerr.precision(HIGH_PRECISION);
	    PRL(time);
	    PRL(timestep);
	    PRL(1.0/timestep);
	    PRL(fmod(time,  timestep));
	    cerr << " Perform correction... \n";
      
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

istream & hdyn::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s, input_line), strcmp(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

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
	    } else {
		cerr << "hdyn::scan_dyn_story: input "
		     << "time data type is real"
		     << endl; 

		real_system_time = system_time = strtod(val, NULL);
	    }
	    // PRC(system_time); xprint(system_time);

	} else {

	    last_real = false;

	    if (!strcmp("t", keyword)) {

		if (read_xreal)
		    time = get_xreal_from_input_line(input_line);
		else {
		    time = strtod(val, NULL);
		}

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
	    else if (!strcmp("R_eff", keyword))
		radius = strtod(val, NULL);
	    else if (!strcmp("steps", keyword))
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

    check_sanity_of_timestep(time, timestep);

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

ostream & hdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0
{
    put_story_header(s, DYNAMICS_ID);

    int precision = 0;
    int use_floats = 0;

    // Short_output options (see also write_unformatted):
    //
    //		0:	normal (long) output
    // 		1:	short output using current time, pos, ...
    //		2:	short output using current time, predicted pos, ...
    //		3:	as for 2, but mark as defunct
    //
    // Note that write_unformatted is only relevant when short_output > 0.
    // Particle name is printed in node/util/tree_io.C

    // **** Output formats are getting out of hand!  (Steve, 10/00) ****

    if (!parent) {				// root node

	// See xreal notes in dyn_io.C...

#ifdef USE_XREAL

	if (print_xreal && !short_output) {

	    put_real_number(s, "  real_system_time  =  ", (real)system_time);
	    put_real_number(s, "  system_time  =  ", system_time);

	} else
		
	    put_real_number(s, "  system_time  =  ", (real)system_time);

#else

	put_real_number(s, "  system_time  =  ", system_time);

#endif

    }

    if (short_output) {

	if (write_unformatted) {

	    // Write time, mass, pos, and vel as unformatted data.
	    // Should allow significantly faster I/O, and converting
	    // doubles to floats will save additional space.

	    // *** Must coordinate with tdyn_io.C. ***

//  	    vector putpos = something_relative_to_root(this, &hdyn::get_pos);
//	    vector putvel = something_relative_to_root(this, &hdyn::get_vel);
	    vector putpos = pos;
	    vector putvel = vel;

	    if(short_output == 2) {
		putpos = pred_pos;
		putvel = pred_vel;
	    }

	    if(use_floats) {
		s << "t64mpv32 =" << endl;
		write_unformatted_real( s, system_time );
		write_unformatted32_real( s, mass );
		write_unformatted32_vector( s, putpos );
		write_unformatted32_vector( s, putvel );
	    } else {
		// Write doubles.
		s << "tmpv =" << endl;
		write_unformatted_real( s, system_time );
		write_unformatted_real( s, mass );
		write_unformatted_vector( s, putpos );
		write_unformatted_vector( s, putvel );
	    }

	    // Use kep as a flag here (careful!).

	    if (kep)
		put_integer(s, "  kep  =  ", 1);

	} else {

	    // Always print high-precision real system time for t.

	    adjust_starlab_precision(HIGH_PRECISION);
	    put_real_number(s, "  t  =  ", (real)system_time);
	    adjust_starlab_precision(-1);
	}

    } else {

	if (print_xreal)
	    put_real_number(s, "  t  =  ", time);	// OK for real or xreal
	else
	    put_real_number(s, "  t  =  ", (real)time);

	put_real_number(s, "  dt =  ", timestep);
    }

    if (!write_unformatted) {

	// (Control output precision with the Starlab "precision" variable.

	if (short_output > 1) {

	    // Use predicted quantities.


	    put_real_number(s, "  m  =  ", mass);
	    put_real_vector(s, "  r  =  ", pred_pos);
	    put_real_vector(s, "  v  =  ", pred_vel);

	    if (kep) put_integer(s, "  kep  =  ", 1);

	    // put_real_vector(s, "  a  =  ",			// not used in
	    //		   acc + (system_time-time)*jerk);	// short output

	} else {

//  	    vector putpos = something_relative_to_root(this, &hdyn::get_pos);
//	    vector putvel = something_relative_to_root(this, &hdyn::get_vel);
  	    vector putpos = get_pos();
  	    vector putvel = get_vel();

	    put_real_number(s, "  m  =  ", mass);
	    put_real_vector(s, "  r  =  ", putpos);
	    put_real_vector(s, "  v  =  ", putvel);

	    if (!short_output) {

		put_real_vector(s, "  a  =  ", acc);		// not used in
								// short output
		put_real_number(s, "  pot  =  ", pot);
		put_real_number(s, "  R_eff  =  ", radius);
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

		if (dyn_story)
		    put_story_contents(s, *dyn_story);
	    }
	}
    }

    // For use with full_dump mode...

    if (short_output == 3)
	put_integer(s, "  defunct  =  ", 1);

    put_story_footer(s, DYNAMICS_ID);

    return s;
}

typedef struct
{
    hdyn* b;
    real  d;
} bd_pair, *bd_pair_ptr;

local int compare(const void * pi, const void * pj)
{
    if (((bd_pair_ptr) pi)->d > ((bd_pair_ptr) pj)->d)
        return 1;
    else if (((bd_pair_ptr)pi)->d < ((bd_pair_ptr)pj)->d)
        return -1;
    else
        return 0;
}

void hdyn::print_perturber_list(ostream & s, char* pre)
{
    // NOTE: 'this' is a binary center of mass.

    s << pre << "perturber list of " << format_label() << flush;

    if (!valid_perturbers || perturber_list == NULL) {
	s << " is invalid" << endl;
	return;
    }

    int np = n_perturbers;

    if (np <= 0) {
	s << " is empty (no perturbers)" << endl;
	return;
    }

    // Sort the perturbers by distance (between top-level nodes).

    hdyn* top = get_top_level_node();

    // s << "creating list" << endl << flush;

    bd_pair_ptr bd_list = new bd_pair[np];
    for (int i = 0; i < np; i++) {
	bd_list[i].b = perturber_list[i];
	if (perturber_list[i] && perturber_list[i]->is_valid()) {
	    bd_list[i].d = square(perturber_list[i]
					  ->get_top_level_node()->pos
				  - top->pos);
	} else
	    bd_list[i].d = VERY_LARGE_NUMBER;
    }

    // s << "sorting list" << endl << flush;

    qsort((void *)bd_list, (size_t)np, sizeof(bd_pair), compare);

    s << "  (" << np << " perturbers):";
    s << endl << pre << "    ";
    for (int i = 0; i < np; i++) {
	if (bd_list[i].b)
	    bd_list[i].b->print_label(s);
	else
	    s << "(null)";
	if ((i+1)%10 == 0) {
	    if (i < np-1)
		s << endl << pre << "    ";
	} else
	    s << " ";
    }
    s << endl;

    delete [] bd_list;

    // s << "leaving print_perturber_list" << endl << flush;
}

void hdyn::find_print_perturber_list(ostream & s, char* pre)
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

#endif
