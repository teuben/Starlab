
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// kira_init.C:  Read in command-line parameters and initialize
//		 the system.
//
// Externally visible functions:
//
//	bool kira_initialize

#include "star/dstar_to_kira.h"
#include "kira_defaults.h"

#define		LONG_BINARY_OUT		10
#define		SEED_STRING_LENGTH	60

#define Rsun_pc 2.255e-8	// R_sun/1 parsec = 6.960e+10/3.086e+18;

local void initialize_index(node * b, bool verbose)
{
    // Give the root node a name if one hasn't already been assigned.

    if (b->is_root() && b->get_index() < 0 && !b->get_name())
	b->set_name("root");

    // Add indices to otherwise unidentified nodes.

    int index = 0;

    for_all_nodes(node, b, bb)
	if (bb != b)
	    if (bb->get_index() > 0)
		index = Starlab::max(index, bb->get_index());

    // Start numbering at some well-separated value.

    index = 1000 * (index / 1000 + 1);

    for_all_nodes(node, b, bb)
	if (bb != b && bb->get_index() < 0 && bb->get_name() == NULL) {
	    if (verbose) {
		cerr << "attached index " << index << " to unnamed ";
		if (bb->is_parent())
		    cerr << "node\n";
		else
		    cerr << "leaf\n";
	    }
	    bb->set_index(index++);
	}
}

local void correct_multiples(hdyn* b,
			     bool verbose)
{
    // Check and fix input problems associated with hierarchical systems.
    // Probably not the best place to do this, but keep it here for now...

    for_all_nodes(hdyn, b, bb) {
	if (bb->get_elder_sister() == NULL
	    && bb->get_kepler() != NULL) {
	    if (bb->is_parent()) {
		if (verbose) {
		    cerr << "unperturbed " << bb->format_label() << endl;
		    cerr << "    multiple system, ";
		}
		if (bb->get_oldest_daughter()->get_kepler() == NULL) {
		    if (verbose)
			cerr << "daughter node "
			     << bb->get_oldest_daughter()->format_label()
			     << " not unperturbed.  Correcting...\n";
		    hdyn *od = bb->get_oldest_daughter();
		    od->reinitialize_kepler_from_hdyn();
		    od->recompute_unperturbed_step();
		} else
		    if (verbose)
			cerr << "daughter node uperturbed.\n";
	    }
	}
    }
}

local void choose_param(hdyn* b, bool verbose,
			real& x, bool x_flag, char* x_id,
			bool zero_OK = false)
{
    // Set kira parameter x, subject to the twin constraints of
    // command-line input and input via the incoming snapshot.

    static bool skip_line = true;	// formatting!

    char kira_x[128];
    strcpy(kira_x, "kira_");
    strcat(kira_x, x_id);

    bool x_in_snap = find_qmatch(b->get_log_story(), kira_x);

    if (x_flag) {

	// Value of parameter has been set on the command line.
	// Command line value takes precedence.
	// Flag the existence of a snapshot version only in verbose mode.

	if (x_in_snap && verbose) {
	    if (skip_line) cerr << endl;
	    cerr << "command-line " << x_id << " = " << x;
	    if (twiddles(x, getrq(b->get_log_story(), kira_x)))
		cerr << " is identical to value";
	    else
		cerr << " overrides value "
		     << getrq(b->get_log_story(), kira_x);
	    cerr << " in input snapshot"
		 << endl;
	    skip_line = false;
	}

    } else {

	if (x_in_snap) {

	    real x_snap = getrq(b->get_log_story(), kira_x);

	    if ((x_snap > 0 || x_snap == 0 && zero_OK)) {

		x = x_snap;
		if (verbose) {
		    if (skip_line) cerr << endl;
		    cerr << "taking " << x_id << " = " << x
			 << " from input snapshot"
			 << endl;
		    skip_line = false;
		}

	    } else {

		if (verbose) {
		    if (skip_line) cerr << endl;
		    cerr << "Error reading " << x_id
			 << " from input snapshot;"
			 << " using default = " << x
			 << endl;
		    skip_line = false;
		}
	    }

	} else {

	    if (verbose) {
		if (skip_line) cerr << endl;
		cerr << "using default " << x_id << " = " << x
		     << endl;
		skip_line = false;
	    }
	}
    }

    if (x < 0 || (x == 0 && !zero_OK)) {
	char tmp[128];
	strcpy(tmp, "Error setting ");
	strcat(tmp, x_id);
	err_exit(tmp);
    }

    // Save the setting in the log story for future use.

    putrq(b->get_log_story(), kira_x, x);
}

local void choose_param(hdyn* b, bool verbose,
			int& ix, bool x_flag, char* x_id,
			bool zero_OK = false)
{
    // Use real version to achieve desired effect.  Rely on I/O
    // formatting to make ints look like ints on output.

    real x = ix;
    choose_param(b, verbose, x, x_flag, x_id, zero_OK);
    ix = (int)(x+0.1);
}

local void set_runtime_params(hdyn *b, bool verbose,
			      real eta, bool eta_flag,
			      real eps, bool eps_flag,
			      real d_min, bool d_min_flag,
			      real lag_factor, bool lag_flag,
			      real gamma, bool gamma_flag,
			      int max_slow, bool max_slow_flag)
{
    // Procedure:  snapshot data may override defaults,
    //		   command-line input will override both.

    // Could read/write the static data as part of get/put_hdyn, but
    // keep it this way for now (Steve, 7/99).

    choose_param(b, verbose, eta, eta_flag, "eta");
    b->set_eta(eta);

    choose_param(b, verbose, eps, eps_flag, "eps", true);
    b->set_eps(eps);				// ^^^^	(0 is OK)

    choose_param(b, verbose, d_min, d_min_flag, "d_min_fac");
    b->set_d_min_fac(d_min);

    // Create a value for d_min_sq().  The value will be revised
    // in evolve_system().

    int nbody = b->n_daughters();	// note scaling with *current* N
					// and assumption that rvirial = 1
    real d_min_sq = square(d_min/nbody);

    // See if a value exists in the snapshot.

    choose_param(b, verbose, d_min_sq, false, "d_min_sq");
    b->set_d_min_sq(d_min_sq);

    // Rudimentary check: flag values of d_min_sq very different
    // from the equilibrium value.

    real ratio = d_min_sq / square(d_min/nbody);
    if (ratio < 0.0002 || ratio > 5.e3)
	cerr << endl
	     << "*** warning: d_min_sq ratio = " << ratio << " ***"
	     << endl;

    choose_param(b, verbose, lag_factor, lag_flag, "lag_factor");
    b->set_lag_factor(square(lag_factor));

    choose_param(b, verbose, gamma, gamma_flag, "gamma");
    b->set_gamma2(square(gamma));

    choose_param(b, verbose, max_slow, max_slow_flag, "log_max_slow", true);
    if (max_slow >= 0) {
	int kappa_max = (int) pow(2.0, max_slow);
	cerr << "setting maximum slowdown factor = " << kappa_max << endl;
	b->set_max_slow_factor(kappa_max);
	if (kappa_max > 1)
	    cerr << "maximum slow perturbation = "
		 << b->get_max_slow_perturbation() << endl;
	putiq(b->get_log_story(), "log_max_slow", max_slow);
    }
}

local real get_scaled_stripping_radius(hdyn* b,
				       bool verbose,
				       real input_stripping_radius,
				       real initial_r_jacobi,
				       real initial_r_virial,
				       real initial_mass)
{
    // Kira has the options of stripping stars beyond some specified
    // radius (referred to as the "stripping" radius in the code).  The
    // stripping radius r is set initially with "-G r" on the command line.
    //
    // Implementation of the stripping is INDEPENDENT of the treatment of
    // the tidal field.  However, the interpretation of r depends on what
    // tidal information is available.
    //
    // 1. If a tidal field has not been specified, r specifies the stripping
    // radius in units of the initial "tidal" radius, if known.  Otherwise,
    // it specifies the stripping radius in units of the virial radius.
    //
    // 2. If a tidal field has been specified, r specifies the stripping
    // radius in units if the jacobi radius.

    real scaled_stripping_radius = 0;	// stripping radius for unit mass
    					// (0 ==> no stripping)

    // See if a kira_scaled_stripping_radius exists in the input file.
    // If it does, it TAKES PRECEDENCE over any other setting.

     if (find_qmatch(b->get_log_story(), "kira_scaled_stripping_radius")) {

	real kira_scaled_stripping_radius
	    = getrq(b->get_log_story(), "kira_scaled_stripping_radius");

	if (kira_scaled_stripping_radius > 0) {

	    scaled_stripping_radius = kira_scaled_stripping_radius;

	    if (verbose) {
		cerr << "using scaled stripping radius ";
		PRL(kira_scaled_stripping_radius);
		cerr << "    from input snapshot" << endl;

		if (input_stripping_radius > 0)
		    cerr << "ignoring \"-G " << input_stripping_radius
			 << "\" found on command line"
			 << endl;
	    }

	} else {

	    if (verbose)
		cerr << "warning: error reading "
		     << "kira_initial_jacobi_radius"
		     << " from input snapshot"
		     << endl;
	    else
		err_exit(
	     "Error reading kira_initial_jacobi_radius from input snapshot");

	}
    }

    if (scaled_stripping_radius <= 0) {

	// See if we can determine the scaled stripping radius from
	// the input snapshot and command-line data.

	if (input_stripping_radius > 0) {

	    if (initial_r_jacobi > 0) {

		if (verbose)
		    cerr << "command-line -G " << input_stripping_radius
			 << " used as scaling factor for initial Jacobi radius"
			 << endl;

		input_stripping_radius *= initial_r_jacobi;

	    } else {

		// See if we have any tidal information in the input data.

		real r_jacobi_over_r_virial
		    = getrq(b->get_log_story(),
			    "initial_rtidal_over_rvirial");

		if (r_jacobi_over_r_virial > 0) {

		    input_stripping_radius *= r_jacobi_over_r_virial;

		    if (verbose)
			cerr << "command-line -G " << input_stripping_radius
			     << " used as scaling factor for tidal radius"
			     << endl;

		} else

		    if (verbose)
			cerr << "command-line -G " << input_stripping_radius
			     << " used as scaling factor for virial radius"
			     << endl;

		input_stripping_radius *= initial_r_virial;

	    }

	    cerr << "rescaled "; PRL(input_stripping_radius);

	    scaled_stripping_radius = input_stripping_radius
					/ pow(initial_mass, 1.0/3.0);

	    // Save scaled_stripping_radius for restart.

	    if (scaled_stripping_radius > 0)
		putrq(b->get_log_story(), "kira_scaled_stripping_radius",
		      scaled_stripping_radius);
	}
    }

    return scaled_stripping_radius;
}

local void get_physical_scales(hdyn* b,
			       bool verbose,
			       real initial_mass,
			       real initial_r_virial,
			       bool M_flag, real m_tot,
			       bool R_flag, real r_vir,
			       bool T_flag, real t_vir,
			       real q_vir,
			       int nbody)
{
    // (Re)determine physical scales and set conversion factors.

    // ***  Redundant in new version...  ***

    // First derive the N-body time unit.

    real initial_t_virial
	= sqrt(2*q_vir*pow(initial_r_virial, 3) / initial_mass);

    // Some or all of the scaling factors were not specified in
    // the input snapshot, or will be overridden by command-line
    // input. The following will be non-negative if the information
    // is already known.

    // Total mass in solar units:

    real old_m_tot
	= b->get_starbase()->conv_m_dyn_to_star(initial_mass);

    // Virial radius in parsecs (note additional scaling, since
    // SeBa prefers to work in solar radii):

    real old_r_vir
	= b->get_starbase()->conv_r_dyn_to_star(initial_r_virial)
	    * Rsun_pc;

    // Virial time scale in Myr:

    real old_t_vir
	= b->get_starbase()->conv_t_dyn_to_star(initial_t_virial);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Note: at this point,
    //
    //	initial_mass	  is the initial system mass, in N-body units
    //	initial_r_virial  is the initial virial radius, in N-body units
    //	initial_t_virial  is the initial virial time scale, in
    //			  N-body units
    //
    // We now determine the physical time units
    //
    //	m_tot		  is the initial system mass, in solar masses
    //	r_vir		  is the initial virial radius, in pc
    //	t_vir		  is the initial virial time scale, in Myr

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    bool M_set = M_flag;

    if (M_set) {
	if (verbose)
	    if (old_m_tot > 0) {
		cerr << "command-line -M " << m_tot;
		if (twiddles(m_tot, old_m_tot))
		    cerr << " is identical to value";
		else
		    cerr << " overrides value " << old_m_tot;
		cerr << " in input snapshot" << endl;
	    } else
		cerr << "taking total mass = " << m_tot
		     << " Msun from command line" << endl;
    } else
	if (old_m_tot > 0) {
	    m_tot = old_m_tot;
	    cerr << "taking total mass = " << m_tot
		 << " Msun from input snapshot" << endl;
	    M_set = true;
	}

    // Note (Steve, 6/99): removed m_flag and m_bar_true options.
    // See 990611kira_init.C for old code.

    if (!M_set) err_exit("Physical mass not specified");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    bool R_set = R_flag;

    if (R_set) {
	if (verbose) {
	    if (old_r_vir > 0) {
		cerr << "command-line -R " << r_vir;
		if (twiddles(r_vir, old_r_vir))
		    cerr << " is identical to value";
		else
		    cerr << " overrides value " << old_r_vir;
		cerr << " in input snapshot" << endl;
	    } else
		cerr << "taking virial radius = " << r_vir
		     << " pc from command line" << endl;
	}
    } else
	if (old_r_vir > 0) {
	    r_vir = old_r_vir;
	    cerr << "taking virial radius = " << r_vir
		 << " pc from input snapshot" << endl;
	    R_set = true;
	}

    if (!R_set && !T_flag && old_t_vir <= 0)
	err_exit("Physical virial radius not specified");

    // (We will compute r_vir from the definition of the time
    // unit below if it is not otherwise specified.)

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (!T_flag) {

	// See if the time unit was specified in the input snapshot.
	// Note, however, that if M and/or R were specified on the
	// command line, then the snapshot version is probably not
	// valid, and should be recomputed from the other quantities.

	if (old_t_vir > 0 && old_m_tot > 0 && old_r_vir > 0
	    && !M_flag && !R_flag) {

	    t_vir = old_t_vir;
	    if (verbose) cerr << "taking virial time scale = " << t_vir
			      << " from input snapshot" << endl;

	} else {

	    // Derive the time unit from the other units.
	    // (Note that r_vir must be known if we reach this point.)

	    // Standard (N-body / Heggie & Mathieu) time unit, in Myr:

	    t_vir = 21.0 * sqrt(q_vir*pow(r_vir, 3)/m_tot);

	    if (verbose) {
		if (old_t_vir > 0) {
		    cerr << "recomputed virial time scale " << t_vir;
		    if (twiddles(t_vir, old_t_vir))
			cerr << " is identical to value";
		    else
			cerr << " overrides value " << old_t_vir;
		    cerr << " in input snapshot" << endl;
		} else
		    cerr << "computed virial time scale = " << t_vir
			 << endl;
	    }
	}

    } else if (verbose) {

	if (old_t_vir > 0) {
	    cerr << "command-line -T " << t_vir;
	    if (twiddles(t_vir, old_t_vir))
		cerr << " is identical to value";
	    else
		cerr << " overrides value " << old_t_vir;
	    cerr << " in input snapshot" << endl;
	} else
	    cerr << "taking virial time scale = " << t_vir
		 << " from command line" << endl;

    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (!R_set || (!R_flag && T_flag)) {

	// Either r_vir has not yet been set, or it was read from
	// the input snapshot and t_vir has since been set on the
	// command line.  In either case, recompute r_vir for
	// consistency.

	r_vir = pow(t_vir / (21.0 * sqrt(q_vir/m_tot)), 2.0/3);

	if (verbose) {
	    if (!R_set)
		cerr << "computed";
	    else
		cerr << "recomputed";
	    cerr << " virial radius = " << r_vir << endl;
	}
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Note from Steve (6/99): use_relax is not presently supported
    // and has been removed.  See 990611kira_init.C for old code.

    if (verbose) {
	cerr << endl;
	cerr << "Units (N-body / physical):" << endl;
	cerr << "    <m>: " << initial_mass/nbody
	     << " / " << m_tot/nbody << " M_sun" << endl;
	cerr << "    [M]: " << initial_mass
	     << " / " << m_tot << " M_sun" << endl;
	cerr << "    [T]: " << initial_t_virial << " / " << t_vir
	     << " Myr" << endl;
	cerr << "    [L]: " << initial_r_virial << " / " << r_vir
	     << " pc"  << endl;
    }

    // Finally, check the consistency of the mass, length, and time
    // scales.  Note that no consistency check is performed if the
    // scales are taken unaltered from the input snapshot (i.e. we
    // never enter this function).

    if (verbose) {
	real tv = 21.0 * sqrt(q_vir*pow(r_vir, 3)/m_tot);
	if (abs(1 - t_vir/tv) > 1.e-4) {
	    cerr << "warning: inconsistent units: ";
	    PRC(t_vir); PRL(tv);
	}
    }

    // Apply the new scaling.

    b->get_starbase()
     ->set_stellar_evolution_scaling(m_tot/initial_mass,
				     r_vir/initial_r_virial,
				     t_vir/initial_t_virial);
}

local char* stredit(char* s, char c1, char c2)	// duplicate in runtime_help.C
{
    if (!s) return NULL;

    char* s1 = new char[strlen(s)+1];
    if (!s1) return NULL;

    strcpy(s1, s);

    char* ss = s1;
    while (*ss != '\0') {
	if (*ss == c1) *ss = c2;
	ss++;
    }

    return s1;
}

#include <time.h>

local void kira_system_id(int argc, char** argv)
{
    // Identify the Starlab version.

    cerr << "Starlab version " << STARLAB_VERSION;

    char* s = stredit(_COMPILE_DATE_, '_', ' ');
    if (s) {
	cerr << "; kira created on " << s;
	delete s;
    }
    cerr << endl;

    // Attempt to identify the user and host.

    cerr << "kira run by user ";
    if (getenv("LOGNAME"))
	cerr << getenv("LOGNAME");
    else
	cerr << "(unknown)";
    cerr << " on host ";
    if (getenv("HOST"))
	cerr << getenv("HOST");
    else if (getenv("HOSTNAME"))
	cerr << getenv("HOSTNAME");
    else
	cerr << "(unknown)";
    if (getenv("HOSTTYPE"))
	cerr << " (host type " << getenv("HOSTTYPE") <<")";
    cerr << endl;

    // Also log the time and date.

    const time_t tt = time(NULL);
    cerr << "    on  " << ctime(&tt);

    cerr << "command line reference:  " << argv[0] << endl;
    cerr << "arguments: ";
    for (int i = 1; i < argc; i++) {
	cerr << " " << argv[i];
	if (i%20 == 0) cerr << endl << "           ";
    }
    cerr << endl;

    if (   (argv[0][0] == '.' && argv[0][1] == '/')	// easier to list
	|| (argv[0][0] == '~' && argv[0][1] == '/')	// excluded strings...
	||  argv[0][0] == '/') {
    } else {

	char tmp[256];
	sprintf(tmp, "which %s > ./KIRA_TEMP", argv[0]);
	// cerr << tmp << endl;

	system(tmp);		// Note: "system" adds newline and insists
				// on writing to stdout, so we must capture
				// the output in a file and read it back...

	ifstream file("./KIRA_TEMP");
	if (file) {

	    file >> tmp;
	    cerr << "system `which " << argv[0] << "` = " << tmp << endl;

	    file.close();
	    sprintf(tmp, "rm ./KIRA_TEMP");
	    system(tmp);	
	}

    }

    cerr << "current directory = " << getenv("PWD") << endl;
    cerr << endl;
}

#define MASS_TOL 1.e-12

local void check_total_mass(hdyn *b, bool reset = true)
{
    // Check that the mass of b is equal to the sum of its leaves,
    // and optionally correct if not the case.

    real total_mass = 0;
    for_all_leaves(hdyn, b, bi) total_mass += bi->get_mass();

    if (total_mass > 0) {
	if (abs(b->get_mass()/total_mass - 1) > MASS_TOL) {
	    cerr << endl
		 << "*** Root mass disagrees with total mass"
		 << endl;
	    real dm = b->get_mass() - total_mass;
	    PRI(4); PRC(b->get_mass()); PRC(total_mass); PRL(dm);
	    if (reset) {
		cerr << "    resetting..." << endl;
		b->set_mass(total_mass);
	    }
	}
    }
}


bool kira_initialize(int argc, char** argv,
		     hdynptr& b,	// hdyn root node
		     real& delta_t,	// time span of the integration
		     real& dt_log,	// output interval--power of 2
		     int&  long_binary_out, // frequency of full binary info
		     real& dt_snap,	// snap output interval
		     real& dt_sstar,	// stellar evolution timestep
		     real& dt_esc,	// escaper removal
		     real& dt_reinit,	// reinitialization interval
		     real& dt_fulldump,	// full dump interval
		     bool& exact,	// no perturber list if true
		     real& cpu_time_limit,
		     bool& verbose,
		     bool& save_snap_at_log,
		     char* snap_save_file,
		     int& n_stop,	// n to terminate simulation
		     bool& alt_flag)	// enable alternative output
{

    // (Probably should create a class/structure to cope with these
    // command-line arguments -- Steve, 7/02.)

    // Establish defaults for parameters (static class members or
    // carried back explicitly).  Probably not really necessary to
    // initialize here, as class variables have already been set
    // in util/hdyn_io.C.

    real eta = DEFAULT_ETA;
    real eps = DEFAULT_EPS;
    real d_min = DEFAULT_D_MIN_FAC;
    real lag_factor = DEFAULT_LAG_FACTOR;
    real gamma = DEFAULT_GAMMA;

    bool eta_flag = false, eps_flag = false, d_min_flag = false,
         lag_flag = false, gamma_flag = false;

    delta_t = DEFAULT_DT;

    dt_log = DEFAULT_DT_LOG;
    dt_snap = VERY_LARGE_NUMBER;
    dt_sstar = DEFAULT_DT_SSTAR;
    dt_esc = VERY_LARGE_NUMBER;
    dt_reinit = DEFAULT_DT_REINIT;
    dt_fulldump = DEFAULT_DT_FULLDUMP;

    bool binary_zero = false;

    // Frequency of full binary output, in units of the log output interval.

    long_binary_out = LONG_BINARY_OUT;

    bool log_flag = false, snap_flag = false, sstar_flag = false,
         esc_flag = false, reinit_flag = false, fulldump_flag = false;

    cpu_time_limit = VERY_LARGE_NUMBER;

    exact = false;
    verbose = true;

    save_snap_at_log = false;
    n_stop = DEFAULT_N_STOP;

    alt_flag = DEFAULT_ALT;

    int max_slow = 0;
    bool max_slow_flag = false;

    // Defaults for ancillary and indirect parameters:

    real input_stripping_radius = 0;
    char *comment;		// comment string

    int  input_seed, actual_seed;

    // For use by stellar evolution (establish units and conversions):

    real m_tot = 1;		// actual total mass, in solar masses
    real r_vir  = 1;            // radius scaling (virial radius in parsecs)
    real t_vir = 1;             // time scaling (system time unit, in Myr)
    real q_vir = 0.5;           // virial ratio (instead of time scaling)
    real T_start = 0;

    real friction_beta = 0;	// "BT" value = 1

    real sec = 0, smc = 0, stc = 0;
    bool sec_flag = false;
    bool smc_flag = false;
    bool stc_flag = false;

    bool B_flag = false;
    bool G_flag = false;
    bool S_flag = false;
    bool c_flag = false;
    bool s_flag = false;

    bool allow_kira_override = true;

    bool ignore_internal = false;

    bool r_virial_set = false;
    real input_r_virial;

    char seedlog[SEED_STRING_LENGTH];

    extern char *poptarg, *poparr[];	// multiple arguments are allowed
					// as of 8/99 (Steve)
    int c;
    char* param_string =
"*:a:Ab.Bc:C:d:D:e:E:f:F.g:G:h:iI:k:K:L:m:M:n:N:oO:q:Qr:R:s:St:T:uUvVW:xX:y:z:Z:";

   // ^	optional (POSITIVE!) arguments are allowed as of 8/99 (Steve)

    kira_system_id(argc, argv);

    // Make a preliminary pass over the argument list to check for
    // an input file:

    char infile[256];
    bool read_from_file = false;

    for (int i = 0; i < argc; i++)
	if (argv[i][0] == '-' && argv[i][1] == 'R') {
	    strcpy(infile, argv[++i]);
	    read_from_file = true;
	}

    // Read in a snapshot:

    if (read_from_file) {
	ifstream s(infile);
	if (!s) {
	    cerr << "data file " << infile << " not found." << endl;
	    exit(1);
	}
	cerr << "input file is " << infile << endl;
	b = get_hdyn(s);
    } else {
	cerr << "reading input from stdin" << endl;
	b = get_hdyn();
    }

    if (b == NULL) err_exit("Can't read input snapshot");

    // NOTE: get_hdyn() reads in the hdyn tree structure using get_node()
    //	     and *also* initializes parameters for unperturbed and slow
    //	     binary motion.  (These parameters require knowledge of the
    //	     tree structure, which is not known until after get_node() is
    //	     complete.

    int nbody = b->n_leaves();
    if (nbody <= 0) err_exit("kira: N <= 0!");

    // ======================================================================
    //
    // New (from Steve, 6/03): the standard external representation (always
    // used in the dyn code) now is for the root node to be at the center of
    // mass of the top-level nodes, and for all center of mass information
    // to reside in root.  This description is consistent with the tree
    // structure at all other levels.  However, within kira it is much more
    // convenient to place the root node at the origin and drop the CM
    // requirement for the daughters (e.g. makes dealing with external fields
    // much simpler).  The member functions offset_com() and reset_com()
    // perform the transformations between these two equivalent descriptions
    // of the system.  We transform root to the origin here.  We must always
    // make the reverse transformation before writing data to disk to maintain
    // the external representation..

    cerr << "offsetting system to place root node at (0,0,0)" << endl;
    b->offset_com();

    // ======================================================================

    // Other preliminaries:

    b->set_kira_counters(new kira_counters);
    b->set_kira_options(new kira_options);
    b->set_kira_diag(new kira_diag);

    check_kira_init(b);			// allow changes to kira defaults

    b->log_history(argc, argv);

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string)) != -1) {
	switch (c) {
	    case 'a':	eta = atof(poptarg);
			eta_flag = true;
			break;
	    case 'A':	alt_flag = true;
	    		cerr << endl
			     << "additional log output enabled"
			     << endl;
			break;
	    case 'b':	if (poptarg)
			   long_binary_out = Starlab::max(0, atoi(poptarg));
			else
			   long_binary_out = 0;
			break;
	    case 'B':   S_flag = true;	    // does B_flag = true make sense
	    		B_flag = true;	    // if S_flag = false?
		        break;
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'C':	b->get_kira_options()->grape_max_cpu = atof(poptarg);
			break;
	    case 'd':	dt_log = atof(poptarg);
	    		if (dt_log < 0) dt_log = pow(2.0, dt_log);
			log_flag = true;
			break;
	    case 'D':	if (streq(poptarg, "all") || streq(poptarg, "full")) {

			    // Turn on full_dump mode.

			    dt_snap = -1;
			    set_write_unformatted(false);

			} else if (*poptarg == 'x') {

			    // Turn on full_dump mode and set its frequency.
			    // Messy...

			    dt_snap = -atoi(poptarg+1);
			    set_write_unformatted(false);

			} else if (*poptarg == 'X') {

			    // Turn on full_dump mode, set its frequency, and
			    // specify unformatted output.
			    // Messier...

			    dt_snap = -atoi(poptarg+1);
			    set_write_unformatted(true);

			} else if (*poptarg == 'B' || *poptarg == 'b') {

			    // Turn on full_dump mode for binary tracking only.

			    dt_snap = 0;
			    binary_zero = true;

			    if (*poptarg == 'b')
				set_write_unformatted(false);
			    else
				set_write_unformatted(true);

			    cerr << endl
				 << "binary tracking on"
				 << endl;

			} else{
			    dt_snap = atof(poptarg);
			    if (dt_snap < 0) dt_snap = pow(2.0, dt_snap);
			}
			snap_flag = true;
			break;
	    case 'e':	eps = atof(poptarg);
			eps_flag = true;
			break;
	    case 'E':	if (atoi(poptarg)) exact = true;
			break;
	    case 'f':	d_min = atof(poptarg);
			d_min_flag = true;
			break;
	    case 'F':	// err_exit("kira: -F option removed: use add_tidal");
			if (poptarg)
			    friction_beta = atof(poptarg);
			else
			    friction_beta = 1;
	    		if (friction_beta < 0) friction_beta = 0;
			break;
	    case 'g':	lag_factor = atof(poptarg);
			lag_flag = true;
			break;
	    case 'G':	G_flag = true;
			input_stripping_radius = atof(poptarg);
			break;
	    case 'h':	dt_sstar = atof(poptarg);
	    		if (dt_sstar < 0) dt_sstar = pow(2.0, dt_sstar);
			sstar_flag = true;
			break;
	    case 'i':	ignore_internal = !ignore_internal;
			break;
	    case 'I':	dt_reinit = atof(poptarg);
	    		if (dt_reinit < 0) dt_reinit = pow(2.0, dt_reinit);
			reinit_flag = true;
			break;
	    case 'J':	err_exit("kira: -J option removed: use add_tidal");
			break;
	    case 'k':	gamma = atof(poptarg);
			gamma_flag = true;
			break;
	    case 'K':	max_slow = atoi(poptarg);
			max_slow_flag = true;
			break;
	    case 'L':	cpu_time_limit = atof(poptarg);
			break;
	    case 'M':	err_exit("kira: -M option removed: use add_star");
			break;
	    case 'n':	n_stop = atoi(poptarg);
			break;
	    case 'N':	b->get_kira_diag()->check_heartbeat = true;
			b->get_kira_diag()->n_check_heartbeat = atoi(poptarg);
			break;
	    case 'o':	allow_kira_override = false;
	    		break;
	    case 'O':	save_snap_at_log = true;
	    		strcpy(snap_save_file, poptarg);
	    		break;
            case 'q':	q_vir = atof(poptarg);
			break;
            case 'Q':	err_exit("kira: -Q option removed: use add_tidal");
	       		break;
	    case 'r':	// err_exit("kira: -r option removed: use scale");

			// Not desirable to use this option -- better to
			// set with scale -- but retain just in case.

			r_virial_set = true;
			input_r_virial = atof(poptarg);	    
			break;
	    case 'R':	// err_exit("kira: -R option removed: use add_star");
			break;
            case 's':	s_flag = true;
			input_seed = atoi(poptarg);
			break;
            case 'S':	S_flag = true;
	       		break;
	    case 't':	delta_t = atof(poptarg);
			break;
            case 'T':	err_exit("kira: -T option removed: use add_star");
			break;
	    case 'u':	toggle_unperturbed(b, 1);
			break;
	    case 'U':	toggle_unperturbed(b, 0);
			break;
	    case 'v':	verbose = !verbose;
			break;

	    case 'V':	print_recalculated_energies(b);
			exit(0);

	    case 'W':	dt_fulldump = atof(poptarg);
	    		if (dt_fulldump < 0)
			    dt_fulldump = pow(2.0, dt_fulldump);
			fulldump_flag = true;
			break;
	    case 'x':	b->get_kira_options()->print_xreal
				= !b->get_kira_options()->print_xreal;
			break;
	    case 'X':	dt_esc = atof(poptarg);
	    		if (dt_esc < 0) dt_esc = pow(2.0, dt_esc);
			esc_flag = true;
			break;
	    case 'y':	sec_flag = true;
	                sec = atof(poptarg);
			break;
	    case 'z':	smc_flag = true;
	                smc = atof(poptarg);
			break;
	    case 'Z':	stc_flag = true;
	                stc = atof(poptarg);
			break;
	    case '*':	b->get_kira_diag()->kira_runtime_flag = atoi(poptarg);
			break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }

    if (c_flag)
	b->log_comment(comment);

    // Perform one-time correction/initialization of various system components.

    correct_multiples(b, verbose);
    initialize_index(b, verbose);

    set_friction_beta(friction_beta);
    cerr << "dynamical friction beta = " << friction_beta << endl;

    check_set_ignore_internal(b, verbose);
    bool snap_ignore_internal = b->get_ignore_internal();

    if (ignore_internal || snap_ignore_internal) {
	if (!snap_ignore_internal) {
	    cerr << "command-line \"-i\" flag forces ignore_internal mode"
		 << endl;
	    putiq(b->get_log_story(), "ignore_internal", 1);

	} else

	    cerr << endl
		 << "ignoring internal forces "
		 << "-- external field only"
		 << endl;

	b->set_ignore_internal(true);
    }

    //----------------------------------------------------------------------

    // Some consistency checks (review from time to time to check if
    // these are really what we want)...  Disable with "-o" on the
    // command line.

    // 1. Turn on stripping if it was previously enabled (performed
    //    below, in the "stripping" section).

    // 2. Turn on stellar/binary evolution if it was turned on in the
    //    previous run:

    bool need_skip = true;	// formatting!!

    if (check_kira_flag(b, "kira_evolve_stars") && !S_flag)
	if (check_allowed(allow_kira_override,
			  "stellar evolution",
			  verbose, need_skip))
	    S_flag = true;

    if (check_kira_flag(b, "kira_evolve_binaries") && !B_flag)
	if (check_allowed(allow_kira_override,
			  "binary evolution",
			  verbose, need_skip))
	    B_flag = true;

    // 3. Turn on unperturbed binaries/multiples if they were previously
    //    enabled.

    if (check_kira_flag(b, "kira_allow_unperturbed")
	&& !b->get_kira_options()->allow_unperturbed)
	if (check_allowed(allow_kira_override,
			  "unperturbed binary motion",
			  verbose, need_skip))
	    set_allow_unperturbed(b);

    if (check_kira_flag(b, "kira_allow_multiples")
	&& !b->get_kira_options()->allow_multiples)
	if (check_allowed(allow_kira_override,
			  "unperturbed multiple motion",
			  verbose, need_skip))
	    set_allow_multiples(b);

    putiq(b->get_log_story(), "kira_allow_unperturbed",
	  b->get_kira_options()->allow_unperturbed);
    putiq(b->get_log_story(), "kira_allow_multiples",
	  b->get_kira_options()->allow_multiples);

    if (verbose) {
	if (need_skip) cerr << endl;
	print_unperturbed_options(b);
    }

    //----------------------------------------------------------------------

    // May be somewhat redundant...

    if (B_flag) S_flag = true;
    if (!S_flag) B_flag = false;

    if (!S_flag) dt_sstar = VERY_LARGE_NUMBER;

    if (S_flag && dt_sstar == VERY_LARGE_NUMBER)
	dt_sstar = DEFAULT_DT_SSTAR;

    b->set_use_sstar(S_flag);
    b->set_use_dstar(B_flag);

    // Need random numbers if stellar evolution is enabled.

    if (S_flag) {
        if (s_flag == FALSE)
	    input_seed = 0;			// default
        actual_seed = srandinter(input_seed);
        sprintf(seedlog,
		"       random number generator seed = %d",actual_seed);
        b->log_comment(seedlog);
	cerr << "initial random seed = " << actual_seed << endl;
    }

    // Establish defaults for time scales:

    if (!snap_flag) dt_snap = delta_t;		// snap output at end
    if (dt_snap == 0 && !binary_zero)
	dt_snap = VERY_LARGE_NUMBER;		// suppress snap output

    if (!esc_flag) dt_esc = dt_reinit;		// check for escapers at each
						//     reinitialization
    if (!fulldump_flag) dt_fulldump = dt_esc;	// full dump on escaper check

    //----------------------------------------------------------------------

    // Update static runtime member data:

    set_runtime_params(b, verbose,
		       eta, eta_flag,
		       eps, eps_flag,
		       d_min, d_min_flag,
		       lag_factor, lag_flag,
		       gamma, gamma_flag,
		       max_slow, max_slow_flag);

    // Get time steps in a similar manner.  Note that the use of the log
    // story effectively allows the user's previous settings to establish
    // a "preference" for the current run, overriding the defaults.  It
    // remains to be seen whether this is really the most convenient
    // choice.						    (Steve, 7/99)

    choose_param(b, verbose, dt_log, log_flag, "dt_log");
    choose_param(b, verbose, dt_sstar, sstar_flag, "dt_sstar");
    choose_param(b, verbose, dt_esc, esc_flag, "dt_esc");
    choose_param(b, verbose, dt_reinit, reinit_flag, "dt_reinit");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Set up merger/encounter criteria, with same precedences as for
    // runtime and timestep parameters.

    choose_param(b, verbose, sec, sec_flag,
		 "stellar_encounter_criterion", true);
    choose_param(b, verbose, smc, smc_flag,
		 "stellar_merger_criterion", true);
    choose_param(b, verbose, stc, stc_flag,
		 "stellar_tidal_capture_criterion", true);

    // (Note that the criterion for merger currently implemented
    // in check_merge_node is:
    //
    //	distance_squared - sum_of_radii_sq
    //		<= (stellar_merger_criterion_sq-1) * sum_of_radii_sq
    //
    // so "-z 1" implies merger on contact.)

    if (stc < smc) stc = smc;

    b->set_stellar_encounter_criterion_sq(pow(sec, 2));
    b->set_stellar_merger_criterion_sq(pow(smc, 2));
    b->set_stellar_capture_criterion_sq(pow(stc, 2));

    //----------------------------------------------------------------------

    // Determine the initial virial radius and initial mass (both in code
    // units), for use in scaling.

    // Note: the parameters set here must give correct results both on
    //       initialization AND on restart...

    if (verbose) cerr << endl;

    real initial_mass = get_initial_mass(b, verbose);

    real initial_r_virial = get_initial_virial_radius(b, verbose,
						      r_virial_set,
						      input_r_virial);

    // NOTE: important that initial_r_virial be set, as several other
    // startup functions depend on it.  Reinstated the '-r' option, in
    // case no previous function has determined the virial radius.
    //						(Steve, 10/01)

    //----------------------------------------------------------------------

    // Check and reset the total mass, just in case...
    // Don't check individual nodes for now (Steve, 4/11/01).

    check_total_mass(b);

    //----------------------------------------------------------------------

    // Check for external fields, including tidal fields.
    // Note that there are NO command-line options -- fields can be
    // specified *only* via the initial snapshot.

    check_set_external(b, verbose);

    if (friction_beta > 0 && !b->get_pl())
	cerr << "warning: dynamical friction, but "
	     << "no external power-law field set."
	     << endl;

    //----------------------------------------------------------------------

    // Stripping of outlying stars (specify with "-G" on the command line).

    real scaled_stripping_radius = 0;	// stripping radius for unit mass
    					// (0 ==> no stripping)

    if (check_kira_flag(b, "kira_remove_escapers") && !G_flag)
	if (check_allowed(allow_kira_override,
			  "escaper removal",
			  verbose, need_skip))
	    G_flag = true;

    if (G_flag) {

	if (verbose) cerr << endl;

	real initial_r_jacobi
	    = get_initial_jacobi_radius(b, initial_r_virial);

	real scaled_stripping_radius
	    = get_scaled_stripping_radius(b, verbose,
					  input_stripping_radius,
					  initial_r_jacobi,
					  initial_r_virial,
					  initial_mass);

	if (scaled_stripping_radius <= 0)
	    err_exit("Unable to determine stripping radius");

	PRL(scaled_stripping_radius);
	cerr << "current stripping radius = "
	     << scaled_stripping_radius * pow(total_mass(b), 1.0/3.0)
	     << endl;

	b->set_scaled_stripping_radius(scaled_stripping_radius);

    } else

	if (verbose) cerr << endl << "no escaper removal" << endl;

    // Save information on whether or not escapers are removed.

    putiq(b->get_log_story(), "kira_remove_escapers", G_flag);

    //----------------------------------------------------------------------

    // (Re)initialize stellar evolution.

    // On restart, pick up scalings from the root star story (if any),
    // impose new scalings (if any), and (re)initialize star parts.

    if (S_flag) {	// Note that B_flag ==> S_flag (Steve 9/19/97)

	if (verbose) cerr << endl;
	if (b->get_starbase() == NULL)
	    err_exit("kira: S_flag and no starbase!");

	// See if any scaling information already exists.

	bool scales_from_snapshot
	    = b->get_starbase()->get_stellar_evolution_scaling();

	// get_stellar_evolution_scaling() will read scaling factors
	// from the input stream if possible.  Return value is true iff
	// all scaling factors are now known.

	if (scales_from_snapshot) {

	    // All scales were specified in the input snapshot.

#if 0
	    // Horrible #ifdef mess here is to accommodate transition
	    // from old to new style in handling stellar properties.

	    // Old:

	    // See if command-line input will override them.

	    if (R_flag || M_flag || T_flag) {

		if (verbose)
		    cerr << "overwriting (some) scale factors from"
			 << " input snapshot"
			 << endl;

		scales_from_snapshot = false;

	    } else {

		if (verbose) {
		    cerr << "scale factors taken from input snapshot"
			 << endl;
		    b->get_starbase()->print_stellar_evolution_scaling(cerr);
		}

	    }

	} else {

	    // (Re)set some or all physical scales.

	    get_physical_scales(b, verbose,
				initial_mass,
				initial_r_virial,
				M_flag, m_tot,
				R_flag, r_vir,
				T_flag, t_vir,
				q_vir,
				nbody);
	}
#else
	    // New:

	    if (verbose) {
		cerr << "scale factors taken from input snapshot"
		     << endl;
		b->get_starbase()->print_stellar_evolution_scaling(cerr);
	    }

	} else {

	    if (verbose)
		cerr << "stellar scaling unavailable -- "
		     << "suppressing stellar evolution." << endl;

	    B_flag = S_flag = false;
	    b->set_use_sstar(S_flag);
	    b->set_use_dstar(B_flag);
	}
#endif

        if (S_flag) {

	    // Add information on the physical initial mass to the
	    // root log story.

	    putrq(b->get_log_story(), "physical_initial_mass",
		  b->get_starbase()->conv_m_dyn_to_star(initial_mass));

	    // Add star parts to nodes.

	    addstar(b,                         // Note that T_start and
		    T_start,                   // Main_Sequence are
		    Main_Sequence,             // defaults, ignored if a
		    true);                     // star story already exists.

	    // Command line specified sec in solar radii.  Convert it to
	    // N-body units for use by check_merge_nodes (which uses d_nn_sq).

	    sec = b->get_starbase()->conv_r_star_to_dyn(sec);
	    b->set_stellar_encounter_criterion_sq(pow(sec, 2));

	    if (verbose) {
		cerr << endl;
		cerr << "stellar_encounter_criterion = "
		     <<  sqrt(b->get_stellar_encounter_criterion_sq()) << endl;
		cerr << "stellar_merger_criterion = "
		     <<  sqrt(b->get_stellar_merger_criterion_sq()) << endl;
		cerr << "stellar_capture_criterion = "
		     <<  sqrt(b->get_stellar_capture_criterion_sq()) << endl;
	    }

	    // Print a diagnostic on tidal parameters in the disk case:

	    if (b->get_tidal_field() == 3) {

		real initial_r_jacobi
		    = get_initial_jacobi_radius(b, initial_r_virial);

		test_tidal_params(b, verbose,
				  initial_r_jacobi,
				  initial_r_virial,
				  initial_mass);
	    }
	}
    }

    // Save information on whether or not stellar evolution is enabled.

    putiq(b->get_log_story(), "kira_evolve_stars", S_flag);
    putiq(b->get_log_story(), "kira_evolve_binaries", B_flag);

    //----------------------------------------------------------------------

    return true;
}


