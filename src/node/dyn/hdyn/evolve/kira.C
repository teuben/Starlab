
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                       
//=======================================================//              /|\ ~


// kira:  This file contains only the main program, evolve_system, and
//	  helper functions.
//
//	  Externally visible functions when compiled as a module:
//
//		void kira
//		void kira_finalize


//#define DUMP_DATA 1	// uncomment to allow detailed TMP_DUMP output


//// Hermite N-body integrator with evolving hierarchical tree structure,
//// stellar and binary evolution, and an arbitrary external field.  The
//// program reads a snapshot from standard input and writes snapshot(s)
//// to standard output.  Periodic log output is sent to standard error.
////
//// Basic options are listed below.  Default values for a new simulation
//// are indicated in square brackets.  For restart (continuation of a
//// previous kira calculation), the defaults for many options [*] are
//// determined from the input snapshot, making it possible to continue
//// a run without having to re-specify all the command-line parameters
//// previously used.  (These new defaults may still be overridden on
//// the command line, of course.)  Some options typically having to do
//// with initial setup may be overridden by data from the input
//// snapshot (if present), as noted below.  Kira may also turn *on*
//// some options (the B, G, S, and u settings) if they were turned
//// on in the previous run, but are not specified on the current
//// command line.  To prevent this, use the "-o" switch.
////
//// The first page of log output gives detailed information on all
//// parameter settings adopted and any modifications made during
//// initialization.  In case of doubt, read the log file!
////
//// Usage: kira [OPTIONS] < infile > output
////
//// Options:
////              -0    force non-GRAPE operation, if relevant
////              -1    suppress density calculation [compute with GRAPE]
////              -2    enable DMA GRAPE access (experimental) [no DMA]
////              -3    enable special treatment of isolated multiples [no]
////              -a    specify accuracy parameter [0.1][*]
////              -A    enable "alternate" output [off]
////              -b    specify frequency of full binary output, in (integer)
////                    units of the log output interval
////                    [10; no argument or 0 ==> no output]
////              -B    turn on binary evolution [off][*]
////              -c    include comment [none]
////              -C    specify GRAPE release interval, in seconds [15]
////              -d    specify log output interval [1][*]
////              -D    specify snapshot interval [end of run]. 
////                    Special values:
////                                xN: formatted full dump, frequency N;
////                                XN: unformatted full dump, frequency N;
////                                full/all: same as x1;
////                                b: track binary changes only (formatted);
////                                B: track binary changes (unformatted);
////              -e    specify softening length [0][*]
////              -E    use exact calculation [false]
////              -f    turn on/off internal dynamical friction on stars [0][*]
////              -F    turn on external dynamical friction on the cluster,
////                    and optionally specify a scaling coefficient [none]
////              -g    specify hysteresis factor [2.5][*]
////              -G    specify initial stripping radius [none][*]
////              -h    specify stellar-evolution time step [0.015625 = 1/64][*]
////              -i    ignore all internal forces (i.e. external only) [false]
////              -I    specify (re)initialization timescale [1][*]
////              -k    specify perturbation factor [1.e-7][*]
////              -K    specify log2(maximum slowdown factor) (integer): [0][*]
////              -l    specify close-encounter distance [0.25 --> 0.25/N][*]
////                    [option renamed from -f, 7/04]
////              -L    specify CPU time limit, in seconds [none]
////              -n    stop at specified number of particles [5]
////              -N    specify frequency of CPU check output [50000]
////              -o    prevent kira from overriding some settings (BGSu)
////                    based on input snapshot data [allow]
////              -O    save (and overwrite) extra snapshot at each output [no]
////              -q    specify initial virial ratio [0.5]
////              -r    specify initial virial radius (may not be
////                    specified in the input snap) [no]
////              -R    specify snapshot file for (re)start [none: use stdin]
////              -s    specify random seed [take from system clock]
////              -S    turn on stellar evolution [off][*]
////              -t    specify time span of calculation [10]
////              -T    enable experimental threading and specify n_threads [0]
////              -u    toggle unperturbed multiple motion [disabled][*]
////              -U    toggle all unperturbed motion [enabled][*]
////              -v    toggle "verbose" mode [on]
////              -W    specify full-dump (worldline) timescale [1]
////              -x    toggle output of extended-precision time [on]
////              -X    specify escaper removal timescale [reinit][*]
////              -y    specify stellar encounter criterion
////                    [0 N-body units or solar radii][*]
////              -z    specify stellar merger criterion [0 stellar radii][*]
////              -Z    specify stellar tidal dissipation criterion
////                    [0 stellar radii][*]
////
//// As a convenient shorthand, any "dt" interval specified less than zero
//// is interpreted as a power of 2, i.e. "-d -3" sets dt_log = 0.125.  In
//// the case of dt_snap, this will also cause snapshot data to be written
//// immediately on restart (usually we wait until time dt_snap).
////
//// Written by J. Makino, S. McMillan, S. Portegies Zwart, and P. Hut, .
////
//// Report bugs to starlab@sns.ias.edu.

// Level-2 help:

//++
//++ Some run-time parameters:
//++
//++ The initial virial radius is read from the input snapshot.  If no
//++ initial virial radius is found there, the value specified by the
//++ "-r" command-line option is used.  If no "-r" option is set, the
//++ value is computed fom the energy.
//++
//++ The initial mass is read from the input snapshot.  If no initial
//++ mass is found there, the mass is computed if t = 0; otherwise, a
//++ default value of 1 is assumed.  The initial mass is saved in the
//++ output snapshot for future use.
//++
//++ If a tidal field is specified ("-Q"), the initial Jacobi radius
//++ is set from the "initial_rtidal_over_rvirial" entry in the input
//++ file, if such an entry exists.  This value may be overridden by
//++ setting the Jacobi radius on the command line ("-J").  If specified,
//++ the command-line value scales the initial tidal radius (if found),
//++ or the virial radius (otherwise).  The initial Jacobi radius is
//++ written to the output snapshot for future use.  If a
//++ "kira-written" initial Jacobi radius is found in the input file,
//++ it is used regardless of the other settings.  Thus, in a typical
//++ situation, the Jacobi radius computed at t = 0 will be used in
//++ all subsequent restarts.
//++
//++ A stripping radius may be specified with the "-G" option.  If a
//++ command-line value specified, it scales the Jacobi radius if one
//++ is known, or the virial radius otherwise.  The scaled stripping
//++ radius (normalized to unit mass) is written to the output
//++ snapshot for future use.  If a scaled stripping radius is found
//++ in the input file, it is used regardless of the other settings.
//++ As with the initial Jacobi radius, in a typical situation, the
//++ scaled stripping radius computed at t = 0 will be used in all
//++ subsequent restarts.
//++
//++ If the input snapshot does not contain physical scales (and if
//++ stellar evolution/interactions are enabled), they may be set on
//++ the command line.  The "-M" option specifies the total system
//++ mass, in solar units.  The "-R" option specifies the value of
//++ the virial radius, in parsecs.  The "-T" option specifies the
//++ (Heggie & Mathieu) system time unit, in Myr.
//++
//++ NOTE from Steve (7/01): the following options have been retired:
//++
//++      -F    specify tidal field type (0 = none, 1 = point-mass,
//++                                      2 = isothermal, 3 = disk,
//++                                      4 = custom)
//++                          [0 if -Q not set; 1 otherwise][*]
//++
//++      -J    specify Jacobi radius [none][*]
//++
//++      -Q    use tidal field [none][*]
//++
//++      -M    specify initial total mass, in solar masses
//++                          [none; take from input snap if present]
//++      -R    specify initial virial radius, in pc
//++                          [none; take from input snap if present]
//++      -T    specify initial virial time scale, in Myr
//++                          [none; take from input snap if present]
//++
//++      -r    specify initial virial radius in code units
//++                          [1; take from input snap if present]

//	J. Makino, S. McMillan, P. Hut, S. Portegies Zwart	12/92-6/04
//
//	Major top-level reorganization to remove compile-time
//	GRAPE selection, S. McMillan				6/04



#include "hdyn.h"

#ifdef USEMPI
#    include <mpi.h>
#    include "kira_mpi.h"
#endif

#ifndef TOOLBOX

#include "star/dstar_to_kira.h"
#include "kira_timing.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira
#   undef T_DEBUG
#endif

#define INITIAL_STEP_LIMIT  0.0625	// arbitrary limit on the first step

static hdyn **next_nodes = NULL;	// *** MAKE THIS GLOBAL HDYN ***
static kira_counters kc_prev; 



local void full_reinitialize(hdyn* b, xreal t, bool verbose,
			     bool init = false)
{
    real cpu_0 = cpu_time();

    b->set_system_time(t);
    b->set_time(t);
    initialize_system_phase1(b, t);
    initialize_system_phase2(b, 3, 0);			// "0" here means
							// timesteps are only
							// set if zero    
    // Reset and save d_min_sq().

    int n = 0;
    real mass = 0;
    real pot = 0;
    for_all_daughters(hdyn, b, bb) {
	n++;
	mass += bb->get_mass();
	pot += bb->get_mass()*bb->get_pot();		// total potential
    }

    if (n <= 1) return;		// sometimes handy to follow just 1 particle

    cerr << endl << "reinitializing system at time " << t << endl;

    // Only want the internal potential, which is counted twice in pot.

    pot -= get_external_pot(b);
    pot /= 2;

    int p = cerr.precision(INT_PRECISION);

    PRC(mass); PRC(pot);
    real r_virial = -0.5*mass*mass/pot;
    PRC(r_virial); PRL(n);

    cerr << "old "; PRC(b->get_d_min_sq());
    real d_min_sq = square(b->get_d_min_fac()*r_virial/n);
    b->set_d_min_sq(d_min_sq);
    cerr << "new "; PRL(b->get_d_min_sq());

    cerr.precision(p);

    putrq(b->get_log_story(), "kira_d_min_sq", d_min_sq);

    b->reconstruct_perturbed_list();

    // Quietly reinitialize all kepler structures.
    // (Repeats actions taken by get_hdyn() on startup.)

    b->initialize_unperturbed();

    // Note from Steve, 7/03: the explanation for the following call to
    // recompute_unperturbed_steps() is that some external agent may
    // have synchronized some unperturbed binaries, in which case
    // the step may no longer end at the correct phase of the orbit.
    // However, this is usually not the case, and if applied routinely,
    // this action will destroy reproducibility.  Could just omit it,
    // except in special cases (command-line flag? -- in which case
    // this should be done in main).  For now, we retain the call, but
    // recompute only if the binary time is the same as the system time
    // (normally, would expect the components to lag the system time).

    if (init) b->recompute_unperturbed_steps();

//#ifdef USEMPI
    recompute_MPI_id_array(b);
    if (!check_MPI_id_array(b)) exit(1);
//#endif

    if (verbose) {
	cerr << "CPU time for reinitialization = "
	     << cpu_time() - cpu_0 << endl;
    }
}

#define N_CHECK 3

local bool check_sync(hdyn* b, int max_count = N_CHECK)
{
    bool need_new_list = false;
    int count = 0;

    for_all_nodes(hdyn, b, bb) {

#if 0
	if (bb->name_is("2009")) {
	    cerr << bb->format_label() << " time = "
		 << bb->get_time() << " (";
	    xprint(bb->get_time(), cerr, false);
	    cerr << ")" << endl;
	    cerr << "system time = "
		 << b->get_system_time() << " (";
	    xprint(b->get_system_time(), cerr, false);
	    cerr << ")" << endl;
	}
#endif

	if (bb->get_time() != b->get_system_time()
	     && bb->get_kepler() == NULL) {

	    if (++count <= max_count) {
	        cerr << "check_sync warning:  node "
		     << bb->format_label() << " not synchronized" << endl
		     << "                     system time = "
		     << b->get_system_time() << " (";
		xprint(b->get_system_time(), cerr, false); cerr << ")" << endl;
		cerr << "                     node time   = "
		     << bb->get_time() << " (";
		xprint(bb->get_time(), cerr, false); cerr << ")" << endl;
		need_new_list = true;
	    }
	}
    }

    if (count > max_count) {
	if (max_count > 0)
	    cerr << ".\n.\n(plus "
		 << count-max_count << " more warnings)" << endl;
	else
	    cerr << "*** " << count << " unsynchronized nodes" << endl;
    }

    return need_new_list;
}

local void backward_step_exit(hdyn* b, real ttmp, real t,
			      hdyn** next_nodes, int n_next)
{

    cerr << "evolve_system: time went backwards!\n";

    cerr.precision(HIGH_PRECISION);
    PRC(ttmp); PRC(b->get_system_time()); PRL(t);

    for (int i = 0; i < n_next; i++) {
	cerr << endl;
	PRC(i); next_nodes[i]->pretty_print_node(cerr);
	cerr << " "<< next_nodes[i]->get_time()<< " ";
	cerr << next_nodes[i]->get_next_time()<<endl;
	pp3(next_nodes[i]->get_top_level_node(), cerr);
    }
    put_node(b, cerr, b->get_kira_options()->print_xreal);
    exit(0);
}

local void check_binary_scheduling(hdyn* b,
				   hdyn** next_nodes,
				   int n_next)
{
    // Force the elder of two binary sisters to be integrated.
    // Is this function necessary?  (Steve, 7/00)

    // Old code (inefficient?):

//      for (int i = 0; i < n_next; i++) {
//  	hdyn *bi = next_nodes[i];
//
//  	if (bi->is_low_level_node()) {
//
//  	    hdyn *od = bi->get_parent()
//  	                 ->get_oldest_daughter();
//
//  	    if (od->get_timestep() != bi->get_timestep()) {
//  		// cerr << "timesteps of sisters are different !! \n";
//  		od->set_timestep(bi->get_timestep());
//  	    }
//
//  	    if (od->get_time() != bi->get_time()) {
//  		od->set_time(bi->get_time());
//  	    }
//
//  	    next_nodes[i] = bi->get_parent()->get_oldest_daughter();
//  	}
//      }

    // Should be better:

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	if (bi->is_low_level_node()) {

	    hdyn *od = bi->get_elder_sister();
	    if (od) {

		// We are a low-level younger sister.  Shouldn't happen!

		if (od->get_timestep() != bi->get_timestep())
		    // cerr << "timesteps of sisters are different !! \n";
		    od->set_timestep(bi->get_timestep());

		if (od->get_time() != bi->get_time())
		    od->set_time(bi->get_time());

		next_nodes[i] = od;
	    }
	}
    }
}

local inline void update_step(real t, real& ts, real dts)
{
    while (ts < t)
	ts += dts;
}



local real internal_potential(hdyn* bb)
{
    real epot, ekin, etot;
    calculate_energies(bb, bb->get_eps2(), epot, ekin, etot);
    return epot;
}

local void print_energy_from_pot(hdyn* b)
{
    // Attempt faster energy computation using existing pot data.

    // NOTE: binary components' pots only include perturbers, *not* the
    // entire system, so we can't use these directly.  However, the pot
    // of the top-level node should be correct, apart from the internal
    // energy.  There will still be a tidal error in the case unperturbed
    // binaries...

    // ***** NOT WORKING YET.  There appears to be a residual error in
    // ***** the energy derived from pots.  Not clear if this is a bug
    // ***** in the computation of pot for a binary CM or if the error
    // ***** is inherent in the approximations made by kira...
    //
    //						Steve (6/99)

    return;

    real top_pot = 0, int_pot = 0;

    for_all_daughters(hdyn, b, bb) {
	top_pot += bb->get_mass()*bb->get_pot();
	if (bb->is_parent())
	    int_pot += internal_potential(bb);
    }

    real kin = 0;

    for_all_leaves(hdyn, b, bb) {
	vec vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	kin += bb->get_mass()*square(vel);
    }

    top_pot *= 0.5;
    kin *= 0.5;

    int p = cerr.precision(INT_PRECISION);
    PRC(top_pot), PRC(int_pot), PRL(kin);
    PRL(get_external_pot(b));
    cerr << "energy_from_pot = " << kin+top_pot+int_pot+0.5*get_external_pot(b)
	 << endl;
    cerr.precision(p);
}

local void kira_alt_output(hdyn *b)
{
    ofstream alt("alt_output", ios::app|ios::out);

    vec cmpos, cmvel;
    get_std_center(b, cmpos, cmvel);	        // use CM for now (note: doesn't
    					        // use pred quantities -- should
						// create a "pred" version...)

    cmpos -= b->get_pos();			// std_center quantities
    cmvel -= b->get_vel();			// include the root node

    for_all_daughters(hdyn, b, bb) {		// top-level nodes only

	real sys_t = bb->get_system_time();
	real m = bb->get_mass();
	vec dr = bb->get_pred_pos() - cmpos;
	vec dv = bb->get_pred_vel() - cmvel;
	real r2 = square(dr);
	real v2 = square(dv);
	real vr2 = square(dr*dv)/r2;
	
	alt << sys_t << " "		       	// wasteful -- fix soon!
	    << m << " " << r2 << " " << v2 << " " << vr2
	    << endl;
    }

    alt.close();
}



#define TMP_SUPPRESS
#undef  TMP_SUPPRESS

// evolve_system:  Main integration loop.
//
//		   Called only from kira() and local to this file.

local void evolve_system(hdyn * b,		// hdyn array
			 real delta_t,		// time span of the integration
			 real dt_log,		// time step of the integration
			 int long_binary_out,
			 real dt_snap,		// snapshot output interval
			 real dt_sstar,		// timestep for single
						//     stellar evolution
			 real dt_esc,		// escaper removal
			 real dt_reinit,	// reinitialization interval
			 real dt_fulldump,	// full dump interval
			 bool exact,		// exact force calculation
			 real cpu_time_limit,
			 bool verbose,
			 bool snap_init,
			 bool save_snap_at_log,	// save snap at log output
			 char* snap_save_file,	// filename to save in
			 int n_stop,		// when to stop (# stars)
			 bool alt_flag)		// enable alternative output

{
    // Modified order in which output/snapshots/reinitialization/etc. are
    // performed -- Steve, 7/01

    // Carried through option to ignore all internal forces
    // and added snapshot "reflection" option		    -- Steve, 12/01.

    // Initialization:

    clean_up_files();
    cpu_init();

    real r_reflect = -1;
    if (find_qmatch(b->get_log_story(), "r_reflect")) {
	r_reflect = getrq(b->get_log_story(), "r_reflect");
	cerr << endl
	     << "*** Reflecting boundary at radius " << r_reflect << " ***"
	     << endl;
    }

    kira_counters *kc = b->get_kira_counters();

#ifdef CPU_COUNTERS
    real cpu = cpu_time(), cpu_prev;
#endif

    initialize_counters_from_log(b);
    kc_prev = *kc;					// (make a local copy)

    kira_options *ko = b->get_kira_options();
    kira_diag *kd = b->get_kira_diag();

    real count = 0, steps = 0;
    real count_top_level = 0, steps_top_level = 0;
    int grape_steps = 0;
    int snaps = 0;

    next_nodes = new hdyn *[2 * b->n_daughters()];	// Overkill?
    bool *next_flag = new bool[2 * b->n_daughters()];
    real *last_dt = new real[2 * b->n_daughters()];

    set_n_top_level(b);

    // Establish limits on time steps.  Note that, other than the static
    // step_limit data, the "dt" data never propagate beyond this function,
    // so there is no particular reason to make them part of the hdyn class.

    real dt_max = Starlab::min(dt_log, dt_sstar); // system will be synchronized
						  // for stellar evolution
    dt_max = Starlab::min(dt_max, dt_reinit);
    dt_max = Starlab::min(dt_max, dt_fulldump);
    b->set_unpert_step_limit(dt_max);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Convenient for now to pass request for full dump and also the
    // dump interval via dt_snap...
    //
    // Full_dump options are becoming more complex.  Now, in addition
    // to specifying the full_dump frequency, we can also opt to produce
    // limited output to follow binary changes, and specify the length
    // of the worldline interval.  Worldline interval (time between
    // complete system dumps) is managed separately; other options are
    // managed via the full_dump variable (now int rather than bool):
    //
    //		full_dump  =  0			// option off
    //			      1			// dump all particles
    //			      2			// track binary changes only
    //
    // Details if how these variables are used, specified and passed around
    // the system remain to be finalized.			(Steve, 9/01)
    //
    // For now, we track all internal binary changes (makes the bookeeping and
    // reconstruction simpler -- read_bundle() should work as is).  Later, we
    // may decide to track only changes involvng the top-level nodes, in which
    // case the reconstruction code will have to become more sophisticated.

    int full_dump = 0;
    real n_dump = 1.0;				// up to at least 30 seems OK!

    bool first_full_dump = false;

    if (dt_snap <= 0) {

	// Turn on full_dump mode, at frequency -dt_snap, or binary dump
	// mode.  Turn off regular snapshot output.  Complete output occurs
	// at intervals determined by dt_fulldump.

	n_dump = -dt_snap;
	dt_snap = VERY_LARGE_NUMBER;

	if (n_dump > 0)
	    full_dump = 1;
	else
	    full_dump = 2;
	    
	// first_full_dump added by SPZ Jan 2000 to assure initial dump 
	// of snapshot to start worldbundle.

	first_full_dump = true;
	PRC(n_dump); PRL(full_dump);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Complete the setup of static class data:

    real initial_step_limit = Starlab::min(INITIAL_STEP_LIMIT, dt_max);

    // Replaced dt_snap by dt_reinit in the next 2 lines.  (Steve, 1/05)

    initial_step_limit = Starlab::min(initial_step_limit, dt_reinit);
    real step_limit = Starlab::min(dt_max, dt_reinit);	// no initial_step_limit

    b->set_mbar(b->get_mass()/b->n_leaves());		// mass scale
    b->set_initial_step_limit(initial_step_limit);
    b->set_step_limit(step_limit);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // The scheduling is such that the system will be synchronized
    // at *every* integer multiple of step_limit.

    // The behavior of the system is governed by the following times:

    real dt_sync = step_limit;

    // Check and flag some basic relationships between time scales.

    if (dt_reinit < dt_sync)
	err_exit("dt_reinit < dt_sync.");	// require dt_reinit >= dt_sync
    if (dt_esc < dt_sync)
	err_exit("dt_esc < dt_sync.");		// require dt_esc >= dt_sync
    if (dt_sstar < dt_sync)
	err_exit("dt_sstar < dt_sync.");	// require dt_sstar >= dt_sync
    if (dt_fulldump < dt_sync)			// require dt_fulldump
	err_exit("dt_fulldump < dt_sync.");	//		>= dt_sync

    if (verbose) {
	if (dt_log < dt_sync)
	    cerr << "warning: dt_log < dt_sync, quoted errors unpredictable."
		 << endl;
	if (dt_snap < dt_sync)
	    cerr << "warning: dt_snap < dt_sync, no restart possible."
		 << endl;
	else if (dt_snap < dt_reinit)
	    cerr << "warning: dt_snap < dt_reinit, restart unpredictable."
		 << endl;
    }

    xreal t = b->get_system_time();	// current time	(was get_time()...)

    real tt = t;			// use to avoid xreal problems with
					// possible VERY_LARGE_NUMBERs below
					// -- could build this into xreal +,
					//    but may be too inefficient for
					//    general use

    real t_end = tt + delta_t;		// final time, at end of integration

    real t_log = tt;			// time of next log output
    //    if (t > (xreal)0) 
    //    	t_log += dt_log;	// uncomment to prevent initial output

    real t_snap = tt + dt_snap;		// time of next snapshot output
    if (snap_init) t_snap = tt;
    real t_sync = tt + dt_sync;		// time of next system synchronization

    // Changes by Steve (7/01):

    real t_esc = tt;    // + dt_esc;	// time of next escaper check
					// 	note: apply IMMEDIATE check
    real t_reinit = tt; // + dt_reinit;	// time of next reinitialization

    // Frequency of internal escaper checks (write to story only):

    real dt_esc_check = Starlab::min(dt_sync, 1.0/16);
    real t_esc_check = tt + dt_esc_check;

    // Time of next complete system dump (Steve, 9/01).

    real t_fulldump = tt;

    // Testing the MPI indexing...

    real t_MPI_check = tt + 0.25 * randinter(0, dt_log);

    //----------------------------------------------------------------------
    // Frequencies of "other" (episodic) output.  Idea is that we will
    // have data blocks of width dt_alt2, sampled at intervals dt_alt1,
    // and separated by intervals dt_alt3.
    //
    // Ultimately these should become command-line options (Steve, 6/02).

    real dt_alt1, dt_alt2, dt_alt3;
    real t_alt1, t_alt2;

    dt_alt1 = 0;			// default: suppress

    if (alt_flag) {
	dt_alt1 = 1./32;		// hardwired options...
	dt_alt2 = 1;
	dt_alt3 = 10;
    }

    t_alt1 = tt;			// assume that these times work
    t_alt2 = tt;			// with the alt output intervals

    //----------------------------------------------------------------------

    if (verbose) {
	cerr << endl;
	PRC(t), PRL(t_end);
	PRL(dt_sync);
	PRL(dt_log);
	PRL(long_binary_out);
	PRL(dt_snap);
	PRL(dt_esc);
	PRL(dt_reinit);
	PRL(dt_sstar);
	PRL(dt_fulldump);

	PRC(t_snap); PRC(t_log); PRC(t_sync); PRC(t_esc); PRL(t_fulldump);

	if (alt_flag) {
	    cerr << endl << "additional output:  ";
	    PRC(t_alt1); PRC(dt_alt1); PRC(dt_alt2); PRL(dt_alt3);
	}

	ko->print();
	kd->print();
    }

#ifdef DUMP_DATA
    ofstream tmp_dump("TMP_DUMP");
#endif

    // Initialize the system.  This is somewhat redundant, as immediate
    // reinitialization is scheduled within the while loop.  However,
    // we need to know time steps for fast_get_nodes_to_move().  This
    // will also compute all perturber lists.

    full_reinitialize(b, t, verbose, true);

    bool tree_changed = true;	// used by fast_get_nodes_to_move; set
    				// by the integration/evolution routines

    // cerr << "check_sync #0 at t = " << b->get_system_time() << " (";
    // xprint(b->get_system_time(), cerr, false);
    // cerr << ")" << endl;
    check_sync(b);

    while (t <= t_end) {

	int n_next;
	xreal ttmp;

	// Create the new time step list.

	fast_get_nodes_to_move(b, next_nodes, n_next, ttmp, tree_changed);

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 1 << endl << flush;

	    if (T_DEBUG_LEVEL > 0) {
		int p = cerr.precision(HIGH_PRECISION);
		cerr << "evolve_system: "; PRC(n_next); PRL(ttmp);
		if (n_next < 4) {
		    PRI(15);
		    for (int ii = 0; ii < n_next; ii++)
			cerr << next_nodes[ii]->format_label() << " ";
		    cerr << endl << flush;
		}
		cerr.precision(p);
	    }
	}
#endif

#ifndef TMP_SUPPRESS
	bool quit_now = false;
	hdyn *bad = NULL;
	for (int ii = 0; ii < n_next; ii++) {
	    if (!next_nodes[ii]->is_valid()) {
		cerr << "next_node[" << ii << "] = " << next_nodes[ii]
		     << " is invalid" << endl << flush;
		bad = next_nodes[ii];
		quit_now = true;
	    }
	}
	if (quit_now) {
	    for_all_nodes(hdyn, b, bb)
		if (bb == bad) cerr << "invalid node contained within the tree"
				    << endl;
	    err_exit("invalid node(s) on timestep list.");
	}
#endif

	// New order of actions (Steve, 7/01):
	//
	//	1. log output
	//	2. flag escapers
	//	3. check STOP
	//	4. snapshot/full dump output
	//	5. remove escapers
	//	6. reinitialize

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	bool save_snap = false;

	if (ttmp > t_log) {

	    // 1. Standard log output.

	    // System should be synchronized here, but it is not essential.
	    // However, energy errors may not be reliable, and the system
	    // may not be restartable (if save_snap_at_log is set).

	    real cpu0 = cpu_time();

	    // Encode binary log data in form suitable for sys_stats.

	    int long_bin = 0;
	    if (long_binary_out > 0 && dt_log > 0) {
		int nlog = (int)rint((real)ttmp/dt_log);
		if (nlog%long_binary_out == 0)
		    long_bin = 2;
		else
		    long_bin = 1;
	    }

	    log_output(b, count, steps, count_top_level, steps_top_level,
		       &kc_prev, long_bin);

//#ifdef USEMPI
	    if (!check_MPI_id_array(b)) {
		recompute_MPI_id_array(b);
		if (!check_MPI_id_array(b)) exit(1);
	    }
//#endif

	    // (Incorporate count, steps, etc. into kira_counters...?)

	    save_snap = save_snap_at_log;

	    cerr << endl << "Total CPU time for log output = "
		 << cpu_time() - cpu0 << endl;

	    cerr << endl; flush(cerr);
	    update_step(ttmp, t_log, dt_log);
	}

	//=================================================================

	if (dt_alt1 > 0 && ttmp > t_alt1) {	// episodic output (time
						// intervals hardcoded)
	    kira_alt_output(b);

	    // Update t_alt1 and check for gaps.

	    while (t_alt1 < ttmp) t_alt1 += dt_alt1;
	    if (t_alt1 > t_alt2 + dt_alt2) {
		t_alt1 += dt_alt3 - dt_alt2 - dt_alt1;
		t_alt2 += dt_alt3;
	    }

	    // PRC(b->get_system_time()); PRC(t_alt1); PRL(t_alt2);
	}

	//=================================================================

	if (full_dump && ttmp > t_esc_check) {

	    // 2. Flag escapers (full_dump mode only).  System may or may not
	    // be synchronized, but dyn functions don't know about pred_pos...

	    refine_cluster_mass(b, 0);
	    update_step(ttmp, t_esc_check, dt_esc_check);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if (ttmp > t_sync) {

	    // 3. System is synchronized.
	    // Check to see if we should stop the run.

	    if (check_file("STOP")) {

		t_end = ttmp - (xreal)1;    // Forces optional snap output 
					    // and end of run in snap_output.

		cerr << endl << "***** Calculation STOPped by user at time "
		     << t << " *****" << endl << endl;
	    }

	    // cerr << endl << "check_sync #1 at t = "
	    // 	 << b->get_system_time() << " (";
	    // xprint(b->get_system_time(), cerr, false);
	    // cerr << ")" << endl;

	    tree_changed |= check_sync(b);
	    update_step(ttmp, t_sync, dt_sync);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// 4. Output a snapshot to cout at the scheduled time, or at end of
	// run. NOTE: This is the *only* place where output to cout can occur.

	// Not essential to have dt_snap >= dt_reinit, but should have this
	// if we want full restart capability.

	// Added special case delta_t = 0 to force snap output.

	bool reg_snap = (ttmp > t_snap || delta_t == 0
			 || (dt_snap < VERY_LARGE_NUMBER && ttmp > t_end)
			 || (fmod(steps, 1000) == 0 && check_file("DUMP")));

	if (reg_snap || save_snap) {

	    // PRC(ttmp), PRC(t_snap), PRL(dt_snap);

	    snap_output(b, steps, snaps,
			reg_snap, save_snap, snap_save_file,
			t, ttmp, t_end, t_snap, dt_snap, verbose);

	    if (reg_snap)
		update_step(ttmp, t_snap, dt_snap);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// 5. Full dumps should precede and follow major changes when the
	// system is properly synchronized.

	// bool full_dump_now = (full_dump
	//		      && (t == 0 || t >= t_reinit
	//			  	 || t>= t_esc || t >= t_sync));

	// No full output at every t_sync.
	// Stellar evolution forces t_sync to be 1/64 or less!

	// bool full_dump_now = (full_dump
	//		      && (t >= t_reinit || t >= t_esc || t >= t_end));
	//
	// Note that t = t_reinit = t_esc now at restart, so the second
	// clause of the if () test is true initially.

	// Now have a separate timescale for full dump.

	bool full_dump_now = (full_dump
			      && (t >= t_fulldump || t >= t_end));
	if (full_dump)
	    update_step(ttmp, t_fulldump, dt_fulldump);

	// This dump ends the current worldbundle, so don't do it initially.
	// First full dump will be written below.

	if (full_dump_now && !first_full_dump) {

	    // We want the last full dump (t_end) to be a complete
	    // snapshot, so we can restart from the data at the end
	    // of the worldbundle file.  (Steve, 7/01)

	    int short_output = 1;
	    if (ttmp > t_end) short_output = 4;		// new (7/01)

	    set_complete_system_dump(true);
	    put_node(b, cout,
		     false,		// don't print xreal
		     short_output);	// short output (uses STARLAB_PRECISION)
	    set_complete_system_dump(false);

	    cerr << "Full dump (";
	    if (short_output == 1) cerr << "t";
	    else cerr << "h";
	    cerr << "dyn format) at time " << t << endl;
	}

	first_full_dump = false;

#if 0
	real tcheck1 = 17.0;
	real tcheck2 = 17.5;
	real tcheck3 = 18.0;
	real tcheck4 = 19.0;
	if (   t <= tcheck1 && ttmp > tcheck1
	    || t <= tcheck2 && ttmp > tcheck2
	    || t <= tcheck3 && ttmp > tcheck3
	    || t <= tcheck4 && ttmp > tcheck4 ) {

	    char name[10];
	    sprintf(name, "TMP_%.1f", (real)t);
	    ofstream f(name);
	    if (f) {

		// put_node(b, f);
		pp3(b, f);
		cerr << "wrote TMP file " << name << endl;

		f.close();
	    }
	}
#endif

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if (ttmp > t_end) return;		// return ==> stop

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	bool reinit = (ttmp > t_reinit);
	if (reinit) update_step(ttmp, t_reinit, dt_reinit);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	bool esc = false;

	if (ttmp > t_esc) {

	    // 5. Check for and remove escapers.

	    // System is reinitialized and time step list is
	    // recomputed if any stars are removed.

	    // *** REQUIRE dt_esc >= dt_sync. ***

	    int n_prev = b->n_leaves();

	    if (b->get_scaled_stripping_radius() > 0)
		check_and_remove_escapers(b, ttmp, next_nodes,
					  n_next, tree_changed);

	    int n_new = b->n_leaves();

	    if (n_new <= n_stop) {
		cerr << "N_stop reached\n";
		exit(0);
	    }

	    if (n_new != n_prev) {
		reinit = true;
		esc = true;
	    }

	    update_step(ttmp, t_esc, dt_esc);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if (reinit) {

	    // 6. Reinitialize the system.

	    // *** REQUIRE dt_reinit >= dt_sync. ***

	    // Offset and reset here mimic what happens to data on output
	    // and restart, where the com pos and vel are subtracted off
	    // by put_node(), then added back by get_node(), which may
	    // cause the original and restarted data to differ in the last
	    // decimal place.

	    b->reset_com();
	    b->offset_com();
	    full_reinitialize(b, t, verbose, false);
	    tree_changed = true;

	    fast_get_nodes_to_move(b, next_nodes, n_next, ttmp, tree_changed);

	    if (esc) cerr << endl << flush;
	}

	// NOTE:  If dt_log and dt_reinit are properly chosen, escapers will
	// always be checked and the system will always be reinitialized
	// immediately after a snapshot or full dump, so continuation should
	// be the same as a restart from that snapshot.

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Second full_dump output marks the start of a new worldbundle.

	if (full_dump_now) {

	    set_complete_system_dump(true);
	    put_node(b, cout, false, 1);
	    set_complete_system_dump(false);

	    cerr << endl << "Full dump (tdyn format) at time " << t << endl;
	}

	//-----------------------------------------------------------------

	// Proceed to the next step.

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 2 << endl << flush;
	}
#endif

	if (ttmp < t)
	    backward_step_exit(b, ttmp, t, next_nodes, n_next);

#ifdef CPU_COUNTERS
	cpu_prev = cpu;
	kc->cpu_time_other += (cpu = cpu_time()) - cpu_prev;
#endif

	// Force the elder of two binary sisters to be integrated.

	check_binary_scheduling(b, next_nodes, n_next);

	xreal t_prev = b->get_system_time();

#if 0
	cerr << "old time = "; xprint(t, cerr, true);
	cerr << "new time = "; xprint(ttmp, cerr, true);
	cerr << "xreal dt = "; xprint((xreal)next_nodes[0]->get_timestep(),
				      cerr, true);
#endif

	b->set_system_time(t = ttmp);
	b->set_time(t);

	// Take a new step to time t (now system_time):

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 3 << endl << flush;
	}
#endif

	if (!b->has_grape()) {
	
	    // We need full prediction if any top-level nodes appear on
	    // the integration list.  Otherwise, we really only need to 
	    // predict the perturbers (done in perturber_acc_and_jerk...).
	    // However, if a perturber list is invalid or overflows, then
	    // we need full prediction again...

	    bool predict_all = false;
	    for (int i = 0; i < n_next; i++) {
		hdyn* n = next_nodes[i];
		if (n->is_top_level_node()
		    || !n->get_top_level_node()->get_valid_perturbers()) {
		    predict_all = true;
		    break;
		}
	    }

	    if (predict_all)
		predict_loworder_all(b, t);	

	    // GRAPE does (top-level) prediction, if present.
	}

	//-----------------------------------------------------------------

	// Integrate all particles in the current block step.
	// Steps counts particle steps; count counts block steps.

	if (0 && full_dump) {
	    cerr << endl << "integration list (n_next = "
		 << n_next << "):" << endl;
	    for (int i = 0; i < n_next; i++) {
		PRI(4); PRC(i); cerr << next_nodes[i]->format_label() << endl;
	    }
	}

	if (full_dump == 1) {

	    // The next_nodes[] array is used to test for changes in
	    // unperturbed status.  The last_dt[] array saves timesteps
	    // at the start of the step.
	    //
	    // Not a particularly clean way to proceed...
	    //
	    //					(Steve, 10/00, 8/01)

	    for (int i = 0; i < n_next; i++) {

		if (next_nodes[i]->get_kepler())
		    next_flag[i] = true;
		else
		    next_flag[i] = false;

		last_dt[i] = next_nodes[i]->get_timestep();
	    }

	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {

	    cerr << "DEBUG: evolve_system " << 4 << endl << flush;

	    cerr << "entering integrate_list: "; PRC(t); PRC(n_next);
	    cerr << next_nodes[0]->format_label();
	    if (n_next > 1) cerr << " " << next_nodes[1]->format_label();
	    cerr << endl;

	    if (T_DEBUG_LEVEL > 0) {
		cerr << endl << "before integrate_list: " << endl;
		pp3("(15,100015)");
		pp3("(21,100021)");
		pp3("(23,100023)");
	    }
	}
#endif

//----------------------------------------------------------------------

#include "kira_pre_step_debug.C"	// for temporary debugging


	// FINALLY, we take the block step....

	int n_list_top_level = 0;
	int ds = integrate_list(b, next_nodes, n_next, exact,
				tree_changed, n_list_top_level,
				full_dump, r_reflect);

//#ifdef USEMPI
	if (tree_changed) recompute_MPI_id_array(b);	// *way* overkill:
							// should just remove
							// the old nodes and
							// add the new...

	// Testing...

//	if (t >= t_MPI_check) {
//	    if (!check_MPI_id_array(b)) exit(1);
//	    t_MPI_check += 0.25 * randinter(0, dt_log);
//	}

//#endif


	if (n_list_top_level > 0) {

	    // Count blocks and steps for top-level nodes only.
	    // (For GRAPE, count actual hardware calls.)

	    count_top_level += 1;
	    steps_top_level += n_list_top_level;
	}


#include "kira_post_step_debug.C"	// for temporary debugging

//----------------------------------------------------------------------

#if 1

	// Check for and repair timestep disparities in near neighbors.
	// Code should eventually move to integrate_list...

	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *bb = next_nodes[ii];
	    if (bb && bb->is_top_level_node()) {

		// Check for a close NN.

		hdyn *nn = bb->get_nn();

		if (nn && nn != bb
		    && bb->get_d_nn_sq() < 0.01*bb->get_d_min_sq()) {

		    // Check and possibly reduce its time step.

		    if (nn->is_top_level_node()
			&& nn->get_timestep() > 16*bb->get_timestep()) {

			cerr << "NN check: synchronized NN "
			     << nn->format_label() << " to time "
			     << bb->get_system_time() << endl;
			cerr << "after step for " << bb->format_label()
			     << endl;

			// Synchronize_node should do the right thing,
			// but may be better to use integrate_list, which
			// will use GRAPE if available.  However, will
			// have to add flags to avoid tree updates, etc.
			// in that case.

			nn->synchronize_node();
		    }
		}
	    }
	}

#endif

#ifdef CPU_COUNTERS
	cpu_prev = cpu;
	kc->cpu_time_integrate += (cpu = cpu_time()) - cpu_prev;
#endif

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp) && T_DEBUG_LEVEL > 0) {
	    cerr << "DEBUG: evolve_system " << 45 << endl << flush;
	    cerr << "after integrate_list: " << endl;
	    pp3("(15,100015)");
	    pp3("(21,100021)");
	    pp3("(23,100023)");
	}
#endif

	bool force_energy_check = false;
	if (ds < 0) {
	    ds = -ds;
	    if (kd->check_heartbeat) force_energy_check = true;
	}

	steps += ds;
	grape_steps += ds;
	count += 1;

	if (force_energy_check)
	    count = 4*kd->n_check_heartbeat
			* (floor(count / (4*kd->n_check_heartbeat)) + 1);

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 5 << endl << flush;
	}
#endif

	if (full_dump == 1) {

	    // Check for dump of worldline information.

	    for (int i = 0; i < n_next; i++) {

		hdyn *curr = next_nodes[i];

		if (curr && curr->is_valid()) {

		    // Always dump new tree information (see hdyn_tree.C),
		    // but allow the possibility of fewer dumps during
		    // normal integration.  Ordinarily, dump every n_dump
		    // steps.

		    bool dump_now = (fmod(curr->get_steps(), n_dump) == 0);

		    // Force dump if node's unperturbed status has changed.

		    dump_now |=
			(next_flag[i] && !curr->get_kepler()
			 || (!next_flag[i] && curr->get_kepler()));

#if 0
		    // *** New treatment of lightly perturbed binaries.
		    // *** What about slow binaries?  Later...

		    if (curr->is_low_level_node()
			&& !curr->get_kepler()
			&& curr->get_perturbation_squared()
				< SMALL_TDYN_PERT_SQ) {

			// Which criterion to use?

			int which = 2;

			if (which == 1) {

			    // Attempt #1: Only dump if we have just crossed
			    //	           the parent step:

			    real t_par = curr->get_parent()->get_time();
			    dump_now = (t >= t_par && t - last_dt[i] < t_par);

			    // However, the parent step for a significantly
			    // perturbed binary may not be significantly greater
			    // than the internal step, so the gain may not be
			    // very great.  The parent timestep is sensitive
			    // to "wiggles" due to perturbers, and increases
			    // significantly as the perturbation increases.

			} else if (which == 2) {

			    // Attempt #2: Simply increase the dump timescale
			    //	           to come closer to a few bound orbits.
			    //	           Possibly could be corrected to take
			    //	           into account whether or not the
			    //	           orbit is bound; however, if we simply
			    //	           interpolate from the start to the end
			    //	           of an unbound segment, should be OK.

			    dump_now = (fmod(curr->get_steps(),
					     10*n_dump) == 0);
			}

			if (dump_now && !curr->get_elder_sister()) {
			    cerr << endl
				 << curr->get_parent()->format_label()
				 << " is lightly perturbed; get_steps = "
				 << curr->get_steps() << endl;
			    PRL(sqrt(Starlab::max(0.0,
					 curr->get_perturbation_squared())));
			    real ratio = curr->get_parent()->get_timestep()
					   / curr->get_timestep();
			    PRL(curr->get_parent()->get_timestep());
			    PRC(curr->get_timestep());
			    PRL(ratio);
			}
		    }
#endif

		    if (dump_now) {

			// cerr << "kira: time " << b->get_system_time();
			// cerr << "  put_node " << i << "/" << n_next << "  "
			//      << curr->format_label() << endl;

			put_single_node(curr, cout, false, 1);

			// Note that we have to use binary_sister here, not
			// younger_sister, because this node may become the
			// younger sister in a new binary.

			if (curr->is_low_level_node()) {

			    // cerr << "      put_node for sister "
			    //	    << curr->get_binary_sister()
			    //	                    ->format_label()
			    //	    << endl;

			    put_single_node(curr->get_binary_sister(),
					    cout, false, 1);
			}
		    }
		}
	    }
	}

#ifdef DUMP_DATA
	if (tmp_dump && t >= 17.0 && t <= 17.5) {

	    // Dump everything to file tmp_dump.

	    tmp_dump << endl << "n_next = " << n_next << endl;
	    for (int i = 0; i < n_next; i++)
		if (next_nodes[i] && next_nodes[i]->is_valid()) {
		    tmp_dump << "i = " << i << " (" << next_nodes[i] << "):"
			     << endl;
		    pp3(next_nodes[i], tmp_dump, -1);
		} else
		    tmp_dump << "i = " << i << " (" << next_nodes[i]
			     << "):  invalid" << endl;
	}
#endif

	// check_slow_consistency(b);

	// Optionally check the heartbeat...

	// Not much point continuing if t - t_prev = 0...
	// Eventually may want more diagnostics.  For now, just quit!

	real dt_true = t - t_prev;

	if (dt_true <= 0
	    || (kd->check_heartbeat
		&& fmod(count, kd->n_check_heartbeat) == 0)) {

	    PRC(count), PRC(t), PRC(t-t_prev), PRL(cpu_time());

#if 0
	    PRL(n_next);
  	    for (int i = 0; i < n_next && i < 5; i++)
  		if (next_nodes[i] && next_nodes[i]->is_valid()) {
  		    PRI(4);
		    cerr << i << " " << next_nodes[i]->format_label()
			 << " " << next_nodes[i]->get_timestep() << endl;
//  		    pp3_minimal(next_nodes[i]);
  		}
#endif

	    if (fmod(count, 4*kd->n_check_heartbeat) == 0)
		print_recalculated_energies(b);
	}

	if (dt_true <= 0) {
	    cerr << endl << "time step error at ";
	    cerr.precision(HIGH_PRECISION);
	    PRL(b->get_system_time());
	    PRC(t); PRL(t_prev);

	    cerr << endl << "integration list (n_next = "
		 << n_next << "):" << endl << endl;

	    for (int ii = 0; ii < n_next; ii++)
		if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
		    PRC(ii); PRL(next_nodes[ii]->format_label());
		    PRI(9); PRL(next_nodes[ii]->get_timestep());
		    hdyn *top = next_nodes[ii]->get_top_level_node();
		    PRL(top->format_label());
		    cerr << "binary properties: " << endl;
		    for_all_nodes(hdyn, top, bi) {
			hdyn *od = bi->get_oldest_daughter();
			if (od) {
			    print_binary_from_dyn_pair(od,
					   od->get_younger_sister());
			    cerr << endl;
			}
		    }

		    cerr << endl << "pp3 output:" << endl;
		    pp3(top);
		}

	    err_exit("zero effective time step");
	}

	// Integrate_list handles tree changes -- allows one per step.

	// Note: if a node is destroyed as part of the step,
	// the corresponding entry in next_nodes returns NULL.
	// (Old convention -- could now simply check node->is_valid().)

	// Finally, evolve the stellar system (stars and binaries).
	// Entire system will be reinitialized if necessary.

	// Note: "if" statement moved from evolve stars() (see note there).
	// Both stars and binaries are currently updated at fixed time
	// intervals (previously updated binaries after every binary step
	// -- pros and cons exist for either choice; details are still
	// under investigation).

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 6 << endl << flush;
	}
#endif

	if (fmod(b->get_system_time(), dt_sstar) == 0.0
	    && b->get_use_sstar())  {

	   // cerr << "pre SE at t = " << b->get_system_time() << endl;
	   // print_recalculated_energies(b);

	    bool tmp = evolve_stars(b, full_dump);

#ifdef T_DEBUG
	    if (IN_DEBUG_RANGE(ttmp)) {
		if (tmp)
		    cerr << "tree change caused by evolve_stars at time "
			 << b->get_system_time() << endl << flush;
	    }
#endif

	    tree_changed |= tmp;

	    // cerr << "post SE at t = " << b->get_system_time() << endl;
	    // print_recalculated_energies(b);

	    // print_energy_from_pot(b);	// should be free, but
						// doesn't quite work...

#if 0

	    real tstellev = VERY_LARGE_NUMBER;
	    for_all_daughters(hdyn, b, bb) {
		real ts = bb->get_starbase()->get_evolve_timestep();
		if (ts < tstellev) tstellev = ts;
	    }
	    PRL(tstellev);


#endif

	}

	if (tree_changed) set_n_top_level(b);

	if (b->has_grape()) {

	    // Check for GRAPE release.

	    if (grape_steps > ko->grape_check_count) {
		grape_steps = 0;
		check_release_grape(b->get_config(), ko, t);
	    }
	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(ttmp)) {
	    cerr << "DEBUG: evolve_system " << 7 << endl << flush;
	}
#endif

	// Miscellaneous checks (see kira_runtime.C):

	if (fmod(count, kd->n_check_runtime) == 0) {

	    if (cpu_time() > cpu_time_limit) {
		cerr << endl
		     << "***** CPU time limit of " << cpu_time_limit
		     << " seconds exceeded"
		     << endl;
		return;				// pretty harsh...
	    }

	    real new_dt_log = 0;
	    real new_dt_snap = 0;
	    char* new_snap_save_file = new char[128];
	    new_snap_save_file[0] = '\0';

	    bool status
		= check_kira_runtime(b, t_end,
				     new_dt_log, new_dt_snap,
				     long_binary_out, new_snap_save_file,
				     tree_changed);

	    if (new_dt_log > 0) {
		dt_log = new_dt_log;
		int i = (int)((real)b->get_system_time() / dt_log);
		t_log = (i+1) * dt_log;
	    }

	    if (new_dt_snap > 0) {
		dt_snap = new_dt_snap;
		int i = (int)((real)b->get_system_time() / dt_snap);
		t_snap = (i+1) * dt_snap;
	    }

	    if (new_snap_save_file[0] != '\0') {
		save_snap_at_log = true;
		snap_save_file = new_snap_save_file;
	    } else
		delete [] new_snap_save_file;

	    if (status) return;

	    ko = b->get_kira_options();		// (just in case...)
	    kd = b->get_kira_diag();
	}
    }
}



// Not currently used...

local void recompute_perturber_lists(hdyn *b, bool verbose = false)
{
    int n_nolist = 0;
    cerr << "checking perturber lists #1" << endl << flush;
    for_all_daughters(hdyn, b, bb)
	if (bb->is_parent() && !bb->get_oldest_daughter()->get_kepler()
	    && !bb->get_perturber_list()) {
	    if (verbose)
		cerr << bb->format_label()
		     << " is perturbed, no perturber list"
		     << endl;
	    n_nolist++;
	}
    PRL(n_nolist);

    cerr << "computing acc and jerk for binary CM nodes" << endl << flush;
    calculate_acc_and_jerk_on_top_level_binaries(b);

    n_nolist = 0;
    cerr << "checking perturber lists #2" << endl << flush;
    for_all_daughters(hdyn, b, bb)
	if (bb->is_parent() && !bb->get_oldest_daughter()->get_kepler()
	    && !bb->get_perturber_list()) {
	    if (verbose)
		cerr << bb->format_label()
		     << " is perturbed, no perturber list"
		     << endl;
	    n_nolist++;
	}
    PRL(n_nolist);
}



// Standalone function to perform kira operations: integrate a system
// for a specified period of time.

void kira(hdyn * b,	       // hdyn array
	  real delta_t,	       // time span of the integration
	  real dt_log,	       // time step of the integration
	  int long_binary_out,
	  real dt_snap,	       // snapshot output interval
	  real dt_sstar,       // timestep for single stellar evolution
	  real dt_esc,	       // escaper removal
	  real dt_reinit,      // reinitialization interval
	  real dt_fulldump,    // full dump interval
	  bool exact,	       // exact force calculation
	  real cpu_time_limit,
	  bool verbose,
	  bool snap_init,
	  bool save_snap_at_log, // save snap at log output
	  char* snap_save_file,  // filename to save in
	  int n_stop,	         // when to stop
	  bool alt_flag)	 // enable alternative output

{

  // Control the behavior of the kepler package.

  set_kepler_tolerance(2);
  set_kepler_print_trig_warning(b->get_kira_diag()
				->report_kepler_trig_error);

  // It is desirable to compute perturber lists to avoid excessive
  // computation in any initial perturbed binaries.  This will be done
  // in evolve_system() when the system is initialized.

  evolve_system(b, delta_t, dt_log, long_binary_out,
		dt_snap, dt_sstar, dt_esc, dt_reinit, dt_fulldump,
		exact, cpu_time_limit,
		verbose, snap_init, save_snap_at_log, snap_save_file,
		n_stop, alt_flag);

}



// Clean up prior to exiting kira.

void kira_finalize(hdyn *b)
{
    cerr << endl << "Finalize kira:"<<endl;
    cerr << "End of run at time " << b->get_system_time()
	 << endl
	 << "Total CPU time for this segment = " << cpu_time()
	 << endl;

    //--------------------------------------------------------------------
    // For unknown reasons, neomuscat sometimes dumps core at end of run...
    //--------------------------------------------------------------------

    // exit(0);			// may still cause a core dump...

    // kill(getpid(), 9);    	// dumb, but brute force works!

    //--------------------------------------------------------------------

    // Clean up static data (to make it easier for ccmalloc to find
    // real memory leaks).

    cerr << endl << "cleaning up hdyn_schedule" << endl << flush;
    clean_up_hdyn_schedule();

    cerr << "cleaning up hdyn_ev" << endl << flush;
    clean_up_hdyn_ev();

    cerr << "cleaning up kira_ev" << endl << flush;
    clean_up_kira_ev();

    if (b->has_grape()) {
	cerr << "cleaning up hdyn_grape" << endl << flush;
	clean_up_hdyn_grape(b->get_config());
    }

    cerr << "cleaning up next_nodes" << endl << flush;
    if (next_nodes) delete [] next_nodes;

    cerr << "cleaning up kira_counters" << endl << flush;
    if (b->get_kira_counters()) delete b->get_kira_counters();

    cerr << "cleaning up perturbed_list" << endl << flush;
    if (b->get_perturbed_list()) delete [] b->get_perturbed_list();

    cerr << "cleanup complete" << endl << flush;

#ifdef USEMPI
    cerr << "Now finalizing mpi:" << endl << flush;
    MPI_Finalize();
#endif

}

#else



// The standalone program.

#include <unistd.h>				// for termination below...
#include <signal.h>

#include <time.h>

// This function has to live in the same file as kira:main(), since it must
// reflect the environment at the time kira was built.  It is referenced by
// kira_initialize.

void kira_system_id(int argc, char** argv)
{
    // Identify the Starlab version.

    cerr << "Starlab version " << VERSION;
    cerr << "; kira created on " << __DATE__ << " at " << __TIME__ << endl;

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
    int jarg = 1;
    for (int i = 1; i < argc; i++) {
	cerr << " " << argv[i];
	if (((jarg++ > 14 && argv[i][0] != '-') || jarg > 18) && i < argc-1) {
	    cerr << " \\" << endl << "           ";
	    jarg = 1;
	}
    }
    cerr << endl;

#ifndef USEMPI

    // wwvv For now, don't try to find out the location of kira
    // when using mpi.

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

#endif

    cerr << "current directory = " << getenv("PWD") << endl;
    cerr << endl;
}

#ifdef USEMPI
void signalhandler(int s)
{
    abort();
}
#endif

main(int argc, char **argv) {

#ifdef USEMPI
    signal(11, signalhandler);
#endif

    check_help();

    // Parameters to be passed directly to kira:

    hdyn *b;			// hdyn root node

    // Time scales -- all should be powers of 2.

    real delta_t;		// time span of the integration
    real dt_log;		// output interval
    real dt_snap;		// snap output interval
    real dt_sstar;		// standard single-star evolution timestep
    real dt_esc;		// timestep to remove escapers
    real dt_reinit;		// timestep to reinitialize entire system
    real dt_fulldump;		// timestep for full system dumps

    int long_binary_out;

    real cpu_time_limit;

    bool exact;			// force calculation using perturber list
    bool verbose;		// Toggle "verbose" mode

    bool snap_init = false;
    bool save_snap_at_log = false;
    char snap_save_file[256];
    int  n_stop;		// n to terminate simulation

    bool alt_flag;

#ifdef USEMPI

    struct rlimit rlim;
    if (getrlimit(RLIMIT_CORE,&rlim)<0)
      perror("getrlimit: ");
    rlim.rlim_cur = rlim.rlim_max;
    if (setrlimit(RLIMIT_CORE,&rlim)<0) 
      perror("setrlimit ");
    cerr << "core limits" << rlim.rlim_max <<" "<< rlim.rlim_cur << endl;

    MPI_Init(&argc,&argv);
    mpi_communicator = MPI_COMM_WORLD;
    MPI_Comm_size(mpi_communicator,&mpi_nprocs);
    cerr << "Number of mpi processes: "<<mpi_nprocs<<endl;
    MPI_Comm_rank(mpi_communicator,&mpi_myrank);

#endif

    if (!kira_initialize(argc, argv,
			 b, delta_t, dt_log, long_binary_out,
			 dt_snap, dt_sstar,
			 dt_esc, dt_reinit, dt_fulldump,
			 exact, cpu_time_limit, verbose,
			 snap_init, save_snap_at_log, snap_save_file,
			 n_stop, alt_flag))
	get_help();

    kira(b, delta_t, dt_log, long_binary_out, 
	 dt_snap, dt_sstar, dt_esc, dt_reinit, dt_fulldump,
	 exact, cpu_time_limit,
	 verbose, snap_init, save_snap_at_log, snap_save_file,
	 n_stop, alt_flag);

    kira_finalize(b);

    // rmtree(b);
    // rmtree(b, false);	// experiment: don't delete root node

    //--------------------------------------------------------------------
    // For unknown reasons, merlot & halley dump core at end of run if
    // rmtree is used.
    //--------------------------------------------------------------------

}

#endif
