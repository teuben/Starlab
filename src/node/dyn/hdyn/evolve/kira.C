
//#define DUMP_DATA 1	// uncomment to allow detailed TMP_DUMP output

       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// kira:  Hermite N-body integrator with evolving hierarchical tree
////        structure, stellar and binary evolution, and an external
////        tidal field.
////
//// Basic options are as follows.  Default values for a new simulation
//// are indicated in square brackets.  For restart (continuation of a
//// previous kira calculation), the defaults for many options [*] are
//// determined from the input snapshot, making it possible to continue
//// a run without having to re-specify all the command-line parameters
//// previously used.  (These new defaults may still be overridden on
//// the command line, of course.)  Some options typically having to do
//// with initial setup may be overridden by data from the input
//// snapshot (if present), as noted below.  Kira may also turn *on*
//// some options (the B, G, Q, S, and u settings) if they were turned
//// on in the previous run, but are not specified on the current
//// command line.  To prevent this, use the "-o" switch.
////
//// The first page of log output gives detailed information on all
//// parameter settings adopted and any modifications made during
//// initialization.  In case of doubt, read the log file!
////
//// Options:     -a    specify accuracy parameter [0.1][*]
////              -b    specify frequency of full binary output, in (integer)
////                    units of the log output interval
////                                  [10; no argument or 0 ==> no output]
////              -B    turn on binary evolution [off][*]
////              -c    include comment [none]
////              -C    specify GRAPE release interval, in seconds [15]
////              -d    specify log output interval [0.25][*]
////              -D    specify snapshot interval [end of run]
////              -e    specify softening length [0][*]
////              -E    use exact calculation [false]
////              -f    specify close-encounter distance [1 --> 1/N][*]
////              -F    specify tidal field type (0 = none, 1 = point-mass,
////                                              2 = isothermal, 3 = disk,
////                                              4 = custom)
////                                  [0 if -Q not set; 1 otherwise][*]
////              -g    specify hysteresis factor [2.5][*]
////              -G    specify initial stripping radius [none][*]
////              -h    specify stellar-evolution time step [0.015625 = 1/64][*]
////              -I    specify (re)initialization timescale [1][*]
////              -J    specify Jacobi radius [none][*]
////              -k    specify perturbation factor [1.e-7][*]
////              -K    specify log2(maximum slowdown factor) (integer): [0][*]
////              -L    specify CPU time limit, in seconds [none]
////              -M    specify initial total mass, in solar masses
////                                  [none; take from input snap if present]
////              -n    stop at specified number of particles [10]
////              -N    specify frequency of CPU check output [50000]
////              -o    prevent kira from overriding some settings (BGQSu)
////                        based on input snapshot data [allow]
////              -O    save (and overwrite) extra snapshot at each output [no]
////              -q    specify initial virial ratio [0.5]
////              -Q    use tidal field [none][*]
////              -r    specify initial virial radius in code units
////                                  [1; take from input snap if present]
////              -R    specify initial virial radius, in pc
////                                  [none; take from input snap if present]
////              -s    specify random seed [take from system clock]
////              -S    turn on stellar evolution [off][*]
////              -t    specify time span of calculation [10]
////              -T    specify initial virial time scale, in Myr
////                                  [none; take from input snap if present]
////              -u    toggle unperturbed multiple motion [disabled][*]
////              -U    toggle all unperturbed motion [enabled][*]
////              -v    toggle "verbose" mode [on]
////              -x    toggle output of extended precision time [on]
////              -X    specify escaper removal timescale [reinit][*]
////              -y    specify stellar encounter criterion
////                                  [0 N-body units or solar radii][*]
////              -z    specify stellar merger criterion [0 stellar radii][*]
////              -Z    specify stellar tidal dissipation criterion
////                                  [0 stellar radii][*]
////
//// As a convenient shorthand, any "dt" interval specified less than zero
//// is interpreted as a power of 2, i.e. "-d -3" sets dt_log = 0.125.

// Level-2 help:

//++ Some run-time parameters:
//++
//++ The initial virial radius is read from the input snapshot.  If no
//++ initial virial radius is found there, the value specified by the
//++ "-V" command-line option is used.  If no "-V' option is set, a
//++ default value of 1 is assumed.
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
//
//      J. Makino, S. McMillan          12/92
//      S. McMillan, P. Hut              9/94
//	J. Makino, S. McMillan          12/95
//	J. Makino, S. McMillan           7/96
//	J. Makino, S. McMillan           1/97
//	S. McMillan			 6/97
//	S. McMillan, S. Portegies Zwart	 3/98
//	S. McMillan, S. Portegies Zwart	 7/98
//	S. McMillan			 8/98
//	S. Portegies Zwart	        12/00
//
// NOTE:  ALL direct references to USE_GRAPE are confined to
//	  this file (referenced externally via the functions below).
//
// Known problems:  First timestep is a bit shaky, in particular 
//                  if the system is cold.
//
// July 1997 (Steve)
//
//   Need the ability to ensure rigorous reproducability on restart.
//   Should restructure evolve_system to allow log output, snap output,
//   escaper removal, stellar evolution, debugging to be performed
//   independently.



#ifdef TOOLBOX

#include "hdyn.h"
#include "star/dstar_to_kira.h"

#define INITIAL_STEP_LIMIT  0.0625	// Arbitrary limit on the first step

//#define CHECK_PERI_APO_CLUSTER          // Follow individual stellar orbits.

static kira_counters kc_prev; 
static bool dbg = false;



// *** All GRAPE-related calls are now defined in kira_grape_include.C ***

#ifdef USE_TREE
  #include "/work6/starlab/Tree++/kira_tree_include.C"
#else
  #include "kira_grape_include.C"
#endif

//========================================================================


// Handy local functions...

// match_label: return true iff b's id (name, index, or #index)
//              matches the specified string.

local bool match_label(hdyn* b, char* label)
{
    if (b->name_is(label))
        return true;
#if 0
    // Function format_label() prepends a # to numbers.
    // Check for the number without the #.

    if (b->get_index() >= 0) {
	char id[64];
	sprintf(id, "%d", b->get_index());
	if (streq(id, label))
	    return true;
    }
#endif
    return false;
}

// match_label_tree: return true iff the label of any
//                   node below b (or b) matches the
//                   specified string.

local bool match_label_tree(hdyn* b, char* label)
{
    for_all_nodes(hdyn, b, bi)
	if (match_label(bi, label)) return true;
    return false;
}

// cond_print_nn: for each b in the list, print out the IDs of the entire
//                subtree containing b if any of those IDs match the
//		  specified string.

local void cond_print_nn(hdyn** list, int n, char* label, char* header = NULL)
{
    for (int i = 0; i < n; i++) {
	hdyn* b = list[i];
	if (b) {
	    if (match_label_tree(b->get_top_level_node(), label)) {
		if (header) cerr << header << endl;
		cerr << "print_nn triggered by node "
		     << b->format_label() << endl;
		for_all_leaves(hdyn, b, bi)
		    print_nn(bi, 2);
	    }
	}
    }
}

local void test_kepler(hdyn *b)
{
    for_all_nodes(hdyn,b,bi) {
	if (bi->get_kepler()) {
	    if (bi->is_top_level_node()) {
		cerr << " test_kepler: top level ";
		bi->pretty_print_node(cerr); cerr << "has kepler \n";
		pp3(bi, cerr);
	    }
	}
    }
}



// Triple debugging...

local void print_triple_stats(hdyn* root, hdyn* b)
{
    if (!b->is_top_level_node()) return;
    if (b->n_leaves() != 3) return;

    hdyn* ss = b->get_oldest_daughter();
    hdyn* bs = ss->get_younger_sister();
    if (ss->is_parent()) {
	hdyn* tmp = ss;
	ss = bs;
	bs = tmp;
    }

    // Triple is (ss, bs), with ss single, bs double.

    real mss = ss->get_mass();
    vector pos_ss = b->get_pos() + ss->get_pos();
    vector vel_ss = b->get_vel() + ss->get_vel();
    hdyn* bs1 = bs->get_oldest_daughter();
    hdyn* bs2 = bs1->get_younger_sister();
    real mb1 = bs1->get_mass();
    real mb2 = bs2->get_mass();
    vector pos_bs1 = b->get_pos() + bs->get_pos() + bs1->get_pos();
    vector vel_bs1 = b->get_vel() + bs->get_vel() + bs1->get_vel();
    vector pos_bs2 = b->get_pos() + bs->get_pos() + bs2->get_pos();
    vector vel_bs2 = b->get_vel() + bs->get_vel() + bs2->get_vel();

    real pot_1 = -mss * (mb1+mb2) / abs(ss->get_pos() - bs->get_pos());
    real pot_2 = -mss * mb1 / abs(pos_ss - pos_bs1)
		 -mss * mb2 / abs(pos_ss - pos_bs2);
    real pot_b = -mb1 * mb2 / abs(pos_bs1 - pos_bs2);

    real pot_ext = 0;
    for_all_daughters(hdyn, root, x) {
	if (x != b)
	    pot_ext -= x->get_mass() * (mss/abs(x->get_pos() - pos_ss)
					+ mb1/abs(x->get_pos() - pos_bs1)
					+ mb2/abs(x->get_pos() - pos_bs2));
    }

    real ke = 0.5 * (mss*square(vel_ss) + mb1*square(vel_bs1) 
		     			+ mb2*square(vel_bs2));
    // PRC(ke), PRL(pot_ext);
    // PRC(pot_2 + pot_b + ke), PRL(pot_ext + pot_2 + pot_b + ke);

    // print_binary_from_dyn_pair(ss, bs, 0, 0, true);
    // print_binary_from_dyn_pair(bs1, bs2, 0, 0, true);

    cerr << endl << "Triple " << b->format_label() << " at time "
	 << b->get_system_time() << ":" << endl
	 << "    Eint = " << pot_2 + pot_b + ke
	 << "  Eint + phi_ext = " << pot_2 + pot_b + ke + pot_ext << endl;
}



local void print_binary_diagnostics(hdyn* bi)
{
    bool diag = true;

    kepler *kep;
    hdyn *s = bi->get_binary_sister();
    real M = bi->get_parent()->get_mass();

    if (bi->get_kepler() == NULL) {
	kep = new kepler;
	kep->set_time(bi->get_time());
	kep->set_total_mass(M);
	kep->set_rel_pos(bi->get_pos() - s->get_pos());
	kep->set_rel_vel(bi->get_vel() - s->get_vel());
	kep->initialize_from_pos_and_vel();
    } else
	kep = bi->get_kepler();

    // if (kep->get_separation() < kep->get_semi_major_axis())
    //    diag = false;

    if (diag) {

#if 0		    
	cerr << "\nLow-level node bi = ", bi->print_label(cerr);
	cerr << endl;
	PRI(4); PRL(bi->get_time());

	int p = cerr.precision(INT_PRECISION);
	PRI(4); cerr << "top-level: ",
	PRL(bi->get_top_level_node()->get_time());
	cerr.precision(p);

	pp2(bi->get_top_level_node(), cerr, 2);
	PRI(4);
	PRL(bi->get_top_level_node()->get_valid_perturbers());
	PRI(4); PRL(bi->get_top_level_node()->get_n_perturbers());
#endif

	if (bi->get_kepler())
	    cerr << "    bi unperturbed";
	else
	    cerr << "    bi perturbed";

	real r = kep->get_separation();
	real E = kep->get_energy();
	real a = r;
	if (E < 0) a = min(r, kep->get_semi_major_axis());

	cerr << ", -E/mu = " << -E
	     << "  P = " << 2*M_PI * sqrt(a*a*a/M)
	     << "  e = " << kep->get_eccentricity()
	     << endl
	     << "    sma = " << kep->get_semi_major_axis()
	     << "  r = " << r
	     << endl;

	PRI(4);
	if (bi->get_kepler()) PRC(bi->get_unperturbed_timestep());
	PRL(bi->get_timestep());
    }

    if (bi->get_kepler() == NULL) delete kep;
}



local void check_unperturbed(hdyn* bi, bool& tree_changed)
{
    // Check for unperturbed motion.

    // This is the ONLY place in kira where unperturbed motion
    // (of any sort) is initiated.

    // The unperturbed criterion is defined in function
    // is_unperturbed_and_approaching()) (see hdyn_unpert.C).

    if (bi->get_kepler() == NULL && bi->get_eps2() == 0
	&& (bi->is_unperturbed_and_approaching())) {

	// Must bring triple components up to date before
	// continuing (about to freeze the entire triple system).

	// *** In general, should bring ANY substructure
	// *** components up to date before proceeding
	// *** with unperturbed startup.

	// cerr << endl << "Unperturbed motion for "
	//      << bi->format_label()
	//      << " at time " << bi->get_time() << endl;

	if (bi->is_parent() || bi->get_binary_sister()->is_parent()) {

	    bool synch = false;

	    if (bi->is_parent()) {
		bool sync = false;
		for_all_nodes(hdyn, bi, bb) {
		    if (bb->get_time() < bi->get_time())
			sync = true;
		}
		if (sync)
		    synchronize_tree(bi);
		synch |= sync;
	    }

	    if (bi->get_binary_sister()->is_parent()) {
		bool sync = false;
		for_all_nodes(hdyn, bi->get_binary_sister(), bb) {
		    if (bb->get_time()
			< bi->get_binary_sister()->get_time())
			sync = true;

		}
		if (sync)
		    synchronize_tree(bi->get_binary_sister());
		synch |= sync;
	    }

	    // If a newly-synchronized node is not on the
	    // integration list (which must be the case, as we
	    // only synchronize if a node is not up to date), then
	    // there is a good chance that the scheduling list will
	    // be corrupted, so force the list to be recomputed.

	    if (synch)
		tree_changed = true;
	}

	// Actual initialization of unperturbed motion:

	bi->startup_unperturbed_motion();

    }
}



// Slow binary functions.

local inline void check_set_slow(hdyn *bi)
{
    // Check for initiation of slow binary motion in a (normal-speed)
    // perturbed binary.

    // Calling function has already checked that kepler, slow,
    // and elder_sister are all NULL.

    // Criteria:	(0) bound!
    //			(1) perturbation less than cutoff
    //			(2) just passed apastron
    //			(3) components are single or unperturbed
    //
    // Apply these (inline) tests before passing control to the real
    // startup function.

    // if (streq(bi->get_parent()->format_label(), "(652a,652b)")) return;

    if (bi->get_max_slow_factor() > 1
	&& (bi->is_leaf() || bi->get_oldest_daughter()->get_kepler())
	&& (bi->get_younger_sister()->is_leaf()
	    || bi->get_younger_sister()->get_kepler())
	&& bi->get_perturbation_squared()
		< bi->get_max_slow_perturbation_sq()/2
	&& bi->passed_apo()
	&& get_total_energy(bi, bi->get_younger_sister()) < 0) {

	// (Energy check shouldn't be necessary, as passed_apo should
	// never return true in an unbound system...)

#if 0
	// Don't start slow motion if there are any binaries on the
	// perturber list.

	if (has_binary_perturbers(bi)) {
	    cerr << "suppressing slow motion for "
		 << bi->get_parent()->format_label()
		 << " because of binary perturbers" << endl;
	    return;
	}
#endif

	bi->startup_slow_motion();
    }
}

local inline void check_extend_slow(hdyn *bi)
{
    // Check for extension or termination of slow binary motion.

    // Calling function has already checked that kepler and elder_sister
    // are all NULL, and slow is currently set.

    // Apply this (inline) test in an inline function before passing
    // control to the real startup function.

    if (bi->passed_apo()) {

	// Use of function passed_apo should be OK most of the time, but
	// it may fail if an orbit is nearly circular.  Best to make sure
	// that we are at the right phase of the orbit before modifying
	// the slow motion.

	// For weakly perturbed binaries, energy should be nearly conserved,
	// so the period should be a reasonably good indicator of the time
	// since the last apocenter.

	real P = get_period(bi, bi->get_younger_sister());

	if (bi->get_time() - bi->get_slow()->get_t_apo() > 0.9 * P)

	    bi->extend_or_end_slow_motion(P);
    }
}



local void merge_and_correct(hdyn* b, hdyn* bi, hdyn* bcoll)
{
    // This intermediate function added mainly to allow
    // deletion of bi and bcoll after leaving merge_nodes.

    cerr << endl << "----------" << endl
	 << "merge_and_correct: merging node "
	 << bi->format_label();
    cerr << " with node " << bcoll->format_label() 
	 << " at time " << bi->get_system_time() << endl;

    PRC(bi), PRC(bcoll), PRC(bi->get_parent()); PRL(cpu_time());

    hdyn* cm = bi->merge_nodes(bcoll);

    delete bi;
    delete bcoll;

    b->get_kira_counters()->leaf_merge++;

    cerr << "after merge_nodes ";
    PRC(cm); PRL(cpu_time());

    cerr << "----------" << endl;
}

local void check_periapo(hdyn * bi) {

  // just book keeping for stellar orbits.

#ifdef CHECK_PERI_APO_CLUSTER
  bi->check_periapo_node();
#endif
}

local hdyn* check_and_merge(hdyn* bi)
{
    hdyn * bcoll;
    if ((bcoll = bi->check_merge_node()) != NULL) {

	// cerr << "check_and_merge: "; PRL(bcoll);

	hdyn* b = bi->get_root();

	merge_and_correct(b, bi, bcoll);

	// Check for multiple mergers.  Note that we check *all*
	// stars for merging, whether or not they are in the
	// current block.
	//
	// This is perhaps more than we really want (Steve, 12/98).

	bool merge_flag = true;
	while (merge_flag) {

	    merge_flag = false;

	    for (hdyn* bb = b;
		 (bb != NULL) && !merge_flag;   
		 bb = (hdyn*) bb->next_node(b)) {

		if (bb->is_leaf()) {
		    hdyn* bcoll2 = bb->check_merge_node();
		    if (bcoll2 != NULL) {
		        // cerr << "check_and_merge (2): "; PRL(bcoll2);
			merge_and_correct(b, bb, bcoll2);
			merge_flag = true;
		    }
		}
	    }
	}
    }
    // cerr << "return from check_and_merge: "; PRL(bcoll);

    return bcoll;
}



// Define TIME_LIST in order to time various parts of integrate_list.

//#define TIME_LIST

local int integrate_list(hdyn * b,
			 hdyn ** next_nodes, int n_next,
			 bool exact, bool & tree_changed,
			 bool full_dump = false)
{
    static bool restart_grape = true;
    static bool reset_force_correction = true;	// no longer used

    int i, steps = 0;
    xreal sys_t = next_nodes[0]->get_system_time();

    //    cerr << "At stars of Integrate_list: sys_t = " << sys_t << endl;

    // Code to time specific force-calculation operations:

    real cpu0, cpu1;

#ifdef TIME_LIST

    int kmax = 1;

    for (int k = 0; k < n_next; k++) {
	hdyn* bi = next_nodes[k];

	if (bi->name_is("(13a,13b)")
	    && bi->get_kepler() == NULL
	    && bi->find_perturber_node()
	    && bi->find_perturber_node()->get_n_perturbers() > 0
	    ) {

	    // Found the particle.  Clean up the rest of the list.

	    for (int kk = 0; kk < n_next; kk++)
		if (kk != k)
		    next_nodes[kk]->set_timestep(VERY_LARGE_NUMBER);
	    n_next = 1;
	    next_nodes[0] = bi;

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl << "timing " << bi->format_label()
		 << " at time " << bi->get_system_time() << endl;
	    if (bi->is_low_level_node()) {
		PRL(bi->get_top_level_node()->format_label());
		if (bi->find_perturber_node()) {
		    PRL(bi->find_perturber_node()->format_label());
		    PRL(bi->find_perturber_node()->get_n_perturbers());
		}
	    }
	    cerr.precision(p);

	    kmax = 25000;
	    cpu0 = cpu_time();
	}
    }
#endif

    // Separate the force calculation from the rest for GRAPE implementation.

    xreal t_next = next_nodes[0]->get_next_time();

#ifdef TIME_LIST
    for (int k = 0; k < kmax; k++) {
#endif

    calculate_acc_and_jerk_for_list(b, next_nodes, n_next, t_next,
				    exact, tree_changed,
				    reset_force_correction,  // no longer used
				    restart_grape);

#ifdef TIME_LIST
    }
    if (kmax > 1) cpu1 = cpu_time();
#endif

    // Apply corrector and redetermine timesteps.

    bool reinitialize = false;

#ifdef TIME_LIST
    for (int k = 0; k < kmax; k++) {
#endif

    bool diag = false;
    for (i = 0; i < n_next; i++) {

	hdyn *bi = next_nodes[i];

	if (bi && bi->is_valid()) {

	    if (!bi->get_kepler()) {

		if (diag) cerr << " perturbed correction for "
		               << bi->format_label() << endl;

		if (!bi->correct_and_update()) {

		    // A problem has occurred during the step.

		    // Recompute acc and jerk on the front end (no bookkeeping
		    // yet) and retry once.  (Steve 9/98)

		    cerr << "retrying force calculation for "
			 << bi->format_label() << endl;

		    // Better do an exact calculation, as we can't (yet) call
		    // correct_acc_and_jerk() to correct a single particle...
		    // This may run into problems with slow binaries, however.

		    bi->clear_interaction();
		    bi->calculate_acc_and_jerk(true);
		    bi->set_valid_perturbers(false);

		    if (bi->is_top_level_node() && b->get_tidal_field() > 0)
			add_external(bi);

		    if (!bi->correct_and_update()) {
			cerr << endl << "failed to correct error for "
			     << bi->format_label() << " at time "
			     << bi->get_system_time() << endl;
			err_exit("Run terminated in integrate_list");
		    }
		}

		bi->init_pred();
		bi->store_old_force();

		// Note that old_acc = acc at the end of a step.

	    } else {

		if (bi->get_eps2() != 0)	// Excessively cautious?
		    err_exit("integrate_list invoked with non-zero softening");

		if (diag) {
		    cerr << "unperturbed motion for "
			 << bi->format_label() << endl;
		    if (bi->get_nn())
			cerr << "nn = " << bi->get_nn()->format_label()
			     << endl;
		}

		// As of 3/99, integrate_unperturbed_motion returns true
		// iff bi is still an unperturbed binary after the step.

		hdyn* parent = bi->get_parent();
		bool top_level = parent->is_top_level_node();

		bool pert = bi->is_perturbed_cpt();   // false iff bi is
						      // fully unperturbed

		if (!bi->integrate_unperturbed_motion(reinitialize)) {

		    // Unperturbed motion is over.  Either the binary
		    // containing bi has become perturbed or it has merged.

		    if (bi->is_valid()) {

			// Parent of bi is a newly perturbed binary.
			// Add it to the list if the binary was top-level
			// and fully perturbed previously (i.e. don't
			// bother with partial unperturbed orbits).

			if (top_level && !pert)
			    parent->add_to_perturbed_list(1);

		    } else {

			// Seem to have had a merger.  Remove binary from
			// the list, if necessary.

			if (top_level && pert)
			    parent->remove_from_perturbed_list(1);

		    }
		}

		// Note that it is possible for integrate_unperturbed_motion
		// to delete bi.  Must take this possibility into account
		// below.

		if (!bi->is_valid())
		    next_nodes[i] = NULL;

	    }
	}
    }

#ifdef TIME_LIST
    }


    if (kmax > 1) {
	cerr << "CPU times: "
	     << (cpu1-cpu0)/kmax << "  "
	     << (cpu_time()-cpu1)/kmax << endl;
	exit(0);
    }
#endif

    // Complete all steps before modifying binary structure...

    diag = false;
    for (i = 0; i < n_next; i++) {

	hdyn *bi = next_nodes[i];

	if (bi && bi->is_valid()) {

	    if (!bi->is_low_level_node()) {

	       // Not much to do for a top-level node -- just update counters.

		if (bi->is_leaf())
		    b->get_kira_counters()->step_top_single++;
		else
		    b->get_kira_counters()->step_top_cm++;

		// Diagnostics:

		if (diag) {
		    cerr << "\nTop-level node bi = " << bi->format_label()
			 << " at time " << bi->get_time() << endl;
		    cerr << "timestep = " << bi->get_timestep() << endl;
		    PRL(bi->get_acc());
		    PRL(bi->get_jerk());
		}

	    } else {

		update_binary_sister(bi);

		b->get_kira_counters()->step_low_level++;

		if (bi->get_slow())
		    b->get_kira_counters()->inc_slow(bi->get_kappa());

		// print_binary_diagnostics(bi);

		// Check for new unperturbed motion.

		if (!bi->get_kepler()) {

		    if (bi->get_eps2() == 0) {

			// Check for unperturbed motion:

			check_unperturbed(bi, tree_changed);

			// Check to see if the binary containing bi has just
			// become unperturbed.

			if (bi->get_parent()->is_top_level_node()
			    && !bi->is_perturbed_cpt())
			    bi->get_parent()->remove_from_perturbed_list(2);
		    }
		}

		// Check for new or modified slow motion.  Note that we
		// recheck the kepler pointer, just in case...

		// There is some redundancy here, since the checks may
		// repeat those in is_unperturbed_and_approaching, and
		// also only want to check immediately after apocenter,
		// but worry about these issues later.
		//					(Steve, 7/99)

		if (!bi->get_kepler()) {

		    // Only check on elder sister (should never see younger
		    // sister in any case).

		    if (!bi->get_elder_sister()) {

			if (bi->get_slow())
			    check_extend_slow(bi);
			else
			    check_set_slow(bi);

			// Need to update counters for slow motion here.
		    }
		}
	    }

	    steps++;
	}
    }

    // Probably makes more sense to check for encounters for all stars
    // before testing for mergers, tree changes, etc. (Steve, 3/24/00).

    if (b->get_stellar_encounter_criterion_sq() > 0)
        for (i = 0; i < n_next; i++) 
	    check_print_close_encounter(next_nodes[i]);

    // ONLY ONE tree reconstruction (following a collision or
    // otherwise) is currently permitted per block step.

    //++ Note from Steve to Steve, 7/98.  Could we relax the
    //++ requirement of only one tree reconstruction per step?

    for (i = 0; i < n_next; i++) {

	// Note somewhat convoluted calling sequence to merge nodes:
	//
	//	check_and_merge
	//	    - calls check_merge_node  to locate bcoll (if any)
	//	    - calls merge_and_correct to merge bi and bcoll
	//		  + calls merge_nodes to do the merging and clean up
	//		  + *deletes* both nodes bi and bcoll
	//	    - attempts to handle multiple mergers.
	//
	// Routine merge_and_correct takes care of all corrections
	// associated with mergers.
	//
	// ** This is now mostly handled by merge_nodes (Steve, 3/9/00). **

	hdyn *bi = next_nodes[i];

	if (bi && bi->is_valid()) {

	    // first check for peri- or apoclustron passage
	    check_periapo(bi);

	    hdyn* bcoll = check_and_merge(bi);

	    // cerr << "integrate_list: "; PRL(bcoll);

	    if (bcoll) {

		// Merger occurred and tree has to be rebuilt.  No further
		// tree reorganization is permitted during this block step.
		// Perturber and other lists should already be up to date,
		// and the merging nodes have already been replaced by
		// their center of mass.

		// If bcoll is non-NULL, then both bi and bcoll have
		// *already* been deleted...

		// PRC(bi), PRL(bcoll);




		// *** Must have check_and_merge take care of full_dump
		// *** output in this case...




		tree_changed = true;
		restart_grape = true;
		reset_force_correction = true;	// no longer used

		synchronize_tree(b);
		steps += b->n_leaves();

		cerr << "call initialize_system_phase2() "
		     << "from integrate_list [1]"
		     << " at time " << b->get_system_time() << endl;

		initialize_system_phase2(b, 1);
		b->reconstruct_perturbed_list();

		PRL(tree_changed);

		// Remove merged star and its merger companion from the
		// integration list.  Note that the list may be still
		// be incomplete on return, as multiple merger companions
		// may remain on it.

		next_nodes[i] = NULL;
		for (int j = 0; j < n_next; j++) {
		    if (next_nodes[j] == bcoll)
			next_nodes[j] = NULL;
		}

		return steps;
	    }
	}
    }

    if (reinitialize) {

	synchronize_tree(b);
	steps += b->n_leaves();

	cerr << "call initialize_system_phase2() from integrate_list [2]"
	     << " at time " << b->get_system_time() << endl;

	initialize_system_phase2(b, 2);
	b->reconstruct_perturbed_list();

	tree_changed = true;
	restart_grape = true;
	reset_force_correction = true;	// no longer used

    } else {

	for (i = 0; i < n_next; i++) {

	    hdyn *bi = next_nodes[i];

	    if (bi && bi->is_valid()) {

		hdynptr* cm_list = NULL;
		int n_list = 0;

		if (bi->is_low_level_node() || bi->is_parent()) {

		    // Make a list of CM nodes in the current clump.

		    for_all_nodes(hdyn, bi, bb)
			if (bb->is_parent()) n_list++;

		    if (n_list > 0) {
			cm_list = new hdynptr[n_list];
			n_list = 0;
			for_all_nodes(hdyn, bi, bb)
			    if (bb->is_parent()) cm_list[n_list++] = bb;
		    }
		}

		// Check (and modify, if necessary) the tree structure.
		// Save some data on bi, in case the tree is restructured.

		hdyn *top_level = bi->get_top_level_node();
		hdyn *od = NULL, *yd = NULL;
		bool pert = false;

		if (top_level->is_parent()) {
		    od = top_level->get_oldest_daughter();
		    yd = od->get_younger_sister();
		    pert = od->is_perturbed_cpt();
		}

		int adjust_tree = bi->adjust_tree_structure(full_dump);

		if (adjust_tree) {

		    b->get_kira_counters()->tree_change++;

		    reset_force_correction = true;	// no longer used
		    restart_grape = true;
		    tree_changed = true;

		    // Check to see if the unperturbed binary list needs
		    // to be updated.

		    // As of 3/99, if adjust_tree_structure returns true (> 0),
		    // its value indicates the type of adjustment:
		    //
		    //	0	no adjustment
		    //	1	low-level combine
		    //	2	top-level combine (direct)
		    //	3	top-level split (direct)
		    //	4	top-level combine (indirect)
		    //	5	top-level split (indirect)
		    //  6	low-level synchronization
		    //
		    // Top-level changes may be driven directly by top-level
		    // nodes (2 and 3), or they may be forced by changes at
		    // lower levels (4 and 5).  The unpleasant logic below
		    // seems unavoidable given the construction of
		    // adjust_tree_structure() and the requirements of
		    // perturbed_list.

		    // PRC(adjust_tree); PRL(pert);

		    if (adjust_tree == 1) {

			// Low-level combine, but still possible that the
			// identity of the top-level node has changed.

			if (bi->get_top_level_node() != top_level
			    && pert) {
			    top_level->remove_from_perturbed_list(3);
			    bi->get_top_level_node()->add_to_perturbed_list(2);
			}

		    } else if (adjust_tree == 2) {

			// Top-level combine.  Node bi is now one component
			// of a new top-level binary.  (In this case, bi is
			// the same as top_level.)

			if (pert) bi->remove_from_perturbed_list(4);

			if (bi->is_perturbed_cpt())
			    bi->get_parent()->add_to_perturbed_list(3);

			yd = bi->get_binary_sister();

			if (yd->is_perturbed_cm())
			    yd->remove_from_perturbed_list(5);

		    } else if (adjust_tree == 3) {

			// Top-level split.  Binary CM node bi = top_level
			// no longer exists, but its components od and yd
			// are now at the top level.

			// Update perturbed_list as necessary.

			if (pert) top_level->remove_from_perturbed_list(6);

			if (od->is_perturbed_cm())
			    od->add_to_perturbed_list(4);

			if (yd->is_perturbed_cm())
			    yd->add_to_perturbed_list(5);

		    } else if (adjust_tree == 4) {

			// Top-level combine was induced by internal
			// reorganization.  Node bi is now part of a new
			// top-level binary, but we don't necessarily
			// know which component...

			top_level = bi->get_top_level_node();
			od = top_level->get_oldest_daughter();
			yd = od->get_younger_sister();

			if (od->is_perturbed_cm())
			    od->remove_from_perturbed_list(7);

			if (yd->is_perturbed_cm())
			    yd->remove_from_perturbed_list(8);

			if (od->is_perturbed_cpt())
			    top_level->add_to_perturbed_list(6);

		    } else if (adjust_tree == 5) {

			// Top-level split was induced by internal
			// reorganization.  Node top_level no longer
			// exists, but its components od and yd are now
			// at the top level.

			if (pert) top_level->remove_from_perturbed_list(9);

			if (od->is_perturbed_cm())
			    od->add_to_perturbed_list(7);

			if (yd->is_perturbed_cm())
			    yd->add_to_perturbed_list(8);

		    }

#ifndef USE_GRAPE

		    if (b->get_kira_diag()->kira_main) {
			cerr << "\nAfter adjusting tree structure... \n";
			cerr << "Time = " << next_nodes[0]->get_system_time()
			     << " single_steps = "
			     << b->get_kira_counters()->step_top_single
			     << endl;
			print_recalculated_energies(b);
			// pp3(bi->get_top_level_node(), cerr);
			flush(cerr);
		    }

#endif

		    // Set next_nodes[i] = NULL for any center of mass in
		    // the clump that has been adjusted, so we can use the
		    // rest of the array in evolve_system on return from
		    // integrate_list.

		    if (cm_list && n_list > 0) {
			for (int j = 0; j < n_next; j++) {
			    if (!next_nodes[j]->is_valid())
				next_nodes[j] = NULL;
			    else
				for (int k = 0; k < n_list; k++)
				    if (next_nodes[j] == cm_list[k]) {
					next_nodes[j] = NULL;
					break;
				    }
			}
		    }
		    if (cm_list) delete [] cm_list;

		    return steps;	// NOTE: we currently return after
					//	 the FIRST tree rearrangement,
					// 	 so we can only have one
					//	 restructuring per block time
					//	 step.
		}
		if (cm_list) delete [] cm_list;
	    }
	}
    }

    return steps;
}



local void full_reinitialize(hdyn* b, xreal t, bool verbose)
{
    cerr << "\nReinitializing system at time " << t << endl;

    real cpu_0 = cpu_time();

    b->set_system_time(t);
    // b->to_com();
    initialize_system_phase1(b, t);
    initialize_system_phase2(b, 3,
			     false);			// "false" here means
							// timesteps are only
							// set if zero

    b->reconstruct_perturbed_list();

    // Quietly reinitialize all kepler structures.
    // (Repeats actions taken by get_hdyn() on startup.)

    b->initialize_unperturbed();

    if (verbose) {
	cerr << "CPU time for reinitialization = "
	     << cpu_time() - cpu_0 << endl;
    }
}

local bool check_sync(hdyn* b)
{
    bool need_new_list = false;

    for_all_nodes(hdyn, b, bb)
	if (bb->get_time() != b->get_system_time()
	     && bb->get_kepler() == NULL) {
	  cerr << "check_sync warning:  node "
	       << bb->format_label() << " not synchronized at time "
	       << b->get_system_time() << "; node time = " << bb->get_time()
	       << endl;
	  need_new_list = true;
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
    put_node(cerr, *b, b->get_kira_options()->print_xreal);
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
	vector vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	kin += bb->get_mass()*square(vel);
    }

    top_pot *= 0.5;
    kin *= 0.5;

    int p = cerr.precision(INT_PRECISION);
    PRC(top_pot), PRC(int_pot), PRL(kin);
    PRL(de_external_pot(b));
    cerr << "energy_from_pot = " << kin+top_pot+int_pot+0.5*de_external_pot(b)
	 << endl;
    cerr.precision(p);
}

// evolve_system:  Main integration loop.

static hdyn **next_nodes = NULL;

local void evolve_system(hdyn * b,	       // hdyn array
			 real delta_t,	       // time span of the integration
			 real dt_log,	       // time step of the integration
			 int long_binary_out,
			 real dt_snap,	       // snapshot output interval
			 real dt_sstar,	       // timestep for single
                                               //     stellar evolution
			 real dt_esc,	       // escaper removal
			 real dt_reinit,       // reinitialization interval
			 bool exact,	       // exact force calculation
			 real cpu_time_limit,
			 bool verbose,
			 bool save_last_snap,  // save snap at log output
			 char* snap_save_file, // filename to save in
			 int n_stop )	       // when to stop
{
    // Initialization:

    clean_up_files();
    cpu_init();

    initialize_counters_from_log(b);
    kc_prev = *b->get_kira_counters();		// (make a local copy)

    kira_options *ko = b->get_kira_options();
    kira_diag *kd = b->get_kira_diag();

    real count = 0;
    real steps = 0;
    int grape_steps = 0;
    int snaps = 0;

    next_nodes = new hdyn *[2 * b->n_daughters()];	// Overkill?

    // bool next_flag[2 * b->n_daughters()];		// fails on SGI
    bool *next_flag = new bool[2 * b->n_daughters()];	// equivalent?

    set_n_top_level(b);

    // Establish limits on time steps.  Note that, other than the static
    // step_limit data, the "dt" data never propagate beyond this function,
    // so there is no particular reason to make them part of the hdyn class.

    real dt_max = min(dt_log, dt_sstar);	// system will be synchronized
						// for stellar evolution occur
    dt_max = min(dt_max, dt_reinit);
    b->set_unpert_step_limit(dt_max);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Convenient for now to pass request for full dump and also the
    // dump interval via dt_snap...

    bool full_dump = false;
    real n_dump = 1.0;				// up to at least 50 seems OK!

    bool first_full_dump = false;
    if (dt_snap < 0) {
	n_dump = -dt_snap;
	dt_snap = VERY_LARGE_NUMBER;
	first_full_dump = full_dump = true;
	PRL(n_dump);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Complete the setup of static class data:

    real initial_step_limit = min(INITIAL_STEP_LIMIT, dt_max);
    initial_step_limit = min(initial_step_limit, dt_snap);

    real step_limit = min(dt_max, dt_snap); 	// no initial_step_limit, note

    b->set_mbar(b->get_mass()/b->n_leaves());	// mass scale
    b->set_initial_step_limit(initial_step_limit);
    b->set_step_limit(step_limit);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // The scheduling is such that the system will be synchronized
    // at *every* integer multiple of step_limit.

    // The behavior of the system is governed by the following times:

    real dt_sync = step_limit;

    // Check and flag some basic relationships between time scales.

    if (dt_reinit < dt_sync)
	err_exit("dt_reinit < dt_sync.");	// Require dt_reinit >= dt_sync
    if (dt_esc < dt_sync)
	err_exit("dt_esc < dt_sync.");		// Require dt_esc >= dt_sync
    if (dt_sstar < dt_sync)
	err_exit("dt_sstar < dt_sync.");	// Require dt_sstar >= dt_sync

    if (verbose) {
	if (dt_log < dt_sync)
	    cerr << "Warning: dt_log < dt_sync, quoted errors unpredictable."
		 << endl;
	if (dt_snap < dt_sync)
	    cerr << "Warning: dt_snap < dt_sync, no restart possible."
		 << endl;
	else if (dt_snap < dt_reinit)
	    cerr << "Warning: dt_snap < dt_reinit, restart unpredictable."
		 << endl;
    }

    xreal t = b->get_time();		// current time

    real tt = t;			// use to avoid xreal problems with
					// possible VERY_LARGE_NUMBERs below
					// -- could build this into xreal +,
					//    but may be too inefficient for
					//    general use

    real t_end = tt + delta_t;		// final time, at end of integration
    real t_log = tt;			// time of next log output
//    if (t > 0) t_log += dt_log;
    real t_snap = tt + dt_snap;		// time of next snapshot output
    real t_sync = tt + dt_sync;		// time of next system synchronization
    real t_esc = tt + dt_esc;		// time of next escaper check
					// Note: do NOT apply immediate check
    real t_reinit = tt + dt_reinit;	// time of next reinitialization

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

	PRC(t_snap); PRC(t_log); PRC(t_sync); PRL(t_esc);

	ko->print();
	kd->print();
    }

#ifdef DUMP_DATA
    ofstream tmp_dump("TMP_DUMP");
#endif

    // Initialize the system.

    full_reinitialize(b, t, verbose);

    bool tree_changed = true;	// Used by fast_get_nodes_to_move.
    				// Set by integration/evolution routines.

    // Note: the system can be reinitialized in several places, resulting
    // in considerable duplication of effort and inefficiency.  It remains
    // to be seen whether this produces any noticeable effect on run time.

    while (t <= t_end) {

	int n_next;
	xreal ttmp;

	// Create the new time step list.

	// get_nodes_to_move(b, next_nodes, n_next, ttmp);
	fast_get_nodes_to_move(b, next_nodes, n_next, ttmp, tree_changed);

	//------------------------------------------------------------
	// Full dumps should precede and follow major changes when the
	// system is properly synchronized.  Tests below should be
	// somewhat redundant.

//	bool full_dump_now = (full_dump
//			      && (t == 0 || t >= t_reinit
//				  	 || t>= t_esc || t >= t_sync));

	// No full output at every t_sync.
	// stellar evolution forces t_sync to be 1/64!

	// May need to be modified, as masses can change and
	// display software may not know that...

	bool full_dump_now = (full_dump
			      && (t == 0 || t >= t_reinit
				  	 || t>= t_esc));

	// This dump ends the current worldbundle, so don't do it
	// at time t = 0.

	if (full_dump_now && t > 0)
	    put_node(cout, *b,
		     false,		// don't print xreal
		     1);		// short output (uses STARLAB_PRECISION)

	// first_full_dump Added by [SPZ Jan 2000] to assure initial dump 
	// of snapshot to start worldbundle.

	if (first_full_dump) {
	  full_dump_now = true;
	  first_full_dump = false;
	}

	//------------------------------------------------------------

	if (ttmp > t_reinit) {

	    // Time to reinitialize the system.

	    // *** REQUIRE dt_reinit >= dt_sync. ***

	    full_reinitialize(b, t, verbose);
	    tree_changed = true;

	    fast_get_nodes_to_move(b, next_nodes, n_next, ttmp, tree_changed);
	    update_step(ttmp, t_reinit, dt_reinit);
	}

	bool last_snap = false;

	if (ttmp > t_log) {

	    // Standard log output.

	    // System should be synchronized here, but it is not essential.
	    // However, energy errors may not be reliable, and the system
	    // may not be restartable (if save_last_snap is set).

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

	    log_output(b, count, steps, &kc_prev, long_bin);

	    // (Incorporate count and steps into kira_counters...?)

	    last_snap = save_last_snap;

	    cerr << endl << "Total CPU time for log output = "
		 << cpu_time() - cpu0 << endl;

	    cerr << endl; flush(cerr);
	    update_step(ttmp, t_log, dt_log);

	}

	if (ttmp > t_esc) {

	    // Check for and remove escapers.

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

		// Should be somewhat redundant, but...

		full_reinitialize(b, t, verbose);
		tree_changed = true;
		fast_get_nodes_to_move(b, next_nodes, n_next, ttmp,
				       tree_changed);

		// log_output(b, steps);

		cerr << endl; flush(cerr);
	    }

	    update_step(ttmp, t_esc, dt_esc);
	}

	if (ttmp > t_sync) {

	    // System is synchronized.
	    // Check to see if we should stop the run.

	    if (check_file("STOP")) {

		t_end = ttmp - 1;	// Forces optional snap output 
					// and end of run in snap_output.

		cerr << "\n***** Calculation STOPped by user *****\n\n";
	    }

	    tree_changed |= check_sync(b);
	    update_step(ttmp, t_sync, dt_sync);
	}

	// Output a snapshot to cout at the scheduled time, or at end of run.
	// NOTE: This is the *only* place where output to cout can occur.

	// Not essential to have dt_snap >= dt_reinit, but should have this
	// if we want full restart capability.

	bool reg_snap = (ttmp > t_snap
			 || (dt_snap < VERY_LARGE_NUMBER && ttmp > t_end)
			 || (fmod(steps, 1000) == 0 && check_file("DUMP")));

	if (reg_snap || last_snap) {

	    // log_output(b, steps);


	    // PRC(ttmp), PRC(t_snap), PRL(dt_snap);

	    snap_output(b, steps, snaps,
			reg_snap, last_snap, snap_save_file,
			t, ttmp, t_end, t_snap, dt_snap, verbose);

	    if (reg_snap)
		update_step(ttmp, t_snap, dt_snap);
	}


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

		pp3(b, f);
//		put_node(f, *b);
		cerr << "wrote TMP file " << name << endl;

		f.close();
	    }
	}
#endif


	if (ttmp > t_end) return;

	//------------------------------------------------------
	// Complete the full_dump output, marking the start of
	// a new worldbundle.

	if (full_dump_now) put_node(cout, *b, false, 1);

	//------------------------------------------------------

	// Proceed to the next step.

	if (ttmp < t)
	    backward_step_exit(b, ttmp, t, next_nodes, n_next);
	    
	// Force the elder of two binary sisters to be integrated.

	check_binary_scheduling(b, next_nodes, n_next);

	xreal t_prev = b->get_system_time();

	b->set_system_time(t = ttmp);
	b->set_time(t);

	// Take a new step to time t (now system_time):

#ifndef USE_GRAPE
	
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

#endif

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

	if (full_dump) {

	    // The only current use of this array is to test for changes
	    // in unperturbed status.  This may expand as we refine the
	    // full_dump output.  Not a clean way of doing this...
	    // (Steve, 10/00)

	    for (int i = 0; i < n_next; i++)
		if (next_nodes[i]->get_kepler())
		    next_flag[i] = true;
		else
		    next_flag[i] = false;

	}

	int ds = integrate_list(b, next_nodes, n_next, exact,
				tree_changed, full_dump);

	steps += ds;
	grape_steps += ds;
	count += 1;

	if (full_dump) {
	    for (int i = 0; i < n_next; i++)
		if (next_nodes[i] && next_nodes[i]->is_valid()) {

		    // Always dump new tree information (see hdyn_tree.C),
		    // but allow the possibility of fewer dumps during
		    // normal integration.

		    // Force dump if node's unperturbed status has changed.

		    bool force_dump =
			(next_flag[i] && !next_nodes[i]->get_kepler()
			 || (!next_flag[i] && next_nodes[i]->get_kepler()));

		    if (force_dump
			|| fmod(next_nodes[i]->get_steps(), n_dump) == 0) {

//			cerr << "kira: time " << b->get_system_time();
//			cerr << "  put_node " << i << "/" << n_next << "  "
//			     << next_nodes[i]->format_label() << endl;

			put_single_node(cout, *next_nodes[i], false, 1);

			// Note that we have to use binary_sister here, not
			// younger_sister, because this node may become the
			// younger sister in a new binary.

			if (next_nodes[i]->is_low_level_node()) {

//			    cerr << "      put_node for sister "
//				 << next_nodes[i]->get_binary_sister()
//				     		 ->format_label()
//				 << endl;

			    put_single_node(cout,
					    *(next_nodes[i]
					        ->get_binary_sister()),
					    false, 1);
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


//  	for (int i = 0; i < n_next; i++)
//  	  if (next_nodes[i] && next_nodes[i]->is_valid()
//  	      && next_nodes[i]->name_is("(1752,101752)")) {
//  	    plot_stars(next_nodes[i]->get_top_level_node());
//  	    pp3(next_nodes[i], cerr, -1);
//  	  }


	//-------------------------------------------------------------------
	// (Removed numerous debugging examples to kira_debug.C -- Steve 8/98)
	//-------------------------------------------------------------------


	// check_slow_consistency(b);

	// Optionally check the heartbeat...

	if (kd->check_heartbeat && fmod(count, kd->n_check_heartbeat) == 0) {
	    PRC(count), PRC(t), PRC(t-t_prev), PRL(cpu_time());
	    if (fmod(count, 4*kd->n_check_heartbeat) == 0)
		print_recalculated_energies(b);
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

	if (fmod(b->get_system_time(), dt_sstar) == 0.0
	    && b->get_use_sstar())  {

	   //cerr << "pre SE at t = " << b->get_system_time() << endl;
	   //print_recalculated_energies(b);

	    tree_changed |= evolve_stars(b, full_dump);

	   //cerr << "post SE at t = " << b->get_system_time() << endl;
	   // print_recalculated_energies(b);

	    // print_energy_from_pot(b);	// should be free, but
						// doesn't quite work...
	}

	if (tree_changed) set_n_top_level(b);

#ifdef USE_GRAPE

	// Check for GRAPE release.

	if (grape_steps > ko->grape_check_count) {
	    grape_steps = 0;
	    check_release_grape(ko, t);
	}

#endif

	// Miscellaneous checks (see kira_runtime.C):

	if (fmod(count, kd->n_check_runtime) == 0) {

	    if (cpu_time() > cpu_time_limit) {
		cerr << endl
		     << "***** CPU time limit of " << cpu_time_limit
		     << " seconds exceeded"
		     << endl;
		return;
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
		save_last_snap = true;
		snap_save_file = new_snap_save_file;
	    } else
		delete [] new_snap_save_file;

	    if (status) return;

	    ko = b->get_kira_options();		// (just in case...)
	    kd = b->get_kira_diag();
	}
    }
}



main(int argc, char **argv)
{
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

    int long_binary_out;

    real cpu_time_limit;

    bool exact;			// force calculation using perturber list
    bool verbose;		// Toggle "verbose" mode

    bool save_last_snap = false;
    char snap_save_file[256];
    int  n_stop;		// n to terminate simulation

    if (!kira_initialize(argc, argv,
			 b, delta_t, dt_log, long_binary_out,
			 dt_snap, dt_sstar,
			 dt_esc, dt_reinit,
			 exact, cpu_time_limit, verbose,
			 save_last_snap, snap_save_file, n_stop))
	get_help();

    // b->to_com();	// don't modify input data -- use tool to do this

    // Control behavior of the kepler package.

    set_kepler_tolerance(2);
    set_kepler_print_trig_warning(b->get_kira_diag()
				   ->report_kepler_trig_error);

    evolve_system(b, delta_t, dt_log, long_binary_out,
		  dt_snap, dt_sstar, dt_esc, dt_reinit,
		  exact, cpu_time_limit,
		  verbose, save_last_snap, snap_save_file, n_stop);


    //--------------------------------------------------------------------

    // Clean up static data (to make it easier for ccmalloc to find
    // real memory leaks).

    clean_up_hdyn_schedule();
    clean_up_hdyn_ev();
    clean_up_kira_ev();

#if defined(USE_GRAPE)
    clean_up_hdyn_grape();
#endif

    if (next_nodes) delete [] next_nodes;

    if (b->get_kira_counters()) delete b->get_kira_counters();
    if (b->get_perturbed_list()) delete [] b->get_perturbed_list();

    // rmtree(b);
    // rmtree(b, false);	// experiment: don't delete root node

    //--------------------------------------------------------------------
    // for unknown reasons, merlot & halley dump core at end of run if
    // rmtree is used.  *** TO BE FIXED... ***
    //--------------------------------------------------------------------

}

#endif
