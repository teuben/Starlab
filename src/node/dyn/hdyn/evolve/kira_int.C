
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                       
//=======================================================//              /|\ ~


// kira_int:  This file contains only integrate_list() and helper
//	      functions.  Integrate_list is called only from
//	      evolve_system() in kira.C.
//
//	      Externally visible function:
//
//		  int integrate_list
//
// Major reorganization to remove compile-time GRAPE selection.
//					      Steve 6/04


#include "hdyn.h"
#include "star/dstar_to_kira.h"

#include "kira_timing.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira
#   undef T_DEBUG
#endif

//#define CHECK_PERI_APO_CLUSTER	// follow individual stellar orbits


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
    vec pos_ss = b->get_pos() + ss->get_pos();
    vec vel_ss = b->get_vel() + ss->get_vel();
    hdyn* bs1 = bs->get_oldest_daughter();
    hdyn* bs2 = bs1->get_younger_sister();
    real mb1 = bs1->get_mass();
    real mb2 = bs2->get_mass();
    vec pos_bs1 = b->get_pos() + bs->get_pos() + bs1->get_pos();
    vec vel_bs1 = b->get_vel() + bs->get_vel() + bs1->get_vel();
    vec pos_bs2 = b->get_pos() + bs->get_pos() + bs2->get_pos();
    vec vel_bs2 = b->get_vel() + bs->get_vel() + bs2->get_vel();

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
	if (E < 0) a = Starlab::min(r, kep->get_semi_major_axis());

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
		    synchronize_tree(bi);	// OK because bi is not root
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
    //			(1) perturber list is valid
    //			(2) perturbation less than cutoff
    //			(3) just passed apastron
    //			(4) components are single or unperturbed
    //
    // Apply these (inline) tests before passing control to the real
    // startup function.

    // if (streq(bi->get_parent()->format_label(), "(652a,652b)")) return;

    if (bi->get_max_slow_factor() > 1
	&& (bi->is_leaf() || bi->get_oldest_daughter()->get_kepler())
	&& (bi->get_younger_sister()->is_leaf()
	    || bi->get_younger_sister()->get_kepler())
	&& bi->get_valid_perturbers()
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



local void merge_and_correct(hdyn* b, hdyn* bi, hdyn* bcoll, int full_dump)
{
    // This intermediate function added mainly to allow
    // deletion of bi and bcoll after leaving merge_nodes.

    cerr << endl << "----------" << endl
	 << "merge_and_correct: merging node "
	 << bi->format_label();
    cerr << " with node " << bcoll->format_label() 
	 << " at time " << bi->get_system_time() << endl;

    PRC(bi), PRC(bcoll), PRL(bi->get_parent());
    PRL(bi->get_parent()->format_label());
    PRL(cpu_time());
    // pp3(bcoll);

    // Note from Steve (9/01): Modified merge_nodes() to accept the
    // full_dump flag.  For a simple merger, could place the put_node()
    // calls here, but merge_nodes() may modify the tree (if bi and bcoll
    // aren't binary sisters), and this must be properly documented.

    hdyn* cm = bi->merge_nodes(bcoll, full_dump);

    delete bi;
    delete bcoll;

    b->get_kira_counters()->leaf_merge++;

    cerr << "after merge_nodes ";
    PRC(cm); PRL(cpu_time());

    cerr << "----------" << endl;
}

local inline void check_periapo(hdyn * bi) {

  // Just bookkeeping for stellar orbits.

#ifdef CHECK_PERI_APO_CLUSTER
    bi->check_periapo_node();
#endif
}

local hdyn* check_and_merge(hdyn* bi, int full_dump)
{
    hdyn * bcoll;
    if ((bcoll = bi->check_merge_node()) != NULL) {

	cerr << "check_and_merge: "; PRL(bcoll->format_label());
	// pp3(bcoll);

	hdyn* b = bi->get_root();

	merge_and_correct(b, bi, bcoll, full_dump);

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
			merge_and_correct(b, bb, bcoll2, full_dump);
			merge_flag = true;
		    }
		}
	    }
	}
    }
    // cerr << "return from check_and_merge: "; PRL(bcoll);

    return bcoll;
}



// integrate_list:  Do the work of actually advancing in time
//		    the nodes on the specified list.
//
//		    Called only from evolve_system().

int integrate_list(hdyn * b,
		   hdyn ** next_nodes, int n_next,
		   bool exact, bool & tree_changed,
		   int& n_list_top_level,
		   int full_dump,
		   real r_reflect)
{
    // Advance the list next_nodes, containing n_next hdyn pointers.
    // Return the total number of steps taken (including retrys and
    // reinitializations), or minus that number if a GRAPE calculation
    // had to be repeated on the front end.

    static bool restart_grape = true;
    static bool reset_force_correction = true;	// no longer used

    int return_fac = 1;

    int i, steps = 0;
    xreal sys_t = next_nodes[0]->get_system_time();

    //    cerr << "At start of integrate_list: sys_t = " << sys_t << endl;

    kira_counters *kc = b->get_kira_counters();

#ifdef CPU_COUNTERS
    real cpu = cpu_time(), cpu_prev;
#endif

#ifdef TIME_LIST

    // Code to time specific force-calculation operations:

    real cpu0, cpu1;
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

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 1 << endl << flush;
    }
#endif

    // NEW FORMULATION (Steve, 4/03): calculate_acc_and_jerk_for_list()
    // will partially sort the list as it goes, placing the top-level
    // nodes at the start.  Return value now is the number of top-level
    // nodes on the list.  Note also the change in arguments.

    n_list_top_level =
	calculate_acc_and_jerk_for_list(next_nodes, n_next, t_next,
					exact, tree_changed,
					reset_force_correction,  // not used
					restart_grape);

    // Next loop through top-level nodes must turn off the
    // on_integration_list() flags...

#ifdef CPU_COUNTERS
    cpu_prev = cpu;
    kc->cpu_time_total_force += (cpu = cpu_time()) - cpu_prev;
#endif

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 2 << endl << flush;
    }
#endif

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

	    bi->clear_on_integration_list();  // cleanup: see note in kira_ev.C

#ifdef CPU_COUNTERS
	    cpu = cpu_time();
#endif

	    if (!bi->get_kepler()) {

		if (diag) cerr << " perturbed correction for "
		               << bi->format_label() << endl;

#ifdef T_DEBUG
		if (IN_DEBUG_RANGE(sys_t)) {
		    cerr << "DEBUG: integrate_list " << 3 << endl << flush;
		}
#endif

		if (!bi->correct_and_update()) {

#ifdef T_DEBUG
		    if (IN_DEBUG_RANGE(sys_t)) {
			cerr << "DEBUG: integrate_list " << 4 << endl << flush;
		    }
#endif

		    // A problem has occurred during the step, presumably
		    // because of a hardware error on the GRAPE.

		    // Recompute acc and jerk on the front end (no bookkeeping
		    // yet) and retry once.  (Steve 9/98)

		    if (bi->get_kira_diag()->grape
			&& bi->get_kira_diag()->grape_level > 0)
			cerr << "retrying force calculation for "
			     << bi->format_label() << endl;

		    // Better do an exact calculation, as we can't (yet) call
		    // correct_acc_and_jerk() to correct a single particle...
		    // This may run into problems with slow binaries, however.

		    bi->clear_interaction();
		    bi->calculate_acc_and_jerk(true);		// painful...
		    bi->set_valid_perturbers(false);

		    if (bi->is_top_level_node()
			&& b->get_external_field() > 0) {
			vec acc, jerk;
			real pot;
			get_external_acc(bi,
					 bi->get_pred_pos(),
					 bi->get_pred_vel(),
					 pot, acc, jerk);
			bi->inc_pot(pot);
			bi->inc_acc(acc);
			bi->inc_jerk(jerk);
		    }

		    if (!bi->correct_and_update()) {

			cerr << endl
			     << "Failed to correct hardware error for "
			     << bi->format_label() << " at time "
			     << bi->get_system_time() << endl << endl;
			err_exit("Run terminated in integrate_list");

		    } else {

			// Should we always print a message here, or only
			// if GRAPE diagnostics are enabled...?

			cerr << endl
			     << "Corrected apparent GRAPE"
			     << " error for "
			     << bi->format_label() << " at time "
			     << bi->get_system_time();

			if (bi->has_grape4())
			    cerr << " (chip " << get_grape4_chip(bi) << ")";

			cerr << endl;

			if (bi->get_kira_diag()->grape) {
			    if (bi->get_kira_diag()->grape_level > 0) {
				cerr << "recomputed  "; PRL(bi->get_acc());
				PRI(12); PRL(bi->get_jerk());
				// PRI(12); PRL(bi->get_pos());
				// PRI(12); PRL(bi->get_vel());
				cerr << endl;
			    }
			    return_fac = -1;
			}
		    }
		}

		bi->init_pred();
		bi->store_old_force();

		// Note that old_acc = acc at the end of a step.

#if 0
		if (bi->is_parent()
		    && (bi->name_is("(5394,21337)")
			|| bi->name_is("(21337,5394)"))) {
		    cerr << endl;
		    pp3(bi);
		    cerr << endl;
		    print_nn(bi, 1);
		    PRC(bi->get_nn()); PRL(bi->get_d_nn_sq());
		    cerr << endl;
		    bi->print_pert();
		    cerr << endl;
		}
#endif

#ifdef CPU_COUNTERS
		cpu_prev = cpu;
		kc->cpu_time_correct += (cpu = cpu_time()) - cpu_prev;
#endif

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

		int unpert = bi->integrate_unperturbed_motion(reinitialize);

		if (unpert == 0) {

		    // Unperturbed motion is over.  Either the binary
		    // containing bi has become perturbed or it has merged.

		    if (bi->is_valid()) {

			// Parent of bi is a newly perturbed binary.
			// Add it to the list if the binary was top-level
			// and fully perturbed previously (i.e. don't
			// bother with partial unperturbed orbits).

			if (top_level && !pert)
			    parent->add_to_perturbed_list(1);

			//--------------------------------------------------
			//
			// New (Steve, 9/03): the end of unperturbed motion
			// may have caused the parent or the parent nn node
			// to have been synchronized out of order.  Check
			// for and correct this, if necessary.  Correction
			// entails recomputing the scheduling list...
			//
			// Indicator: par and pnn are up to date, but one
			// or both aren't on the integration list.
			//
			// New (Steve, 2/04): this may fail if the neighbor
			// has changed while the parent was being updated!
			// Also, the logic becomes more tricky if par->nn
			// and pnn are not the same thing...  Best to rely
			// on the resched flag.

			if (!pert) {

			    story *s = b->get_root()->get_dyn_story();
			    if (find_qmatch(s, "resched")) {
				
				tree_changed = true;	// force a new
							// timestep list

				cerr << "kira: recomputing scheduling list "
				     << "at time " << sys_t << endl;

				rmq(s, "resched");
			    }
			}

			//--------------------------------------------------

		    } else {

			// Seem to have had a merger.  Remove binary from
			// the list, if necessary.

			if (top_level && pert)
			    parent->remove_from_perturbed_list(1);

		    }

		} else if (unpert == 2)

		    // Still unperturbed, but rescheduling is needed.

		    tree_changed = true;

		// Note that it is possible for integrate_unperturbed_motion
		// to delete bi, or to rearrange the clump containing bi.
		// Must take these possibilities into account below.

		if (!bi->is_valid())
		    next_nodes[i] = NULL;

#ifdef CPU_COUNTERS
		cpu_prev = cpu;
		kc->cpu_time_unperturbed += (cpu = cpu_time()) - cpu_prev;
#endif

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

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 5 << endl << flush;
    }
#endif

    if (r_reflect > 0) {
	for (i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    if (bi && bi->is_valid()) {
		real ri2 = square(bi->get_pos());

		// Apply reflection if necessary.  Don't worry about error
		// introduced in the jerk -- particle is at the edge of the
		// system, so jerk change should be small (and no effect on
		// other particles, as this option should only be used when
		// ignore_internal is true.

		if (ri2 > 1) {
		    real vr = bi->get_pos() * bi->get_vel();
		    bi->set_vel(bi->get_vel() - 2*vr*bi->get_pos()/ri2);
		}

		// cerr << endl << "reflected " << bi->format_label()
		//      << " at time " << bi->get_time() << endl << flush;
	    }
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 6 << endl << flush;
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
		    kc->step_top_single++;
		else
		    kc->step_top_cm++;

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

		kc->step_low_level++;

		if (bi->get_slow())
		    kc->inc_slow(bi->get_kappa());

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

#ifdef CPU_COUNTERS
    cpu_prev = cpu;
    kc->cpu_time_final_step += (cpu = cpu_time()) - cpu_prev;
#endif

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 7 << endl << flush;
    }
#endif

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

	    // First check for peri- or apoclustron passage.

	    check_periapo(bi);

	    hdyn* bcoll = check_and_merge(bi, full_dump);

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

		kira_synchronize_tree(b);
		steps += b->n_leaves();

#if 0
		cerr << "call initialize_system_phase2(b, 1) "
		     << "from integrate_list [1]"
		     << " at time " << b->get_system_time() << endl;
		pp3("(21,100021)");
		pp3("(23,100023)");
#endif

		initialize_system_phase2(b, 1);		// default set_dt
		b->reconstruct_perturbed_list();

#if 0
		cerr << "after initialize_system_phase2(b, 1):" << endl;
		pp3("(21,100021)");
		pp3("(23,100023)");
#endif

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

#ifdef CPU_COUNTERS
		cpu_prev = cpu;
		kc->cpu_time_tree_check += (cpu = cpu_time()) - cpu_prev;
#endif

		return return_fac* steps;
	    }
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: integrate_list " << 8 << endl << flush;
    }
#endif

    if (reinitialize) {

	kira_synchronize_tree(b);
	steps += b->n_leaves();

#if 0
	cerr << "call initialize_system_phase2(b, 2) from integrate_list [2]"
	     << " at time " << b->get_system_time() << endl;
	pp3("(21,100021)");
	pp3("(23,100023)");
#endif

	initialize_system_phase2(b, 2, 2);	// always set dt
	b->reconstruct_perturbed_list();

#if 0
	cerr << "after initialize_system_phase2(b, 2):" << endl;
	pp3("(21,100021)");
	pp3("(23,100023)");
#endif

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

		    kc->tree_change++;

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

		    // PRC(adjust_tree); PRC(tree_changed); PRL(pert);

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
			     << kc->step_top_single
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
			    if (!next_nodes[j]->is_valid()) {

// 				cerr << "next_nodes[" <<j<< "] = "
// 				     << next_nodes[j]
// 				     << " is no longer valid" << endl;

				next_nodes[j] = NULL;
			    } else
				for (int k = 0; k < n_list; k++)
				    if (next_nodes[j] == cm_list[k]) {
					next_nodes[j] = NULL;
					break;
				    }
			}
		    }
		    if (cm_list) delete [] cm_list;

#ifdef CPU_COUNTERS
		    cpu_prev = cpu;
		    kc->cpu_time_tree_check += (cpu = cpu_time()) - cpu_prev;
#endif

		    return return_fac*steps;
					// NOTE: we currently return after
					//	 the FIRST tree rearrangement,
					// 	 so we can only have one
					//	 restructuring per block time
					//	 step.
		}
		if (cm_list) delete [] cm_list;
	    }
	}
    }

#ifdef CPU_COUNTERS
    cpu_prev = cpu;
    kc->cpu_time_tree_check += (cpu = cpu_time()) - cpu_prev;
#endif

    return return_fac*steps;
}
