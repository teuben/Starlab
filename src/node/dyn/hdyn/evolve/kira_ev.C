
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions associated with integration within kira.
//
// Externally visible functions:
//
//	void calculate_acc_and_jerk_for_list
//	void calculate_acc_and_jerk_on_all_top_level_nodes
//	void initialize_system_phase2
//	void clean_up_kira_ev

#include "hdyn.h"
#include "kira_timing.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira_ev
#   undef T_DEBUG
#endif

int calculate_acc_and_jerk_for_list(hdyn *b,
				    hdyn **next_nodes,
				    int  n_next,
				    xreal time,
				    bool exact,
				    bool tree_changed,
				    bool &reset_force_correction, // no longer
								  // used...
				    bool &restart_grape)
{
    // Note that this function uses the same function calls as
    // calculate_acc_and_jerk_on_top_level_node, but in a different
    // order: all prologue calls are performed first, then all
    // top-level (i.e. GRAPE) force calculations, then all epilogue
    // calls.

    kira_counters *kc = b->get_kira_counters();
    kira_diag *kd = next_nodes[0]->get_kira_diag();
    kira_options *ko = next_nodes[0]->get_kira_options();

#ifdef CPU_COUNTERS
    real cpu = cpu_time(), cpu_prev;
#endif

    // Return the number of top-level nodes on the list.

    int n_list_top_level = 0;

    bool ignore_internal = next_nodes[0]->get_ignore_internal();

#ifdef TIME_INTERNAL

    // Code to time specific force-calculation operations:

    real cpu0, cpu1, cpu2, cpu3, cpu4;
    int kmax = 1;

#if 0
	for (int k = 0; k < n_next; k++) {
	    hdyn* bi = next_nodes[k];

	    if (bi->name_is("13a")
		&& bi->get_kepler() == NULL
		&& bi->get_time() > 6.6
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

    for (int k = 0; k < kmax; k++) {

#endif

    xreal sys_t = next_nodes[0]->get_system_time();

    // Note that time and system_time should be the same...

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 1 << endl;
	int p = cerr.precision(HIGH_PRECISION);
	PRI(7); PRC(n_next); PRC(exact); PRC(time); PRL(time - sys_t);
	cerr.precision(p);
    }
#endif

    for (int i = 0; i < n_next; i++) {

	hdyn *bi = next_nodes[i];

	// Predict the entire clump containing node bi (Steve, 8/98):

	predict_loworder_all(bi->get_top_level_node(), sys_t);

	if (bi->is_top_level_node()) {
	    bi->clear_interaction();
	    bi->top_level_node_prologue_for_force_calculation(exact);
	    n_list_top_level++;
	}
    }

#ifdef CPU_COUNTERS
    cpu_prev = cpu;
    kc->cpu_time_predict += (cpu = cpu_time()) - cpu_prev;
#endif

#ifdef TIME_INTERNAL
    }
    if (kmax > 1) cpu1 = cpu_time();
#endif

    if (tree_changed)
	restart_grape = true;

#ifdef TIME_INTERNAL
    for (int k = 0; k < kmax; k++) {
#endif

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 2
	     << endl << flush;
    }
#endif

    if (!exact) {

	// Top-level forces.

	if (!ignore_internal)
	    kira_calculate_top_level_acc_and_jerk(next_nodes, n_next,
						  time, restart_grape);

	// (Uses grape_calculate_acc_and_jerk or
	//  top_level_node_real_force_calculation, as appropriate.)

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 3
		 << endl << flush;
	}
#endif

	int n = get_n_top_level();

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];

	    if (ignore_internal) {
		bi->set_nn(bi);
		bi->set_coll(bi);
		bi->set_d_nn_sq(VERY_LARGE_NUMBER);
	    }

	    if (bi->is_top_level_node()) {

		bi->inc_direct_force(n-1);	// direct force counter

		bi->top_level_node_epilogue_force_calculation();
	    }
	}

#ifdef CPU_COUNTERS
	cpu_prev = cpu;
	kc->cpu_time_top_level_force += (cpu = cpu_time()) - cpu_prev;
#endif
    }

#ifdef TIME_INTERNAL
    }
    if (kmax > 1) cpu2 = cpu_time();
#endif

#ifdef TIME_INTERNAL
    for (int k = 0; k < kmax; k++) {
#endif

    bool print = kd->kira_ev;

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 4
		 << endl << flush;
#if 0
	    pp3("(21,100021)");
	    pp3("(23,100023)");
	    print = true;
#endif
	}
	int ppp = cerr.precision(HIGH_PRECISION);
#endif

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];

	// Compute forces on low-level nodes.

	if (bi->is_low_level_node() && bi->get_kepler() == NULL) {

	    if (print) {
	        cerr << "\nComputing force on low-level node "
		     << bi->format_label() << endl;
		pp3(bi);
		PRC(bi->get_system_time());
		PRC(bi->get_time()); PRL(bi->get_t_pred());
		PRL(bi->get_pred_pos());
		PRL(bi->get_pred_vel());
	    }

#ifdef KIRA_DEBUG
	    hdyn *sister = bi->get_binary_sister();
#else
	    hdyn *sister = bi->get_younger_sister();	// assume binary tree
	    if (!sister) {
		sister = bi->get_elder_sister();
		if (!sister)
		    continue;				// really an error...
	    }
#endif

	    bi->clear_interaction();

	    // Doing these steps here allows us to bypass calculate_acc_and_jerk
	    // and go directly to calculate_acc_and_jerk_on_low_level_node.

	    bi->set_d_coll_sq(VERY_LARGE_NUMBER);
	    bi->set_coll(NULL);
	    sister->set_d_coll_sq(VERY_LARGE_NUMBER);
	    bi->calculate_acc_and_jerk_on_low_level_node();

	    kc->pert_step++;

	    if (print) {
		cerr << "after..."<<endl;
		pp3(bi);
		PRL(bi->get_acc());
		PRL(bi->get_binary_sister()->get_acc());
		PRL(bi->get_pos());
		PRL(bi->get_binary_sister()->get_pos());
	    }
	}

	// Update all step counters (direct force counters are
	// incremented above; indirect force counters are handled
	// in hdyn_ev.C)

#ifdef TIME_INTERNAL
	if (k == 0)
#endif

	bi->inc_steps();
    }

#ifdef CPU_COUNTERS
    cpu_prev = cpu;
    kc->cpu_time_low_level_force += (cpu = cpu_time()) - cpu_prev;
#endif

#ifdef TIME_INTERNAL
    }
    if (kmax > 1) cpu3 = cpu_time();
#endif

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 5
	     << endl << flush;
#if 0
	pp3("(21,100021)");
	pp3("(23,100023)");
#endif
    }
    cerr.precision(ppp);
#endif

    // Complete calculation of accs and jerks by correcting for C.M.
    // interactions.

#ifdef TIME_INTERNAL
    for (int k = 0; k < kmax; k++) {
#endif

    if (!exact) {

	// Note: correct_acc_and_jerk() now checks list membership to
	// determine if correction is needed.

	// The new version of correct_acc_and_jerk appears to work, but
	// retain the possibility of reverting to the old version until we
	// are sure there are no problems.  Results of the two versions
	// are similar, but *not* identical.

	// On_integration_list flags are for use by correct_acc_and_jerk(),
	// to ensure that we only correct interactions between objects on
	// the list.

	for (int i = 0; i < n_next; i++)
	    next_nodes[i]->set_on_integration_list();

	if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list)

	    correct_acc_and_jerk(b,			// old version
				 reset_force_correction);
	else

	    correct_acc_and_jerk(next_nodes, n_next);	// new version

	// Cautiously unset the "on_integration_list" flags...

	for (int i = 0; i < n_next; i++) {
	    hdyn *n = next_nodes[i];
	    if (n && n->is_valid())
		n->clear_on_integration_list();
	}

#ifdef CPU_COUNTERS
	cpu_prev = cpu;
	kc->cpu_time_top_level_force += (cpu = cpu_time()) - cpu_prev;
#endif
    }

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 6
		 << endl << flush;
	}
#endif

#ifdef TIME_INTERNAL
    }
    if (kmax > 1) cpu4 = cpu_time();

    for (int k = 0; k < kmax; k++) {
#endif

    if (b->get_external_field() > 0) {

        // Add external forces.

        for (int i = 0; i < n_next; i++)
	    if (next_nodes[i]->is_top_level_node()) {
		hdyn *bb = next_nodes[i];
		real pot;
		vector acc, jerk;
	        get_external_acc(bb, bb->get_pred_pos(), bb->get_pred_vel(),
				 pot, acc, jerk);
		bb->inc_pot(pot);
		bb->inc_acc(acc);
		bb->inc_jerk(jerk);
	    }

#ifdef CPU_COUNTERS
	cpu_prev = cpu;
	kc->cpu_time_external_force += (cpu = cpu_time()) - cpu_prev;
#endif
    }

#ifdef TIME_INTERNAL
    }
    if (kmax > 1) {
	cerr << "CPU times:  "
	     << (cpu1-cpu0)/kmax << "  "
	     << (cpu2-cpu1)/kmax << "  "
	     << (cpu3-cpu2)/kmax << "  "
	     << (cpu4-cpu3)/kmax << "  "
	     << (cpu_time()-cpu4)/kmax << endl;
	exit(0);
    }
#endif

    return n_list_top_level;
}

void calculate_acc_and_jerk_on_all_top_level_nodes(hdyn * b)
{
    int n_top = b->n_daughters();
    hdynptr * list = new hdynptr[n_top];

    int i_top = 0;
    for_all_daughters(hdyn, b, bb)
	list[i_top++] = bb;

    bool reset_force_correction = true; 	// (no longer used)
    bool restart_grape = true;

    calculate_acc_and_jerk_for_list(b, list, n_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    false,	// called before tree changes
				    reset_force_correction, // no longer used
				    restart_grape);

    delete [] list;
}

void calculate_acc_and_jerk_on_top_level_binaries(hdyn * b)
{
    int n_top = b->n_daughters();
    hdynptr * list = new hdynptr[n_top];

    int i_top = 0;
    for_all_daughters(hdyn, b, bb)
	if (bb->is_parent()) list[i_top++] = bb;

    bool reset_force_correction = true; 	// (no longer used)
    bool restart_grape = true;

    calculate_acc_and_jerk_for_list(b, list, i_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    true,
				    reset_force_correction, // no longer used
				    restart_grape);

    delete [] list;
}



// initialize_system_phase2:  Calculate acc, jerk, timestep, etc for all nodes.

// System should be synchronized prior to calling this function.

// Static data:

static int work_size = 0;
static hdyn ** nodes = NULL;
static int nnodes = -1;

// Allow possibility of cleaning up if necessary:

void clean_up_kira_ev() {if (nodes) delete [] nodes;}

void initialize_system_phase2(hdyn * b,
			      int call_id,	// default = 0
			      int set_dt)	// 0 ==> set only if zero
						// 1 ==> set with limit (def)
						// 2 ==> always set
{
//    cerr << "initialize_system_phase2: "; PRC(call_id), PRL(set_dt);

    dbg_message("initialize_system_phase2", b);

    if (!b->is_root()) {
	err_exit("initialize_system_phase2 called with non-root");
    }

    int n = 0;
    for_all_nodes(hdyn, b, b1) n++ ;

    if  (work_size < n) {
	if (!nodes) delete [] nodes;
	work_size = n + 10;
	nodes = new hdynptr[work_size];
    }

    // (New variable names necessary here because of DEC C++...)

    for_all_nodes(hdyn, b, b2) {
	if (!b2->get_kepler() && (b2->get_unperturbed_timestep() > 0)) {

	    // When is this necessary? (SLWM 3/98)

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl
		 << "initialize_system_phase2: "
		 << "creating kepler for unperturbed binary"
		 << endl
		 << "    " << b2->get_parent()->format_label()
		 << " at system time " << b->get_system_time()
		 << "  call_id = " << call_id
		 << endl;
	    cerr.precision(p);

	    b2->update_kepler_from_hdyn();
	}
    }

    // Make a list of all nodes except unperturbed binary components.

    n = 0;
    for_all_nodes(hdyn, b, b3) {
	if ((b3 != b) && (!b3->get_kepler())) {
	    nodes[n] = b3;
	    n++ ;
	}
    }

    xreal time = nodes[0]->get_system_time();

    bool tree_changed = true;
    bool reset_force_correction = true;	// no longer used
    bool restart_grape = true;
    bool exact = false;

    calculate_acc_and_jerk_for_list(b, nodes, n, time,
				    exact,
				    tree_changed,
				    reset_force_correction,  // no longer used
				    restart_grape);

    // (All perturber lists have been redetermined...)

    real min_dt = VERY_LARGE_NUMBER;

    for_all_nodes(hdyn, b, b4) {
	if ((b4 != b) && (!b4->get_kepler())) {

	    b4->store_old_force();

	    if (set_dt || (b4->get_time() <= 0 || b4->get_timestep() <= 0) ) {

		real dtlim = b4->get_timestep()/2;
		b4->set_first_timestep();
		if (set_dt == 1 && b4->get_timestep() < dtlim)
		    b4->set_timestep(dtlim);

		if (b4->get_time() + b4->get_timestep()
		    < b->get_system_time()) {

		    // Note from Steve and Simon, 26Feb, 1999:

		    // This bug can only occur if the function is
		    // improperly called.

		    cerr << endl << "warning: initialize_system_phase2: "
			 << "time will go backwards!" << endl;
		    cerr << "function should not be called with "
			 << "unsynchronized system." << endl << endl;
		    pp3(b4);

		    // Fix would be to force timestep to go past system
		    // time, but this is not in general possible while
		    // maintaining the block step structure.  Could
		    // simply terminate here, but for now we let the
		    // code die in kira.C.

		}

	    }

	    min_dt = Starlab::min(min_dt, b4->get_timestep());
	}
    }

    // Check for perturbed binaries in the input data...

    if (set_dt && b->get_system_time() <= 0) {
	for_all_nodes(hdyn, b, b5) {
	    if ((b5 != b) && (!b5->get_kepler())) {
		if (b5->is_parent()
		    && b5->get_oldest_daughter()
			 ->get_perturbation_squared() > 1)
		    b5->set_timestep(min_dt);
	    }
	}
    }

    // cerr << "leaving initialize_system_phase2()\n\n";
}
