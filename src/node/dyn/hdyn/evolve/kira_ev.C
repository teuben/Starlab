
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

// Substantially rewrote code and removed timing and some other
// debugging #ifdefs -- too confusing, and no longer relevant
// after the rewrite.
//					      Steve 4/03, 8/03

#include "hdyn.h"
#include "kira_timing.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira_ev
#   undef T_DEBUG
#endif

int calculate_acc_and_jerk_for_list(hdyn **next_nodes,
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

    if (!next_nodes[0] || !next_nodes[0]->is_valid()) return 0;	  // unnecessary

    hdyn *b = next_nodes[0]->get_root();

    kira_counters *kc = b->get_kira_counters();
    kira_diag *kd = b->get_kira_diag();
    kira_options *ko = b->get_kira_options();

    bool ignore_internal = b->get_ignore_internal();

    xreal sys_t = b->get_system_time();

    // Note that time and system_time should be the same...

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 1 << endl;
	int p = cerr.precision(HIGH_PRECISION);
	PRI(7); PRC(n_next); PRC(exact); PRC(time); PRL(time-sys_t);
	cerr.precision(p);
    }
#endif

    // Explicitly split the calculation into top-level and low-level
    // nodes (Steve, 4/03).  Top-level nodes run from 0 to n_top-1;
    // low-level nodes from n_top to n_next-1.
    //
    // The first loop here does the splitting, at the same time
    // performing the first part of the force calculation.

    bool print = kd->kira_ev;
    int n_top = 0;

    for (int i = 0; i < n_next; i++) {

	hdyn *bi = next_nodes[i];
	bool top = false;

	if (bi->get_parent() == b) {

	    // This is a top-level node.

	    top = true;

	    predict_loworder_all(bi, sys_t);

	    // Prologue operations:

	    bi->clear_interaction();

	    // Note that top_level_node_prologue_for_force_calculation()
	    // does the *entire* calculation in the case exact = true.
	    // For exact = false it does almost nothing...

	    if (!ignore_internal)
		bi->top_level_node_prologue_for_force_calculation(exact);

	} else {

	    // This is a low-level node.
	    // Predict the entire clump containing node bi (Steve, 8/98):

	    if (!ignore_internal)
		predict_loworder_all(bi->get_top_level_node(), sys_t);
	}

	if (ignore_internal && (!top || !exact)) {

	    // Set flags in case of no internal forces...

	    bi->set_nn(bi);
	    bi->set_coll(bi);
	    bi->set_d_nn_sq(VERY_LARGE_NUMBER);
	}

	bi->inc_steps();

	// All initial operations on bi are complete.
	// Reorder the list if necessary.

	if (top) {
	    next_nodes[i] = next_nodes[n_top];
	    next_nodes[n_top++] = bi;
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 2
	     << endl << flush;
    }
#endif

    // Complete the top-level force calculation (now top-level nodes
    // run from 0 to n_top-1):

    if (n_top > 0) {

	if (tree_changed)
	    restart_grape = true;

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 3
		 << endl << flush;
	}
#endif

	int n_force = 1;

	// Calculate top-level forces.  Note new return value from
	// kira_calculate_top_level_acc_and_jerk() (Steve, 4/03).

	// Top_level force calculation uses
	//
	//	grape_calculate_acc_and_jerk()
	//
	// or
	//
	//	top_level_node_real_force_calculation()
	//
	// as appropriate.

	if (!exact && !ignore_internal)
	    n_force = kira_calculate_top_level_acc_and_jerk(next_nodes,
							    n_top, time,
							    restart_grape);

//	if (sys_t >= 44.15329 && sys_t <= 44.1533) {
//	    cerr << "after top-level acc_and_jerk" << endl;
//	    pp3(next_nodes[0]->get_top_level_node());
//	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 4
		 << endl << flush;
	}
#endif

	// Complete calculation of top-level accs and jerks by correcting
	// for C.M. interactions.

	if (!exact) {

	    // Note: correct_acc_and_jerk() now checks list membership to
	    // determine if correction is needed.

	    // The new version of correct_acc_and_jerk appears to work, but
	    // retain the possibility of reverting to the old version until we
	    // are sure there are no problems.  Results of the two versions
	    // are similar, but *not* identical.

	    // On_integration_list flags are for use by correct_acc_and_jerk(),
	    // to ensure that we only correct interactions between objects on
	    // the list.  Should be unset immediately on return, but see the
	    // note below.

	    n_force--;

	    // Epilogue operations.

	    for (int i = 0; i < n_top; i++) {

		hdyn *bi = next_nodes[i];

		// Epilogue force calculation mostly performs CM corrections
		// and cleans up the perturber list.

		bi->inc_direct_force(n_force);		// direct force counter
		bi->top_level_node_epilogue_force_calculation();
		bi->set_on_integration_list();
	    }

	    if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list)

		correct_acc_and_jerk(b,			// old version
				     reset_force_correction);
	    else

		correct_acc_and_jerk(next_nodes,	// new version
				     n_top);

	    // Don't unset the "on_integration_list" flags here, because
	    // the loop through memory may cost more than it is worth.
	    // *** Must unset these flags in the calling function. ***

	    // for (int i = 0; i < n_top; i++)
	    //     next_nodes[i]->clear_on_integration_list();

	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 5
		 << endl << flush;
	}
#endif

	// Add external forces, if any, to top-level nodes.

	if (b->get_external_field() > 0) {

	    for (int i = 0; i < n_top; i++) {

		hdyn *bb = next_nodes[i];

		real pot;
		vec acc, jerk;
		get_external_acc(bb, bb->get_pred_pos(), bb->get_pred_vel(),
				 pot, acc, jerk);
		bb->inc_pot(pot);
		bb->inc_acc(acc);
		bb->inc_jerk(jerk);
	    }
	}
    }

    // Everything is done for top-level nodes.  Complete the low-level
    // force calculations.  Do this last, so any recomputation of
    // top-level perturber lists has already occurred.

    if (!ignore_internal) {

	for (int i = n_top; i < n_next; i++) {

	    hdyn *bi = next_nodes[i];
	    if (!bi->get_kepler()) {

		// Compute the forces (perturbed motion).

		if (print) {
		    cerr << "\nComputing force on low-level node "
			 << bi->format_label() << endl;
		    pp3(bi);
		    PRC(bi->get_system_time());
		    PRC(bi->get_time()); PRL(bi->get_t_pred());
		    PRL(bi->get_pred_pos());
		    PRL(bi->get_pred_vel());
		}

		hdyn *sister = bi->get_younger_sister();
		if (!sister) {
		    sister = bi->get_elder_sister();
		    if (!sister) continue;		  // really an error...
		}

		bi->clear_interaction();

		// Doing these steps here allows us to bypass
		// calculate_acc_and_jerk and go directly to
		// calculate_acc_and_jerk_on_low_level_node.

		bi->set_d_coll_sq(VERY_LARGE_NUMBER);
		bi->set_coll(NULL);
		if (sister) sister->set_d_coll_sq(VERY_LARGE_NUMBER);

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
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 6
	     << endl << flush;
    }
#endif

//   if (sys_t >= 44.15329 && sys_t <= 44.1533) {
//     cerr << "after low-level acc_and_jerk" << endl;
//     pp3(next_nodes[0]->get_top_level_node());
//   }

    return n_top;
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

    calculate_acc_and_jerk_for_list(list, i_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    false,	// called before tree changes
				    reset_force_correction, // obsolete
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

    calculate_acc_and_jerk_for_list(list, i_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    true,
				    reset_force_correction, // no longer used
				    restart_grape);
    delete [] list;
}



// initialize_system_phase2:  Calculate acc, jerk, timestep, etc for all nodes.

// NOTE: System should be synchronized prior to calling this function.

// Static data:

static int work_size = 0;
static hdyn ** nodes = NULL;
static int nnodes = -1;

// Allow possibility of cleaning up if necessary:

void clean_up_kira_ev() {if (nodes) delete [] nodes;}

void initialize_system_phase2(hdyn *b,
			      int call_id,	// default = 0
			      int set_dt)	// 0 ==> set only if zero
						// 1 ==> set with limit (def)
						// 2 ==> always set
{
    // cerr << "initialize_system_phase2: "; PRC(call_id), PRL(set_dt);
    dbg_message("initialize_system_phase2", b);

    if (!b->is_root())
	err_exit("initialize_system_phase2 called with non-root");

    xreal time = b->get_system_time();

    int n = 0;
    for_all_nodes(hdyn, b, bb) n++;

    if  (work_size < n) {
	if (!nodes) delete [] nodes;
	work_size = n + 10;
	nodes = new hdynptr[work_size];
    }

    for_all_nodes(hdyn, b, bb) {
	if (!bb->get_kepler() && (bb->get_unperturbed_timestep() > 0)) {

	    // When is this necessary? (SLWM 3/98)

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl
		 << "initialize_system_phase2: "
		 << "creating kepler for unperturbed binary"
		 << endl
		 << "    " << bb->get_parent()->format_label()
		 << " at system time " << time
		 << "  call_id = " << call_id
		 << endl;
	    cerr.precision(p);

	    bb->update_kepler_from_hdyn();
	}
    }

    // Make a list of all nodes except unperturbed binary components.

    n = 0;
    for_all_nodes(hdyn, b, bb) {
	if ((bb != b) && (!bb->get_kepler())) {
	    nodes[n] = bb;
	    n++ ;
	}
    }

    bool tree_changed = true;
    bool reset_force_correction = true;	// no longer used
    bool restart_grape = true;
    bool exact = false;

    calculate_acc_and_jerk_for_list(nodes, n, time,
				    exact,
				    tree_changed,
				    reset_force_correction,  // no longer used
				    restart_grape);

    // (All perturber lists have been redetermined...)

    real min_dt = VERY_LARGE_NUMBER;

    for_all_nodes(hdyn, b, bb) {
	if ((bb != b) && (!bb->get_kepler())) {

	    bb->store_old_force();

	    if (set_dt || (bb->get_time() <= 0 || bb->get_timestep() <= 0) ) {

		real dtlim = bb->get_timestep()/2;
		bb->set_first_timestep();
		if (set_dt == 1 && bb->get_timestep() < dtlim)
		    bb->set_timestep(dtlim);

		if (bb->get_time() + bb->get_timestep() < time) {

		    // Note from Steve and Simon, Feb 26, 1999:

		    // This bug can only occur if the function is
		    // improperly called.

		    cerr << endl << "warning: initialize_system_phase2: "
			 << "time will go backwards!" << endl;

		    int p = cerr.precision(HIGH_PRECISION);
		    PRL(bb->format_label());
		    PRC(time); PRL(bb->get_time() + bb->get_timestep());
		    cerr.precision(p);

		    cerr << "function should not be called with "
			 << "unsynchronized system." << endl << endl;
		    pp3(bb);

		    // Fix would be to force timestep to go past system
		    // time, but this is not in general possible while
		    // maintaining the block step structure.  Could
		    // simply terminate here, but for now we let the
		    // code die in kira.C.
		}
	    }

	    min_dt = Starlab::min(min_dt, bb->get_timestep());
	}
    }

    // Check for perturbed binaries in the input data...

    if (set_dt && time <= 0) {
	for_all_nodes(hdyn, b, bb) {
	    if ((bb != b) && (!bb->get_kepler())) {
		if (bb->is_parent()
		    && bb->get_oldest_daughter()
			 ->get_perturbation_squared() > 1)
		    bb->set_timestep(min_dt);
	    }
	}
    }

    // cerr << "leaving initialize_system_phase2()\n\n";
}
