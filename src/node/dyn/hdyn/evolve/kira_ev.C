
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// Functions associated with force evaluations within kira.
//
// Externally visible functions:
//
//	void calculate_acc_and_jerk_for_list
//	void calculate_acc_and_jerk_on_all_top_level_nodes
//	void kira_synchronize_tree
//	void initialize_system_phase2
//	void clean_up_kira_ev

// Substantially rewrote code and removed timing and some other
// debugging #ifdefs -- too confusing, and no longer relevant
// after the rewrite.
//					      Steve 4/03, 8/03
//
// Major reorganization to remove compile-time GRAPE selection.
//					      Steve 6/04

#include "hdyn.h"
#include "kira_timing.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira_ev
#   undef T_DEBUG
#endif


//-------------------------------------------------------------------------
//
// NOTE:  Currently, kira_calculate_top_level_acc_and_jerk() is a simple
// compile-time switch between the various (GRAPE and non-GRAPE) methods
// available to compute the acc and jerk (due to internal forces only)
// on the particles listed in next_nodes[].
//
// Alternatively, it could be replaced by a function pointer leading
// directly to the functions involved, set at startup based on the
// configuration of the system.  Implemented and then removed, but
// most of the commented-out code can still be found...
//
// The function declaration is:
//
//	int kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
//						  int n_next,
//					  	  xreal time,
//					  	  bool restart_grape);
//
// where restart_grape may be ignored or used for other purposes in future
// applications.  It indicates a change in top-level tree structure that
// requires an internal reset in the GRAPE arrays.  However, the GRAPE-6
// code now sets an internal flag to convey the information, so only the
// GRAPE-4 code actually uses this variable.  The functions's return
// value is simply the number of top-level nodes in the system.
//
// Currently, the list of possible functions that can be used is
//
//	top_level_acc_and_jerk_for_list()	below
//	grape4_calculate_acc_and_jerk()		in hdyn_grape4.C
//	grape6_calculate_acc_and_jerk()		in hdyn_grape6.C
//
// Note the acc_function_ptr typedef in hdyn.h.
//
// Currently, the only call comes from calculate_acc_and_jerk_for_list(),
// and restart_grape is always set false there after the call.
//
//-------------------------------------------------------------------------


local inline int _kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
							int n_next,
							xreal time,
							bool restart_grape)
{
    if (n_next <= 0) return 0;

    int n_top = 0;
    unsigned int config = next_nodes[0]->get_config();

    switch (config) {

	case 0:	n_top = top_level_acc_and_jerk_for_list(next_nodes, n_next,
							time);
		break;

	case 1:	n_top = grape4_calculate_acc_and_jerk(next_nodes, n_next,
						      time, restart_grape);
		break;

	case 2:	n_top = grape6_calculate_acc_and_jerk(next_nodes, n_next,
						      time, restart_grape);
		break;
    }

    return n_top;
}

// (Use the inline version within this file...)

int kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
					  int n_next,
					  xreal time,
					  bool restart_grape)
{
    return _kira_calculate_top_level_acc_and_jerk(next_nodes, n_next,
						  time, restart_grape);
}

int top_level_acc_and_jerk_for_list(hdyn **next_nodes,
				    int n_next,
				    xreal time)
{
    int n_top = 0;

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	if (bi->is_top_level_node())
	    n_top = bi->top_level_node_real_force_calculation(); // in hdyn_ev.C
    }

    return n_top;
}



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

	    // Prologue operations.  DO NOT clear the interaction here,
	    // as acc and jerk may be needed to update the GRAPE before
	    // the new force calculation.

	    // bi->clear_interaction();

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

    // Retain restart_grape for compatibility with other functions, but
    // in the GRAPE-6 version we actually transmit the information via
    // the static restart_grape flag.  The reason is that (as of 3/05)
    // there now are two functions which may need to reset the internal
    // data structures, and it is too complicated to carry the flag
    // through the many levels before we get to the new call.  The flag
    // will be unset as soon as it is tested and action taken.
    //
    // In general, the use of static flags is probably a better way to
    // handle this sort of communication (Steve, 3/05).

    // Synchronize the old and new flags.  Retain the old version until we
    // are sure the new one works, then we can remove all restart_grape
    // references from the GRAPE version of calculate_acc_and_jerk()...
    // 							(Steve 3/05)

    if (tree_changed) restart_grape = true;
    if (restart_grape) b->set_restart_grape_flag();

    // Complete the top-level force calculation (now top-level nodes
    // run from 0 to n_top-1):

    if (n_top > 0) {

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 3
		 << endl << flush;
	}
#endif

	int n_force = 1;

	// Calculate top-level forces.  Note new return value from
	// kira_calculate_top_level_acc_and_jerk() (Steve, 4/03).

	// Top_level force calculation actually uses
	//
	//	grape*_calculate_acc_and_jerk()
	//
	// or
	//
	//	top_level_node_real_force_calculation(), via
	//	top_level_acc_and_jerk_for list()
	//
	// as appropriate.

	if (!exact && !ignore_internal) {


// 	    acc_function_ptr get_acc_and_jerk
// 		= b->get_kira_calculate_top_level_acc_and_jerk();
//
// 	    n_force = get_acc_and_jerk(next_nodes,
// 				       n_top, time,
// 				       restart_grape);


	    n_force = _kira_calculate_top_level_acc_and_jerk(next_nodes,
							     n_top, time,
							     restart_grape);

	    // Note that we now clear restart_grape explicitly here.

	    restart_grape = false;
	    b->clear_restart_grape_flag();	// should be unnecessary
	    
	}

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



void kira_synchronize_tree(hdyn *b,
			   bool sync_low_level)		// default = false
{
    // GRAPE replacement for synchronize_tree().  Synchronize all
    // top-level nodes.  Called from integrate_list() in kira_ev.C and
    // hdyn::merge_nodes().  For now, at least, the entire algorithm
    // from hdyn_ev.C is repeated here.

    if (b->has_grape()) {

	// Code is similar to that in integrate_list(), but only top-level
	// nodes are considered and we don't check for errors in function
	// correct_and_update.  Possibly should merge this with (part of)
	// integrate_list() and drop synchronize_tree() completely.
	//						     (Steve, 1/02)

	// Make a list of top-level nodes in need of synchronization.
	// Generally interested in recomputation of acc and jerk, so
	// probably don't need to treat low-level nodes.  Default is
	// not to touch them.  Note that, even if sync_low_level is
	// true, we still won't synchronize unperturbed binaries.

	xreal sys_t = b->get_system_time();

	cerr << endl
	     << "synchronizing tree using GRAPE at time " << sys_t
	     << endl << flush;

	// Note: no need to set time steps here, and in fact GRAPE will
	// complain if j-particle times and timesteps are not consistent.

	int n_next = 0, n_top = 0;
	for_all_daughters(hdyn, b, bi) {
	    n_top++;
	    if (bi->get_time() < sys_t) {
		// unnecessary and bad:
		// bi->set_timestep(sys_t - bi->get_time());
		n_next++;
	    }
	}

	hdyn **next_nodes = new hdynptr[n_next];
	n_next = 0;
	for_all_daughters(hdyn, b, bi)
	    if (bi->get_time() < sys_t) next_nodes[n_next++] = bi;

	// Integrate all particles on the list.  Start by computing forces.
	// (Assume exact = false and ignore_internal = false.)

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    predict_loworder_all(bi, sys_t);
	    bi->clear_interaction();
	    bi->top_level_node_prologue_for_force_calculation(false);
	}

	bool restart_grape = false;

//	b->get_kira_calculate_top_level_acc_and_jerk()(next_nodes, n_next,
//						       sys_t, restart_grape);

	kira_calculate_top_level_acc_and_jerk(next_nodes, n_next,
					      sys_t, restart_grape);

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    bi->inc_direct_force(n_top-1);
	    bi->top_level_node_epilogue_force_calculation();
	    bi->inc_steps();
	}

	// Complete calculation of accs and jerks by correcting for C.M.
	// interactions and applying external fields.

	kira_options *ko = b->get_kira_options();

	for (int i = 0; i < n_next; i++)
	    next_nodes[i]->set_on_integration_list();

	if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list) {
	    bool reset = false;
	    correct_acc_and_jerk(b, reset);		// old version
	} else
	    correct_acc_and_jerk(next_nodes, n_next);	// new version

	for (int i = 0; i < n_next; i++)
	    next_nodes[i]->clear_on_integration_list();

	if (b->get_external_field() > 0) {

	    // Add external forces.

	    for (int i = 0; i < n_next; i++) {
		hdyn *bi = next_nodes[i];
		real pot;
		vec acc, jerk;
		get_external_acc(bi, bi->get_pred_pos(), bi->get_pred_vel(),
				 pot, acc, jerk);
		bi->inc_pot(pot);
		bi->inc_acc(acc);
		bi->inc_jerk(jerk);
	    }
	}

	// Apply corrector and redetermine timesteps.

	real st = sys_t;

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    bi->correct_and_update();
	    bi->init_pred();
	    bi->store_old_force();

	    // As in synchronize_tree, make sure time step is consistent with
	    // system_time (= time).

	    real timestep = bi->get_timestep();

	    int iter = 0;
	    while (fmod(st, timestep) != 0) {
		if (iter++ > 30) break;
		timestep *= 0.5;
	    }

	    if (iter > 20) {
		cerr << "kira_synchronize_tree: " << bi->format_label() << " ";
		PRL(iter);
		PRI(4); PRC(fmod(st, timestep)); PRL(timestep);
	    }

	    bi->set_timestep(timestep);
	}

	delete [] next_nodes;

	// Top-level nodes are all synchronized.  Deal with low-level nodes.  

	if (sync_low_level) {
	    for_all_daughters(hdyn, b, bi) {
		hdyn *od = bi->get_oldest_daughter();
		if (od && !od->get_kepler()) {

		    // Synchronizing od will also synchronize its sister,
		    // but call synchronize_tree explicitly to ensure that
		    // substructure is properly handled.

		    synchronize_tree(od);

		    hdyn *yd = od->get_younger_sister();
		    if (yd) synchronize_tree(yd);		// else error?

		}
	    }
	}

	cerr << endl
	     << "end of synchronization"
	     << endl << flush;

    } else {

	cerr << endl
	     << "synchronizing tree without GRAPE at time "
	     << b->get_system_time() << endl;

	synchronize_tree(b);
    }
}



// initialize_system_phase2:  Calculate acc, jerk, timestep, etc for all nodes.
//
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

	    if (set_dt || ((real)bb->get_time() <= 0
			   || bb->get_timestep() <= 0) ) {

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

    if (set_dt && (real)time <= 0) {
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
