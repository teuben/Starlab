
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //           _\|/_
//=======================================================//             /|\ ~

//  hdyn_ev.C: functions related to orbit integration within the hdyn class.
//.............................................................................
//    version 1:  Nov 1994   Piet Hut, Jun Makino, Steve McMillan
//    version 2:
//.............................................................................
//  This file includes:
//    - node synchronization functions
//    - functions related to time steps
//    - function for finding the next node to integrate the orbit of
//    - functions for initialization
//    - predictor and corrector steps for the Hermite integration scheme
//    - functions for calculating acceleration & jerk, for two given nodes
//    - functions for determining which nodes should exert forces on each other
//    - driver function for a one-time-step integration of a node
//.............................................................................

// Externally visible functions:
//
//	void hdyn::synchronize_node
//	void hdyn::set_first_timestep
//	void hdyn::update
//	bool hdyn::correct
//	void hdyn::flat_calculate_acc_and_jerk
//	void hdyn::perturber_acc_and_jerk_on_leaf
//	void hdyn::tree_walk_for_partial_acc_and_jerk_on_leaf
//	void hdyn::calculate_partial_acc_and_jerk_on_leaf
//	void hdyn::calculate_partial_acc_and_jerk
//	void hdyn::create_low_level_perturber_list
//	void hdyn::create_low_level_perturber_lists
//	void hdyn::calculate_acc_and_jerk_on_low_level_node
//	void hdyn::calculate_acc_and_jerk_on_top_level_node
//	void hdyn::top_level_node_prologue_for_force_calculation
//	void hdyn::top_level_node_real_force_calculation
//	void hdyn::top_level_node_epilogue_force_calculation
//	void hdyn::calculate_acc_and_jerk
//	void hdyn::integrate_node
//
//	void update_binary_sister
//	hdyn* node_to_move
//	void get_nodes_to_move
//	void initialize_system_phase1
//	void correct_acc_and_jerk
//	void clean_up_hdyn_ev
//	void synchronize_tree

//  More memo...
//
//  The procedure to determine nn/coll is somewhat sorted out. For nn,
//  the following functions clear nn and d_nn_sq as initialization
//  --- flat_calculate_acc_and_jerk
//  --- perturber_acc_and_jerk
//
//  The following functions do not clear
//  --- calculate_partial_acc_and_jerk(_on_leaf)
//  --- tree_walk_for_...
//  These functions still update nn and d_nn_sq (which are supplied as
//  arguments, not the class members)
//  Thus, if these functions are called to determine nn without preceding
//  call to (e.g.) flat_calculate_acc_and_jerk, d_nn_sq must be initialized
//  to some large number.
//
//  The pointer coll and  d_coll_sq are initialized in the following three
//  functions.
//   --- flat_calculate_acc_and_jerk
//   --- perturber_acc_and_jerk
//   --- calculate_acc_and_jerk (The top-level routine for force calculation)
//   In practice, it might be okay to initialize only at the top-level routine.
//   The pointer coll is different from nn not just by the definition but
//   because of the fact that it is meaningfully defined only for leafs
//   (physical particles) now. Because of this difference, coll is
//   changed always in place, without propagating the pointer to nodes in
//   higher levels.
//   This means that the pointer coll of a leaf can be overwritten during
//   the force calculation of a CM node. This cannot cause any problem, I
//   believe.

// MEMO 8 Aug 1996
//
// collision radius is calculated only for leaves (physical particles)
//
// the actual quantity stored in d_coll_sq is
//
//		rij^2 - (d_i + d_2)^2
//

// MEMO from JM (23 Dec 1995):
//
// Changes over the last few days:
//
// calculate_partial_acc_and_jerk(_on_leaf) has changed
// interface (and algorithm!!!).  It used to take two
// node pointers, top and mask.  Top was assumed to be
// the ancestor of the mask, and mask the ancestor of this.
// It calculated the force from all leaves under top except
// for those under mask.
//
// Now it takes three arguments, top, common and mask.  Common
// must be the common ancestor of top and mask (common = top is OK)
// and the routine calculates the force from the subtree under top.
//
// This change was necessary to implement perturbation correction.

//.........................................................................

// NOTE from Steve to Steve, 7/2/97...
//
// Perturber lists are associated with TOP-LEVEL and parent nodes.
// Kepler pointers are associated with the elder binary component...
// Perturbations   are associated with the binary components...
//
// Perturber lists are recomputed at each top-level CM step
//
// Problem: a binary in a deeply nested hierarchy may be handled
//	    very inefficiently, as it uses the top-level perturber
//	    list, which may be quite large.

// ****  Perturber lists added to all CM nodes (Steve, 8/98).  ****
// ****  Updated at CM step, based on parent perturber list.   ****
// ****  Binary component perturber lists are INDEPENDENT.     ****

// Extra problem: top-level perturber list may be invalid (too big),
// but we can't afford to compute perturbation due to the entire system
// on low-level binaries.  Fix by reusing invalid top-level perturber
// list for low-level binaries (flag: valid_perturbers_low).
// Better (3/05): use GRAPE to compute perturbation, if possible.

// #define ALLOW_LOW_LEVEL_PERTURBERS		true	// now defined
							// in hdyn.h

// Also, no longer resolve an unperturbed binary into components for
// purposes of computing its center of mass motion, and optionally allow
// unperturbed centers of mass to remain on a perturber list without
// resolving them into their components.  (See the ~warning below.)

// #define RESOLVE_UNPERTURBED_PERTURBERS	false	// now defined
							// in hdyn.h

//.........................................................................

#include "hdyn.h"

#ifdef USEMPI
#    include <mpi.h>
#    include "kira_mpi.h"
#endif

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_hdyn_ev
#   undef T_DEBUG
#endif

#include "hdyn_inline.C"

#define INITIAL_STEP_LIMIT	0.0625	// Arbitrary limit on the first step
#define USE_POINT_MASS		true

// General debugging flag:

static bool dbg = false;

extern bool temp_debug_flag;			// defined in kira.C


//=============================================================================
//  node synchronization functions
//=============================================================================

//-----------------------------------------------------------------------------
// update_binary_sister: make sure that the companion of binary component
//			 bi is consistent by updating its time and mirroring
//			 bi's pos, vel, acc, and jerk.
//-----------------------------------------------------------------------------

void update_binary_sister(hdyn * bi)
{
#ifdef KIRA_DEBUG
    hdyn *sister = bi->get_binary_sister();
#else
    hdyn *sister = bi->get_younger_sister();		// assume binary tree
    if (!sister) {
	sister = bi->get_elder_sister();
	if (!sister)
	    return;
//	    err_exit("calculate_acc_and_jerk_on_low_level_node: no sister!");
    }
#endif

    sister->set_time(bi->get_time());
    sister->set_timestep(bi->get_timestep());

    real factor = -bi->get_mass() / sister->get_mass();
    sister->set_pos(factor * bi->get_pos());
    sister->set_vel(factor * bi->get_vel());
    sister->set_acc(factor * bi->get_acc());

    // Note from Steve, 10/12/98: It is apparently possible for an
    // unperturbed binary component to have jerk = Infinity if it
    // has never taken a perturbed step, and this can lead to overflow
    // here.  Check for this explicitly.  (Not entirely clear how the
    // infinity arises...  May be due to GRAPE-4 overflow.)

    real j0 = bi->get_jerk()[0];

    if (j0 < VERY_LARGE_NUMBER && j0 > -VERY_LARGE_NUMBER)
	sister->set_jerk(factor * bi->get_jerk());
    else {
	cerr << "update_binary_sister():  setting jerk = 0 for "
	     << bi->format_label() << " at time " << bi->get_system_time()
	     << endl;
	bi->set_jerk(vec(0,0,0));
	sister->set_jerk(vec(0,0,0));
    }

    sister->set_pot(-factor * bi->get_pot());
    sister->store_old_force();

    sister->set_perturbation_squared(bi->get_perturbation_squared());
}

//-----------------------------------------------------------------------------
// synchronize_node -- force the time of a node to a specified value,
//                     by integrating it forward in time to system_time.
//
//                     No prediction here -- do in calling function.
//
//                     The only difference now between this function and
//		       integrate_node() is that here we make sure explicitly
//		       that the time step is consistent with the current
//		       system time, that is, we don't assume that the node
//		       was necessarily scheduled to advance to this time.
//
//							(Steve, 4/03)
//
//-----------------------------------------------------------------------------

void hdyn::synchronize_node(bool exact)		// default = true, so
						// integrate_node does an
						// exact calculation in
						// updating positions

						// Hmmm...  exact = false
						// may not be very useful...
						// Just added -- maybe cut?
						// 		Steve, 8/03
						// But see opposing note at
						// the end of this file....
{
  bool do_diag = false; 	// time > 0.0492;


#ifdef KIRA_DEBUG
    do_diag = (diag->ev && diag->check_diag(this));
    if (do_diag)
	cerr << "Enter synchronize_node for " << format_label() << endl;
#endif

    if (time == system_time || kep) return;	// don't sync unperturbed
						// binary motion

    if (is_low_level_node()
	&& (get_parent()->get_oldest_daughter() != this)) {
	if (do_diag)
	    cerr << "Low-level node: "<<format_label()<<endl<<flush;
	get_elder_sister()->synchronize_node();
	update_binary_sister(get_elder_sister());
	return;
    }

    // No need to adjust timestep before taking the step -- step is now
    // determined by system_time (Steve, 4/03).

    // real old_timestep = timestep;
    // timestep = system_time - time;
    //
    // if (timestep < 0.0) {
    //
    // if (do_diag)
    //     cerr << "synchronize_node: negative timestep..."
    //  	<< system_time << " " << time << flush << endl;
    //     put_node(this, cerr, options->print_xreal);
    //     exit(-1);
    // }

    if (do_diag) {
	cerr << "synchronize_node for " << format_label()
	     << " at time " << system_time << endl << flush;
	PRL(system_time-time);
	PRC(timestep); PRL(is_low_level_node());
	if (is_low_level_node()) pp3(get_parent());
	cerr << flush;
    }

    // NOTE: the "false" here means that we do NOT attempt to
    //	     synchronize unperturbed binaries...

    // cerr << "synchronize_node " << format_label() << endl;
    // PRC(system_time); PRL(time);

    // Use integrate_node() to do the actual work.  Don't integrate
    // unperturbed binary motion.

    // cerr << "Energy before "; print_recalculated_energies(get_root());

    if (do_diag) PRL(timestep);
    integrate_node(exact,
		   false,		// don't integrate unperturbed binaries
		   false);		// don't force unperturbed binary time
    if (do_diag) PRL(timestep);

    // cerr << "Energy after "; print_recalculated_energies(get_root());

    if (do_diag) {
	cerr <<"After integrate_node "<<endl;
	PRC(timestep); PRC(system_time); PRL(time);
	PRL(kep);
    }

    if (!kep) {

        // Must ensure that the new particle time steps are consistent
        // with system_time (which should now be identical to time).
        // However, as in kira_ev.C, we skip this if the current time
        // is in a very high block number, and rely on later
        // synchronization to fix the problem.

	int kb = get_effective_block(time);
	if (kb < 35) {				// ~ arbitrary limit

	    real old_timestep = timestep;
	    int iter = 0;
	    while (fmod(time, timestep) != 0) {
	      if (iter++ > 30) break;		// also arbitrary
		timestep *= 0.5;
	    }

	    if (do_diag || iter > 20) {
		cerr << "synchronize_node: " << format_label() << " "; PRL(iter);
		int p = cerr.precision(15);
		PRI(4); PRC(time); PRL(system_time);
		PRI(4); PRC(old_timestep); PRL(fmod(time, old_timestep));
		PRI(4); PRC(timestep); PRL(fmod(time, timestep));
		cerr.precision(p);
	    }
	}

	if (is_low_level_node())
	    update_binary_sister(this);
    }

    if (do_diag)
	cerr << "Leave synchronize_node for " << format_label() << endl;
}

//=============================================================================
//  functions related to time steps
//=============================================================================

//-----------------------------------------------------------------------------
// set_first_timestep -- calculate the first timestep using acc, jerk,
//                       and pot (acceleration, jerk, potential).
//                       The expression used here is:
//                                                 ________
//               1              /   | acc |       V | pot |   \         ~
//   timestep = ---  eta * min (   ---------  ,  -----------   )        ~
//               16             \   | jerk |       | acc |    /         ~
//
//  The first expression gives the time scale on which the acceleration is
//  changing significantly.
//  The second expression gives an estimate for the free fall time
//  in the case of near-zero initial velocity.
//-----------------------------------------------------------------------------

void hdyn::set_first_timestep(real additional_step_limit) // default = 0
{
    real eta_init = eta / 16;		// initial accuracy parameter (arbitrary)

    real j2 = jerk * jerk;
    real a2 = acc * acc;

    real step_limit = eta_init;		// another arbitrary limit...

    real dt_adot = VERY_LARGE_NUMBER;	// time step based on rate of change of acc
    real dt_ff   = VERY_LARGE_NUMBER;	// time step based on free-fall time

//    if (j2 > 0)
//	dt_adot = eta_init * sqrt(a2 / j2);
//
//    if (pot < 0 && a2 > 0)
//	dt_ff = eta_init * sqrt(-pot / a2);
//
//    real dt = min(dt_adot, dt_ff);

    if (j2 > 0) dt_adot = a2 / j2;
    if (pot < 0 && a2 > 0) dt_ff = -pot / a2;

    real dt = eta_init * sqrt(Starlab::min(dt_adot, dt_ff));

    // Apply an upper limit to the time step:

    dt = Starlab::min(dt, step_limit);

    real true_limit = Starlab::min(initial_step_limit, step_limit);
    if (additional_step_limit > 0) true_limit = Starlab::min(true_limit,
						    additional_step_limit);

    timestep = adjust_number_to_power(dt, true_limit);

    if ((real)time != 0)
	while (fmod(time, timestep) != 0)
	    timestep *= 0.5;

    xreal tnext = time + timestep;

    if (timestep == 0 || time == tnext) {

	cerr << endl << "set_first_timestep: dt = 0 at time "
	     << get_system_time() << endl;

	PRL(additional_step_limit);
	PRC(dt_adot); PRL(dt_ff);
	PRC(a2); PRC(j2); PRL(pot);
	pp3(this);

	cerr << endl << endl << "System dump:" << endl << endl;

	pp3(get_root());
	exit(0);
    }
}



local inline real local_kepler_step(hdyn *b,
				    real correction_factor = 1)
{
    // Simple "kepler" time step criterion for elder binary component.
    // Used at low level, so minimal checks -- be careful!!

    hdyn* s = b->get_younger_sister();

    if (s) {					// excessively cautious?

	real dist, dtff2, dtv2;

	// Use predicted quantities here, as sister may not have been
	// updated yet.

	dist = abs(b->get_nopred_pos() - s->get_nopred_pos());
	dtff2 = dist*dist*dist / (b->get_parent()->get_mass());
	dtv2  = square(b->get_nopred_pos()) / square(b->get_nopred_vel());

	// PRC(dist); PRC(square(b->get_vel())); PRC(dtff2); PRL(dtv2);

	// Current choice takes ~100 steps for a circular orbit with no
	// correction and eta = 0.1.  Calibrate against the Aarseth step.

	return 0.6 * correction_factor
	    	   * b->get_eta()
		   * sqrt(Starlab::min(dtff2, dtv2));

    } else if (false) {

        // Experimental code, in case s isn't a binary sister someday.

	real dist, dtff2, dtv2, mtot = b->get_mass() + s->get_mass();

	dist = abs(b->get_nopred_pos() - s->get_nopred_pos());
	dtff2 = dist*dist*dist / mtot;
	dtv2  = square(b->get_nopred_pos() - s->get_nopred_pos())
		  / square(b->get_nopred_vel() - s->get_nopred_vel());

	return 0.6 * correction_factor
	    	   * b->get_eta()
		   * sqrt(Starlab::min(dtff2, dtv2));

    } else

	return 0;
}

real kepler_step(hdyn *b,
		 real correction_factor)	// default = 1
{
    // Globally accessible version of the above.

    real keplstep = local_kepler_step(b, correction_factor);

#if 0
    if (keplstep > 0) {

	// Debugging repeats earlier calculations...

	hdyn* s = b->get_younger_sister();
	real dist  = abs(b->get_pos() - s->get_pos());
	real dtff = sqrt(dist*dist*dist / (2*b->get_parent()->get_mass()));
	real dtv  = sqrt(square(b->get_pos()) / square(b->get_vel()));
	real fac = 0.5*correction_factor*b->get_eta();

	cerr << endl << "kepler_step: ";
	PRC(dist); PRC(dtff); PRL(dtv);
	PRC(correction_factor); PRC(fac); PRL(keplstep);

    }
#endif

    return keplstep;
}

//-----------------------------------------------------------------------------
// new_timestep:  Calculate the next timestep following the Aarseth
//                formula (Aarseth 1985), including necessary adjustments
//                for the hierarchical timestep scheme.
//-----------------------------------------------------------------------------

static int mycount = 0;

local inline real new_timestep(hdyn *b,			// this node
			       vec& at3,		// 3rd order term
			       vec& bt2,		// 2nd order term
			       vec& jerk,		// 1st order term
			       vec& acc,		// 0th order term
			       real dt,			// current step
			       xreal time,		// present time
			       real correction_factor,
			       real pert_sq)
{
    // On entry, dt is the actual step taken to reach the current
    // time.  The current particle step is timestep, which usually
    // will be dt.  However, in case of a forced synchronization, dt
    // will be less than timestep.  We must use dt for computation of
    // the derivatives below, but timestep is still the proper
    // starting point for step modification below.

    // N.B. "dt" may really be dtau (if b->get_kappa() > 1).

    if (dt <= 0) {

	// Shouldn't happen now, as this will be caught in
	// correct_and_update().

	cerr << "new_timestep: actual time step <= 0" << endl;
	PRC(b->format_label()); 
	int p = cerr.precision(HIGH_PRECISION);
	PRL(b->get_system_time());
	PRC(time); PRL(b->get_time());
	cerr.precision(p);
	PRC(dt); PRL(b->get_timestep());

	if (b->get_timestep() > 0)
	    return b->get_timestep();
	else
	    exit(0);
    }

    // Too many debugging flags in this file!  This one is strictly local.

    bool local_debug = false;		// (b->get_time() > 0.0492);

    //----------------------------------------------------------------------
    //
    // Determine the particle's new "natural" step.  Normally use the
    // Aarseth criterion, but that may be replaced by a simpler kepler
    // step for lightly perturbed binaries.  Optionally compare the
    // Aarseth and kepler steps for debugging purposes.

    real newstep, altstep;
    bool keplstep = (pert_sq >= 0 && pert_sq < 0.0001
		     && b->is_low_level_node()
		     && b->get_kira_options()->allow_keplstep);

    bool timestep_check = false;

#ifdef KIRA_DEBUG
    timestep_check = (b->get_kira_diag()->timestep_check
			   && b->get_kira_diag()->check_diag(b));
#endif

    if (0 && time > 400.0 && b->name_is("3277")) {
      local_debug = true;
      // keplstep = true;		// force a Kepler step
      timestep_check = true;		// force Kepler/Aarseth comparison
    }

    if (local_debug) {
	cerr << endl;
	PRC(keplstep); PRL(correction_factor);
	PRC(dt); PRL(b->get_timestep());
    }

    if (keplstep) {

	// Use a simple "kepler" time step criterion if b is a binary
	// and the perturbation is fairly small.

	newstep = altstep = local_kepler_step(b, correction_factor);
	if (newstep == 0) keplstep = false;

	if (local_debug) {
	    PRC(keplstep);
	    PRL(newstep);
	}

#if 0
	if (newstep < 1.e-11) { // && name_is(b, "14")) {
	    PRL(11111);
	    PRL(b->format_label());
	    PRL(newstep);
	    PRL(correction_factor);
	    PRL(abs(b->get_pos()));
	}
#endif

    }

    real a2, j2, k2, l2;
    real tmp1, tmp2, tmp3, tmp4;

    if (!keplstep || timestep_check) {

	// Simple criterion:
	// real newstep = b->get_eta() * sqrt(sqrt(a2 / k2));

	// Use the Aarseth criterion here (use dt to determine derivatives).

	real dt2 = dt * dt;
	real dt4 = dt2 * dt2;
	real dt6 = dt4 * dt2;

	// A note on terminology: these scaled quantities are simply the
	// derivatives of the acceleration -- a, j = a', k = a'', l = a'''.
	// Old notation called j2 a1scaled2, k2 a2scaled2, and l2 a3scaled2.
	// I prefer the new names... (Steve, 9/99)

	// real a3scaled2 = 36 * at3 * at3 / dt6;	// = (a''')^2
	// real a2scaled2 = 4 * bt2 * bt2 / dt4;	// = (a'')^2
	// real a1scaled2 = jerk * jerk;		// = (a')^2

	l2 = 36 * at3 * at3 / dt6;			// = (a''')^2
	k2 = 4 * bt2 * bt2 / dt4;			// = (a'')^2
	j2 = jerk * jerk;				// = (a')^2
	a2 = acc * acc;

	real aarsethstep = 0;

	// Not completely clear why this option should be tied to
	// diag->grape...  (But the most likely place for a problem to
	// occur probably is in the GRAPE acc and jerk.)

	if (!b->get_kira_diag()->grape) {

	    // The simple way:

	    aarsethstep = b->get_eta() * sqrt((sqrt(a2 * k2) + j2)
					       / (sqrt(j2 * l2) + k2));

	    // Note that overflow may conceivably occur here or below,
	    // as the higher derivatives may become very large...

#if 0
	    if (b->get_system_time() > 361.1566) {

	      // Handy to be able to check details (e.g. in case of
	      // inconsistencies in external fields).

	      cerr << endl;
	      PRL(b->format_label());
	      PRL(acc);
	      PRL(jerk);
	      PRL(bt2);
	      PRL(at3);
	      PRC(abs(acc)); PRL(abs(jerk));
	      PRC(a2); PRL(j2);
	      PRC(abs(bt2)); PRL(abs(at3));
	      PRC(k2); PRL(l2);
	    }
#endif

	} else {

	    // Break the calculation up into pieces, for debugging purposes.

	    // Numerator:

	    real sqrt1 = 0;
	    tmp1 = a2 * k2;
	    //	if (tmp1 > 0 && tmp1 < VERY_LARGE_NUMBER) sqrt1 = sqrt(tmp1);
	    if (tmp1 > 0) sqrt1 = sqrt(tmp1);
	    tmp2 = sqrt1 + j2;

	    // Denominator:

	    real sqrt3 = 0;
	    tmp3 = j2 * l2;
	    //	if (tmp3 > 0 && tmp3 < VERY_LARGE_NUMBER) sqrt3 = sqrt(tmp3);
	    if (tmp3 > 0) sqrt3 = sqrt(tmp3);
	    tmp4 = sqrt3 + k2;

	    //	if (tmp2 > 0 && tmp2 < VERY_LARGE_NUMBER
	    //	    && tmp4 > 0 && tmp4 < VERY_LARGE_NUMBER)
	    if (tmp2 > 0 && tmp4 > 0)
		aarsethstep = b->get_eta() * sqrt(tmp2/tmp4);

	    // Stop and flag an error if any of the intermediates
	    // exceeds VERY_LARGE_NUMBER.

	    if (tmp1 >= 0 && tmp2 >= 0 && tmp3 >= 0 && tmp4 >= 0
		&& tmp1 < VERY_LARGE_NUMBER
		&& tmp2 < VERY_LARGE_NUMBER
		&& tmp3 < VERY_LARGE_NUMBER
		&& tmp4 < VERY_LARGE_NUMBER) {	// test negation to catch NaN,
						// which always tests false
	    } else {

		cerr.precision(HIGH_PRECISION);
		cerr << endl << "new_timestep:  found errors at time "
		     << b->get_system_time() << endl;

		PRL(acc);
		PRL(jerk);
		PRL(bt2);
		PRL(at3);
		PRC(abs(acc)); PRL(abs(jerk));
		PRC(a2); PRL(j2);
		PRC(abs(bt2)); PRL(abs(at3));
		PRC(k2); PRL(l2);
		PRC(tmp1); PRC(tmp2); PRC(tmp3); PRL(tmp4);

		plot_stars(b);

		cerr << endl;
		pp3(b, cerr, -1);

		cerr << endl << endl
		     << "Top-level node dump:"
		     << endl << endl;

		pp3(b->get_top_level_node());
		exit(0);
	    }

#if 0
	    if (b->get_system_time() > 361.1566) {

	      // Handy to look at timestep details (e.g. inconsistent acc
	      // and jerk in an external potential...).

	      cerr << endl;
	      PRL(b->format_label());
	      PRL(acc);
	      PRL(jerk);
	      PRL(bt2);
	      PRL(at3);
	      PRC(abs(acc)); PRL(abs(jerk));
	      PRC(a2); PRL(j2);
	      PRC(abs(bt2)); PRL(abs(at3));
	      PRC(k2); PRL(l2);
	      PRC(tmp1); PRC(tmp2); PRC(tmp3); PRL(tmp4);

	    }
#endif

	}

	newstep = aarsethstep * correction_factor;

	if (local_debug) {
	  PRC(aarsethstep); PRL(correction_factor); PRL(newstep);
	}

#if 0
	if (temp_debug_flag && b->is_top_level_node()) {
	    PRC(aarsethstep); PRL(correction_factor);
	    PRC(dt); PRC(b->get_timestep()); PRL(newstep);
	}
#endif
    }

    if (keplstep && timestep_check) {

	if (newstep < 1.e-13 || altstep < 1.e-13) {

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl << "in new_timestep:" << endl;

	    PRC(b->format_label());
	    PRC(b->get_posvel()); PRL(pert_sq);
	    PRC(newstep); PRC(altstep); PRL(altstep/newstep);

	    if (mycount > 0 || altstep/newstep > 10) {
		PRL(b->format_label());
		PRC(b->get_system_time()); xprint(b->get_system_time());
		PRC(b->get_time()); xprint(b->get_time());
		PRC(b->get_t_pred()); xprint(b->get_t_pred());
		PRC(dt); PRL(b->get_timestep());
		xreal xdt = b->get_t_pred() - b->get_time();
		PRC(xdt); xprint(xdt);
		PRL(acc);
		PRL(jerk);
		PRL(bt2);
		PRL(at3);
		PRC(abs(acc)); PRL(abs(jerk));
		PRC(abs(bt2)); PRL(abs(at3));
		PRC(a2); PRL(j2);
		PRC(k2); PRL(l2);
		PRC(tmp1); PRL(tmp2);
		PRC(tmp3); PRL(tmp4);
		if (altstep/newstep > 10) mycount++;
		if (mycount > 100) exit(1);
	    }
	    cerr.precision(p);
	}

	if (local_debug) {
	    PRC(b->format_label()); PRC(altstep); PRL(newstep);
	}

	newstep = altstep;	// comment out to retain Aarseth step
    }

    // Reduce the time step of a strongly perturbed binary component.
    // Exponent in the pow() comes from a simple comparison of time
    // scales...

    if (b->is_low_level_node() && b->get_perturbation_squared() > 1)
	newstep *= pow(b->get_perturbation_squared(), -0.25);

    // Arbitrarily reduce low-level steps by some factor,
    // to facilitate timing checks of low-level steps.
    //
    // if (b->is_low_level_node()) newstep /= 32;

    //----------------------------------------------------------------------
    //
    // Now newstep is the "natural" time step for the particle.  Force
    // it into an appropriate scheduling block.  Ordinarily this will
    // mean that the time step should be a power of two, but in some
    // circumstances (usually if the system took a very short step in
    // the past, or synchronization is being forced), it may not be
    // possible to find a reasonable step that is commensurate with
    // the current time.  In that case, we may want to take a
    // non-power-of-two step to restore the scheduling (NOT yet
    // implemented).

    real timestep = b->get_timestep();
    if (timestep == 0) cerr << "new_timestep: particle timestep = 0" << endl;

    // If dt is less than timestep, it means that the current step has
    // been forced, usually by a system synchronization.  Timestep is
    // probably still OK as a starting point for determining the new
    // step (dt could in principle be very short), but cautiously
    // reduce timestep below dt, if possible.  However, don't let the
    // time step decrease by too much.

    // Find the first power of 2 below timestep.  This is redundant if
    // timestep is already a power of 2, but will become necessary if
    // non-power-of-two time steps are allowed.

    int exponent;
    real timestep2 = timestep/(2*frexp(timestep, &exponent));

    if (local_debug) {
      PRC(1001); PRC(dt); PRL(timestep); PRL(timestep2);
      int kb = get_effective_block(time);
      PRC(kb); PRL(pow(2.0,-kb));
    }

    while (dt < timestep2 && timestep2 > 0.01*timestep) timestep2 /= 2;

    if (local_debug) {PRC(1002); PRC(dt); PRL(timestep); PRL(timestep2);}

    // Note: At this stage it is possible that timestep2 is not
    // commensurate with the current time.  DON'T force it (likely to
    // drop the step to an unacceptable value):
    //
    // while (fmod(time, timestep2) != 0) timestep2 /= 2;

    // Choose the final time step.  Start with timestep2, then adjust
    // it toward newstep if necessary.

    real final_step = timestep2;

    if (newstep < timestep2) {

	// Halving the timestep is always OK.  Note that we don't
	// necessarily drop the step below newstep -- we assume that a
	// single halving is sufficient.

	final_step = 0.5 * timestep2;

    } else {

	real t2 = 2 * timestep2;
	if (newstep >= t2) {

	    //To preserve the synchronization, a doubling is OK only
	    // if the current time could be reached by an integral
	    // number of doubled time steps (use of fmod below).

	    if (fmod(time, t2 * b->get_kappa()) == 0.0
		&& t2 * b->get_kappa() <= b->get_step_limit()) {

		// Added by Steve 7/98 to deal with pathological
		// ICs from SPZ...  Do not double if this is a
		// strongly perturbed center of mass node.

		if (b->is_leaf()
		    || b->get_oldest_daughter()
			->get_perturbation_squared() < 1) final_step = t2;
	    }
	}
    }

    // There is no guarantee that final_step as just defined is
    // commensurate with time, so if we are currently in an
    // undesirably low block we will remain there.  Correcting this
    // would require a non-power-of-two step.  If that is desired, do
    // it here (and reinstate the "timestep2" code above).

    if (local_debug) PRL(final_step);
    return final_step;
}

real timestep_correction_factor(hdyn *b)
{
    // Compute a correction factor to reduce the timestep
    // in a close binary, in order to improve energy errors.

    real m = b->get_parent()->get_mass()/b->get_mbar();
    real pot_sq = m*m*b->get_d_min_sq()/square(b->get_pos());

    // Steve 8/98:  1. Local energy error is fifth order; scale
    //		       the step accordingly
    //		    2. Could rewrite to replace pow() if necessary...
    //		       *** see kira_approx.C ***

    // Hmmm...  This pow() seems to fall victim to the strange math.h
    // bug in Red Hat Linux 5.  In that case, use Steve's approximate
    // version instead.

    real correction_factor = 1;

    if (pot_sq > 100)		// threshold is ~arbitrary
        correction_factor = pow(pot_sq/100, -0.1);

    // May be desirable to place limits on the correction...

    if (correction_factor < 0.2) correction_factor = 0.2;

#if 0
    if (b->get_perturbation_squared() < 0.0001 && pot_sq > 1
	&& timestep < 1.e-11) {
	PRC(b->format_label()); PRL(correction_factor);
    }
#endif

    return correction_factor;
}

//-----------------------------------------------------------------------------
// update:  Update the time and the timestep.
//-----------------------------------------------------------------------------

void hdyn::update(vec& bt2, vec& at3)	// pass arguments to
					// avoid recomputation
{
    // Note from Steve (4/03):  We now draw a distinction betwen dt, the
    // step actually taken, and timestep, the particle's natural step.
    // This is inportant in new_timestep, as the derivatives need dt, but
    // the synchronization needs timestep...

    // time += timestep;
    // real dt = timestep;

    bool local_debug = false;	// time > 0.614;

    real dt = system_time - time;
    time = system_time;

    if (slow) {					// wrong if dt != timestep,
	dt = slow->get_dtau();			// but shouldn't be permitted
	slow->inc_tau(dt);
    }

    // vec at3 = 2 * (old_acc - acc) + dt * (old_jerk + jerk);
    // vec bt2 = -3 * (old_acc - acc) - dt * (2 * old_jerk + jerk);

    // Reminder:	at3  =  a''' dt^3 / 6		(a' = j)
    //			bt2  =  a''  dt^2 / 2

#if 0
    //if (time > 361.1566) {
    if (temp_debug_flag && is_top_level_leaf()) {

	// Handy to be able to check details (e.g. in case of
	// inconsistencies in external fields).

	cerr << endl;
	PRC(time); PRL(format_label());
	PRI(4); PRC(index); PRL(dt);
	int p = cerr.precision(10);
	PRI(4); PRL(pred_pos);
	PRI(4); PRL(pos);
	PRI(4); PRL(pos-pred_pos);
	PRI(4); PRL(pred_vel);
	PRI(4); PRL(vel);
	PRI(4); PRL(old_acc);
	PRI(4); PRL(acc);
	PRI(4); PRL(old_jerk);
	cerr.precision(p);
	PRI(4); PRL(jerk);
	PRI(4); PRL((acc-old_acc)/dt);
	PRI(4); PRL(2*bt2/(dt*dt));
	PRI(4); PRL((jerk-old_jerk)/dt);
	PRI(4); PRL(6*at3/(dt*dt*dt));
	PRI(4); PRL(-3 * (old_acc - acc));
	PRI(4); PRL(-dt * (2 * old_jerk + jerk));
	PRI(4); PRL(bt2);
	PRI(4); PRL(2 * (old_acc - acc));
	PRI(4); PRL(dt * (old_jerk + jerk));
	PRI(4); PRL(at3);
    }
#endif

    bool is_low = is_low_level_node();

    // Define correction factor for use in new_timestep.

    real correction_factor = 1;
    if (is_low) correction_factor = timestep_correction_factor(this);

    if (local_debug) {PRC(101); PRL(timestep);}

    real new_dt = new_timestep(this, at3, bt2, jerk, acc, dt, time,
			       correction_factor, perturbation_squared);

    if (local_debug) {PRC(102); PRL(timestep);}

#if 0
    if (temp_debug_flag && is_top_level_node()) {
	cerr << "clearing temp_debug_flag" << endl << endl;
	temp_debug_flag = false;
    }
#endif


    if (new_dt <= 0) {
	cerr << "update:  dt = " << new_dt
	     << " at time " << get_system_time();
	pp3(this,cerr);
    }

    // Computed new_dt will be timestep if slow not set; dtau otherwise.

    if (slow) {
	slow->set_dtau(new_dt);
	timestep = get_kappa() * new_dt;
    } else
	timestep = new_dt;

    // Update posvel = pos*vel (for use in binary handling).

    if (is_low) {

	prev_posvel = posvel;
	posvel = pos * vel;

    } else

	// Finally, store a'' for possible use in prediction of
	// (top-level, for now) nodes.

	k_over_18 = bt2 / (9*dt*dt);

#if 0
	if (new_dt < 1.e-11) {
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl << "in update:" << endl;;
	    PRC(time);
	    PRC(dt); PRL(abs(pos));
	    PRL(pred_pos);
	    PRL(pos);
	    PRL(pred_vel);
	    PRL(vel);
	    PRL(old_acc);
	    PRL(acc);
	    PRL(old_jerk);
	    PRL(jerk);
	    PRL(-3 * (old_acc - acc));
	    PRL(-dt * (2 * old_jerk + jerk));
	    PRL(bt2);
	    PRL(2 * (old_acc - acc));
	    PRL(dt * (old_jerk + jerk));
	    PRL(at3);
	    PRC(new_dt), PRL(timestep);
	    cerr << endl;
	    cerr.precision(p);
	}
#endif

}

//=============================================================================
//  function for finding the next node to integrate
//=============================================================================

//-----------------------------------------------------------------------------
// node_to_move:  Return the node with minimum t+dt.
//-----------------------------------------------------------------------------

hdyn *node_to_move(hdyn * b,
		   real & tmin)
{
    hdyn *bmin = NULL;
    if (b->is_parent())
	for_all_daughters(hdyn, b, d) {
	    hdyn *btmp;
	    real ttmp = VERY_LARGE_NUMBER;
	    btmp = node_to_move(d, ttmp);
	    if (ttmp < tmin) {
		tmin = ttmp;
		bmin = btmp;
	    }
	}

    if (b->get_parent() != NULL)
	if (tmin > b->get_next_time()) {
	    bmin = b;
	    tmin = b->get_next_time();
	}
    return bmin;
}

//-----------------------------------------------------------------------------
// get_nodes_to_move:  Return the node(s) with minimum t+dt.  Use a recursive
//		       loop to traverse the system.
//-----------------------------------------------------------------------------

void get_nodes_to_move(hdyn * b,
		       hdyn * list[],
		       int &nlist,
		       real & tmin)
{
    hdyn *bmin = NULL;

    if (b->get_parent() == NULL) {			// b is root: initialize
	nlist = 0;
	tmin = VERY_LARGE_NUMBER;
    }

    if (b->is_parent())
	for_all_daughters(hdyn, b, d) {
	    get_nodes_to_move(d, list, nlist, tmin);
	}

    if (b->get_parent() != NULL) {

	if (tmin > b->get_next_time()) {
	    nlist = 1;
	    tmin = b->get_next_time();
	    list[0] = b;

	} else if (tmin == b->get_next_time()) {	// NOTE: we are assuming
	    list[nlist] = b;				// that equality is OK
	    nlist++;					// here because all
	    if (b->is_low_level_node()) {		// times are powers of 2
		hdyn *bs = b->get_binary_sister();
		for (int j = nlist - 2; j >= 0; j--) {
		    if (list[j] == bs) {
			nlist--;
			j = -1;
		    }
		}
	    }
	}
    }
}

//=============================================================================
//  functions for initialization
//=============================================================================

//-----------------------------------------------------------------------------
// initialize_system_phase1:  Set current time and "predicted" position
//                            and velocity for all nodes.
//-----------------------------------------------------------------------------

void initialize_system_phase1(hdyn * b,  xreal t)
{
    if (b->get_oldest_daughter() != NULL)
	for_all_daughters(hdyn, b, d)
	    initialize_system_phase1(d, t);

    if (b->get_kepler() == NULL && b->get_unperturbed_timestep() <= 0) {

	b->set_time(t);
	b->set_system_time(t);
	b->init_pred();

    }
}

//=============================================================================
//  Predictor and corrector steps for the Hermite integration scheme
//=============================================================================

//-----------------------------------------------------------------------------
// predict_loworder_all:  Predict the positions of all nodes at time t.
//			  Now defined in _dyn_/_dyn_ev.C (Steve, 6/99).
//-----------------------------------------------------------------------------

//void predict_loworder_all(hdyn * b,
//			    xreal t)
//{
//    if (!b) return;			// shouldn't happen -- flag as error?
//
//    if (b->is_parent())
//	for_all_daughters(hdyn, b, d)
//	    predict_loworder_all(d, t);
//
//    b->predict_loworder(t);
//}

//-----------------------------------------------------------------------------
// correct_and_update:  Apply the corrector for 3rd-order hermite scheme
//			and update time and timestep.
//
// See Makino & Aarseth (PASJ 1992) for derivation of corrector terms.
//
//-----------------------------------------------------------------------------

local inline void get_derivs(vec& acc, vec& jerk,
			     vec& old_acc, vec& old_jerk,
			     real dt, vec& bt2, vec& at3)
{
    bt2 = -3 * (old_acc - acc) - dt * (2 * old_jerk + jerk);
    at3 =  2 * (old_acc - acc) + dt * (old_jerk + jerk);
}

local inline void update_derivs_from_perturbed(vec& acc_p,
					       vec& jerk_p,
					       vec& old_acc_p,
					       vec& old_jerk_p,
					       real kdt,
					       real ki2,	// not used
					       real ki3,	// not used
					       vec& acc,
					       vec& jerk,
					       vec& bt2,
					       vec& at3,
					       bool print = false)
{
    // Return here to neglect all perturbed terms and compute slow binary
    // CM motion in the center of mass approximation.

    // if (CM_ONLY) return;

    vec bt2_p, at3_p;
    get_derivs(acc_p, jerk_p, old_acc_p, old_jerk_p,
	       kdt, bt2_p, at3_p);

    // Don't scale here -- should be implicit in acc and jerk.

    // bt2_p *= ki2;
    // at3_p *= ki3;

    acc += acc_p;
    jerk += jerk_p;
    bt2 += bt2_p;
    at3 += at3_p;

    if (print) {
	cerr << endl;
	PRC(abs(acc-acc_p)); PRL(abs(acc_p));
	PRC(abs(jerk-jerk_p)); PRL(abs(jerk_p));
	PRC(abs(bt2-bt2_p)); PRL(abs(bt2_p));
	PRC(abs(at3-at3_p)); PRL(abs(at3_p));
    }
}

local inline void correct_slow(slow_binary *s, real dt,
			       vec &acc, vec &jerk,
			       vec &bt2, vec &at3,
			       int id)
{

    // Internal slow motion operates on timescale dtau, not dt.
    // This effectively reduces vel by kappa, acc_p by kappa^2,
    // jerk_p by kappa^3, etc.

    real kap = s->get_kappa();

    real ki = 1/kap;
    real ki2 = ki*ki;
    real ki3 = ki2*ki;

    // Note from Steve, 9/99: the manipulations here would be much
    // simpler if we had direct (friend?) access to the acc_p (etc.)
    // data in the slow_binary class -- can this be done?
    //
    // For now, simply work with copies and copy back changes as
    // necessary.

    vec acc_p = s->get_acc_p();
    vec jerk_p = s->get_jerk_p();

    // Scale acc_p and jerk_p first because the old_ versions are
    // already scaled.

    acc_p *= ki2;
    jerk_p *= ki3;

    // Propogate the changes back to the class.

    s->set_acc_p(acc_p);
    s->set_jerk_p(jerk_p);

#if 0

    // Critical point for consistency is that the computed d(acc)
    // should be approximately kap*dt*jerk before scaling, or
    // dt*jerk afterwards...

    PRL(acc_p-old_acc_p);
    PRL(0.5*dt*(jerk_p+old_jerk_p));

#endif

    vec old_acc_p = s->get_old_acc_p();
    vec old_jerk_p = s->get_old_jerk_p();

    update_derivs_from_perturbed(acc_p, jerk_p, old_acc_p, old_jerk_p,
				 dt, ki2, ki3,
				 acc, jerk, bt2, at3);
}

#define ONE12 0.0833333333333333333333
#define ONE3  0.3333333333333333333333

bool hdyn::correct_and_update()
{
#if 0
    if (is_low_level_node()) {
	cerr << "old_acc  " << old_acc << endl;
	cerr << "    acc  " << acc << endl;
	cerr << "old_jerk " << old_jerk << endl;
	cerr << "    jerk " << jerk << endl;
	cerr << " dt " << timestep << endl;
	cerr << "pos " << pos << endl;
	cerr << "vel " << vel << endl;
	cerr << "pos " << pred_pos << endl;
	cerr << "vel " << pred_vel << endl;
	cerr << "pos " << get_younger_sister()->pos << endl;
	cerr << "vel " << get_younger_sister()->vel << endl;
	cerr << "pos " << get_younger_sister()->pred_pos << endl;
	cerr << "vel " << get_younger_sister()->pred_vel << endl;
    }
#endif

    // As of 4/03, the step we are taking is not necessarily the time step
    // for this particle.  Rather, the step runs from time to system_time,
    // and the actual time step may be less than timestep.  This is
    // necessary to allow synchronization of nodes without explicitly
    // modifying their timestep, so timestep remains a good measure of the
    // natural step, and also a power of 2.

    // Determine the actual step (to be used in calculating derivatives).

    // real dt = timestep;
    real dt = system_time - time;	// not necessarily timestep

    // if (dt != timestep) {
    //	cerr << "correct_and_update: " << format_label() << ": ";
    // 	PRC(dt); PRL(timestep);
    // }

    if (dt <= 0) {

	// Nothing to be done (may happen e.g. if a low-level node
	// forces a parent to be integrated ahead of schedule...).

	return true;
    }

    if (slow) dt = slow->get_dtau();

    // For slow binaries and their perturbers, acc and jerk are computed
    // in the center of mass approximation, with the perturbative terms
    // saved separately for special treatment, but old_acc and old_jerk
    // necessarily include the (adjusted) perturbative terms, as they are
    // used in prediction.  Must adjust them here before proceeding with
    // the correction.  Note that this is the only function where such
    // modification occurs.

    hdyn *od = get_oldest_daughter();
    bool is_top_level = is_top_level_node();

    // For treatment of low-level slow binaries...

    int low_slow_comp = 0;		// indicates which component, if
					// any, is slow (0: neither, 1: od,
					// 2: ys, 3: both)
    hdyn *sb = NULL;

    if (is_top_level) {

	if (od && od->slow) {			// slow binary CM

	    old_acc -= od->slow->get_old_acc_p();
	    old_jerk -= od->slow->get_old_jerk_p();

	} else if (sp) {			// slow binary perturber

	    slow_perturbed *s = sp;
	    while (s) {				// correct for all perturbees
		old_acc -= s->get_old_acc_p();
		old_jerk -= s->get_old_jerk_p();
		s = s->get_next();
	    }
	}

    } else {

	if (od && od->slow) {
	    low_slow_comp = 1;		// this is a low-level slow binary CM
	    sb = od;
	}

	hdyn *od2 = get_younger_sister()->get_oldest_daughter();

	if (od2 && od2->slow) {

	    low_slow_comp += 2;		// sister is a low-level slow binary CM

	    if (!sb || sb->get_kappa() < od2->get_kappa()) sb = od2;
	}


	low_slow_comp = 0;		// *** not implemented, for now... ***


	if (low_slow_comp) {

	    // Low-level slow binaries are resolved into components
	    // for purposes of computing the sister interaction...

	    // Need to remove the perturbative portions of both the old
	    // and current accs and jerks.  Compute the acc and jerk
	    // terms here.  Old info is stored in the slow data structure
	    // of the (elder) slow binary sb.

	    old_acc -= sb->slow->get_old_acc_p();
	    old_jerk -= sb->slow->get_old_jerk_p();

	    // Compute relative acceleration and jerk due to sister
	    // in the point-mass approximation.
	    //
	    // Hmmm.....
	    //
	    // What we need is to compute the 2-body force in the point-mass
	    // approximation, but taking the external perturbation into account.
	    // Then subtracting the non-point-mass force leaves just the
	    // perturbative part of the interaction, for correction below.
	    //
	    // However, it looks as though calculate_partial_acc_and jerk()
	    // computes only the force on this node, not the relative force
	    // between the components; it also seems to resolve the second
	    // component if it is a binary, which is not what we want here
	    // (However, specifying point mass mode will do the right thing
	    // if the first component is a binary...)
	    //
	    // Function calculate_acc_and_jerk_on_low_level_node() uses
	    // calculate_partial_acc_and jerk(), so it isn't what we want
	    // either...
	    //
	    // If we redo this calculation by hand, we end up computing the
	    // perturber interaction twice...


	    vec a_2b, j_2b;
	    real p_2b;
	    hdyn *p_nn_sister;
	    real d_min_sister = VERY_LARGE_NUMBER;


	    calculate_partial_acc_and_jerk(get_parent(), get_parent(), this,
					   a_2b, j_2b, p_2b,
					   d_min_sister, p_nn_sister,
					   USE_POINT_MASS,
					   get_parent()->find_perturber_node(),
					   this);	      // node to charge


	    // Save CM pieces only in acc and jerk; store the rest in acc_p
	    // and jerk_p.

	    PRL(old_acc);
	    PRL(acc);
	    PRL(a_2b);

	    vec dx = pred_pos - get_younger_sister()->pred_pos;
	    real r2 = dx*dx;
	    vec a = -get_younger_sister()->mass * dx / pow(r2, 1.5);

	    PRL(a);

	    // pp3(get_parent());

	    sb->slow->set_acc_p(acc-a_2b);
	    sb->slow->set_jerk_p(jerk-j_2b);

	    acc = a_2b;
	    jerk = j_2b;
	}
    }

    vec bt2, at3;
    get_derivs(acc, jerk, old_acc, old_jerk, dt, bt2, at3);

    // Reminder:	bt2  =  a''  dt^2 / 2		(a' = j)
    //			at3  =  a''' dt^3 / 6

    // For slow systems or perturbers, the derivatives at this point
    // are in the center of mass approximation.  We now include the
    // perturbative terms separately.  There are three circumstances in
    // which a perturbative term must be included: (1) this is the CM
    // of a slow binary, (2) this is (the top_level_node of) a perturber
    // of one or more slow binaries, and (3) this is part of a multiple
    // system containing a slow binary.  In each case, acc and jerk are
    // the now center-of-mass part only; the relevant correction terms
    // are saved in slow->acc_p and slow->jerk_p, or in the corresponding
    // entries in the slow_perturbed list.

    // As currently coded, a slow binary sees any other slow binaries on
    // its perturber list in the center-of-mass approximation only, so the
    // issue of how to scale the slow-slow interaction doesn't arise.

#if 0

    PRL(acc-old_acc);
    PRL(0.5*dt*(jerk+old_jerk));

#endif

    if (is_top_level) {

	if (od && od->slow)				// case (1)

	    correct_slow(od->slow, dt, acc, jerk, bt2, at3, 1);

	else if (sp) {					// case (2)

	    // Loop over all slow binaries known to be perturbed by this node.
	    // Scale the derivatives, but do not apply any correction if kappa
	    // has changed or the old_ quantities are zero -- the default is
	    // thus to use the CM approximation in those cases.  Function
	    // store_old_force will properly set the old_ quantities later.

	    slow_perturbed *s = sp;
	    bool cleanup = false;

	    // cerr << "correct: checking slow_perturbed list of "
	    //      << format_label() << endl;

	    while (s) {

		hdyn *bj = (hdyn*)s->get_node();

		od = NULL;
		if (bj && bj->is_valid())
		    od = bj->get_oldest_daughter();

		if (od && od->get_slow()) {

		    // We have a valid slow binary node.  Only compute
		    // the correction if kappa is what we expect and
		    // old_acc_p is nonzero.

		    bool corr = true;

		    int kap = s->get_kappa();
		    if (kap != od->get_kappa()) {

			kap = od->get_kappa();
			s->set_kappa(kap);

			corr = false;

#if 0
			cerr << "correct: skipping correction for "
			     << bj->format_label()
			     << " because kappa doesn't match: " << endl;
			PRC(kap); PRL(od->get_kappa());
#endif
		    }

		    real ki = 1.0/kap;
		    real ki2 = ki*ki;
		    real ki3 = ki2*ki;

		    vec acc_p = s->get_acc_p();     // as above, better to
		    vec jerk_p = s->get_jerk_p();   // have _dyn_ as friend?

		    // Scale acc_p and jerk_p (see comment above).

		    acc_p *= ki2;
		    jerk_p *= ki3;

		    // Propogate the changes back to the class.

		    s->set_acc_p(acc_p);
		    s->set_jerk_p(jerk_p);

		    if (corr) {

			// Apply the correction if the old_ quantities are set.

			vec old_acc_p = s->get_old_acc_p();

			if (square(old_acc_p) > 0) {

			    vec old_jerk_p = s->get_old_jerk_p();
			    real kdt = dt; // kap*dt;

			    update_derivs_from_perturbed(acc_p, jerk_p,
							 old_acc_p, old_jerk_p,
							 kdt, ki2, ki3,
							 acc, jerk, bt2, at3);

#if 0
			    cerr << "correct: corrected " << format_label();
			    cerr << " for slow binary " << bj->format_label()
				 << endl;
#endif

			} else {

#if 0
			    cerr << "correct: skipping correction for "
				 << bj->format_label()
				 << " because old_acc_p isn't set" << endl;
#endif

			}
		    }

		} else {

		    // This shouldn't happen, but may conceivably occur when
		    // the perturber lists have been disrupted.  Do nothing,
		    // but flag a warning for now.

		    int p = cerr.precision(INT_PRECISION);
		    cerr << "correct: warning: inconsistency in slow_perturbed "
			 << "correction" << endl
			 << "         for node " << format_label()
			 << " at time " << get_time()
			 << endl;
		    cerr.precision(p);

		    PRC(bj); PRL(bj->format_label());
		    PRL(od);
		    if (od) PRL(od->get_slow());

		    cleanup = true;

		}

		s = s->get_next();
	    }

	    if (cleanup)
		check_slow_perturbed(diag->slow_perturbed);

	}

    } else if (low_slow_comp)				// case (3)

	correct_slow(sb->slow, dt, acc, jerk, bt2, at3, 2);


    vec new_pos = pred_pos + (0.05 * at3 + ONE12 * bt2) * dt * dt;
    vec new_vel = pred_vel + (0.25 * at3 + ONE3 * bt2) * dt;


//    if (name_is("11") && system_time > 1.48 && system_time < 1.484377) {
//  	PRL(pred_pos);
//  	PRL(pred_vel);
//  	PRL(at3);
//  	PRL(bt2);
//  	PRL(dt);
//  	PRL(new_pos);
//  	PRL(new_vel);
//    }

    if (has_grape4()) {					// GRAPE-4 only

	// Note from Steve (10/01):
	//
	// Need to check for possible hardware errors (at least on the GRAPE-4).
	// By construction, no quantity should change much from one step to the
	// next.  Thus, a large change in acc or jerk (on in new_vel relative to
	// vel, which combines the two in a natural way) probably indicates a
	// problem.  One difficulty with checking new_vel is that supernovae may
	// result in legal but large kick velocities.  However, these should
	// be applied *after* the dynamical step, in which case both the old and
	// the new velocities should be large.
	//
	// The GRAPE-4 error of 10/01 occurs exclusively in acc and in only one
	// component, usually (but not always) z.  Thus, another possible check
	// might be to see if the velocity, acc, or jerk (whichever is "large")
	// is directed predominantly along one coordinate axis.
	//
	// Probably best to use vel as an indicator, then apply successively
	// finer criteria before declaring a problem.  Once a problem is found,
	// we simply quit this function, returning false.  It is up to the
	// calling function to take corrective action.  Currently, this normally
	// consists of recomputing the acc and jerk on the front end and calling
	// this function again (just one retry, so we should be sure not to flag
	// high-speed neutron stars...).
	//
	// Note also that, if the acc/jerk error is sufficiently small that it
	// doesn't register significantly in vel, then there is no need to take
	// action, as the problem seems to disappear from one GRAPE call to the
	// next.

	real old_v = abs1(vel);

	// Numbers here are somewhat arbitrary (but see above note).

	// real new_v = abs1(new_vel);
	// if (new_v/old_v > 1000 || new_v > 1.e6) {	// too loose...

	real dv = abs1(new_vel-vel);
	if (dv > old_v && dv > 0.5) {

	    // Possible runaway -- speed has changed significantly.

	    // Refine the possibilities before flagging an error.
	    // Neutron star shouldn't show a large delta(vel), and the acc
	    // or jerk should be good indicators of problems in any case.

	    if (abs1(acc-old_acc) > abs1(old_acc)
		|| abs1(jerk-old_jerk) > 5*abs1(old_jerk)) {

	      	if (diag->grape && diag->grape_level > 0) {

		    cerr << endl << "correct: possible hardware error at time "
		         << get_system_time() << endl;

		    if (has_grape4())
		        PRL(get_grape4_chip(this));	// direct access to data
							// in hdyn_grape4.C

#if 0
		    cerr << endl << "pp3 with old pos and vel:" << endl;
		    pp3(this);
#else
		    PRL(old_acc);
		    PRL(old_jerk);
		    PRL(acc);
		    PRL(jerk);
#endif
		}

		// cerr << endl << endl << "System dump:" << endl << endl;
		// pp3(get_root());

		// Options:	quit
		//		restart
		//		flag and continue
		//		discard new_pos and new_vel, retain
		//		    old acc and jerk and continue
		//		flag, recompute acc and jerk, and continue  <--

		// Flag the problem internally.

		char tmp[128];
		sprintf(tmp, "runaway in correct at time %f", (real)time);
		log_comment(tmp);

		int n_runaway = 0;
		if (find_qmatch(get_log_story(), "n_runaway"))
		    n_runaway = getiq(get_log_story(), "n_runaway");

		n_runaway++;
		if (diag->grape && diag->grape_level > 0)
		    PRL(n_runaway);

		if (n_runaway > 10)
		    exit(0);		// pretty liberal, as errors are
					// usually associated with pipes,
					// not stars...

		putiq(get_log_story(), "n_runaway", n_runaway);

		return false;		// trigger a retry on return; don't
					// even bother with update
	    }
	}
    }

#if 0
    if (name_is("3007")) {
      int p = cerr.precision(12);
      cerr << endl; PRC(format_label()); PRL(time);
      PRL(new_pos-pred_pos);
      PRL(new_pos-pos);
      PRL(new_vel-pred_vel);
      PRL(new_vel-vel);
      PRL(2*pos);
      PRL(2*vel);
      PRL(2*pred_pos);
      PRL(2*pred_vel);
      PRL(2*new_pos);
      PRL(2*new_vel);
      PRL(2*old_acc);
      PRL(2*old_jerk);
      PRL(2*acc);
      PRL(2*jerk);
      PRL(2*bt2*2/pow(timestep,2));
      PRL(2*at3*6/pow(timestep,3));
      PRL(-3 * (old_acc - acc));
      PRL(dt * (2 * old_jerk + jerk));
      cerr.precision(p);
    }
#endif

    pos = new_pos;
    vel = new_vel;

    //--------------------------------------------------------
    // Call update here to avoid recomputation of bt2 and at3.
    //--------------------------------------------------------

    update(bt2, at3);		// sets time

#if 0
    if (name_is("3007")) {
      PRL(timestep);
    }
#endif

    return true;		// normal return value
}


//
// Some notes on the force-evaluation functions (Steve, 7/98).
// ----------------------------------------------------------
//
// Functions (note: "force" <--> "acc, jerk, and potential"):
//
// accumulate_acc_and_jerk
// 	compute force due to a specified particle
//
// hdyn::flat_calculate_acc_and_jerk
// 	compute force on top-level node 'this' due to all other
// 	top-level nodes in the system (masking allowed)
//
// hdyn::perturber_acc_and_jerk_on_leaf
// 	compute force on 'this' leaf due to its top-level
// 	perturber list
//
// hdyn::tree_walk_for_partial_acc_and_jerk_on_leaf
// 	accumulate force on 'this' leaf while traversing the portion
// 	of the tree under a specified node, excluding everything
// 	under a mask node.  Currently, the descent STOPS at an
// 	unperturbed center of mass.
//
// hdyn::calculate_partial_acc_and_jerk_on_leaf
// 	compute the force on `this' leaf from all particles under
// 	a specified node, excluding everything under a mask node
//
// hdyn::calculate_partial_acc_and_jerk
// 	compute the force on `this' node or leaf from all particles
// 	under a specified node, excluding everything under a mask node
//
// hdyn::calculate_acc_and_jerk_on_low_level_node
// 	compute the force on `this' low-level node or leaf
//
// hdyn::calculate_acc_and_jerk_on_top_level_node
// 	compute the force on `this' top-level node or leaf
//
// hdyn::calculate_acc_and_jerk
// 	compute the force on `this' node or leaf
//
// correct_acc_and_jerk
// 	correct computed force in situations where GRAPE-style functions
// 	lead to incorrect results
//
//
// To allow for efficient utilization of the GRAPE hardware from kira,
// function calculate_acc_and_jerk_on_top_level_node is divided into
// calls to the following functions:
//
// hdyn::top_level_node_prologue_for_force_calculation
// 	performs entire calculation in exact case; setup otherwise
// hdyn::top_level_node_real_force_calculation
// 	top-level force calculation only
// 	    calls flat_calculate_acc_and_jerk
// 	REPLACED by grape4/6_calculate_acc_and_jerk in kira if GRAPE flag set
// hdyn::top_level_node_epilogue_force_calculation
// 	performs remainder (non-GRAPE) of force calculation
// 	    calls calculate_partial_acc_and_jerk
//

//
// Calling sequences (schematic, "exact" = false)
//
//                               calculate_acc_and_jerk
//
//                                          |
//                                          |
//                     ------------------------------------------
//                    |                                          |
//                    v                                          v
//
// calculate_acc_and_jerk_on_top_level_node   calculate_acc_and_jerk_on_low_level_node
//
//                    |                                          |
//                    v                                          v
//
//     flat_calculate_acc_and_jerk               calculate_partial_acc_and_jerk
//      (GRAPE-style, or GRAPE itself)           (relative motion, perturbation)
//                  +
//     calculate_partial_acc_and_jerk                            |
//       (perturbations, etc.)                                   |
//                                                               |
//                   |                                           |
//                    -------------------------------------------
//                                          |
//                                          |
//                                          v
//
//                           (calculate_partial_acc_and_jerk)
//
//                                          |
//                                          |
//                    -------------------------------------------
//                   |                                           |
//                   v                                           v
//
// calculate_partial_acc_and_jerk_on_leaf        calculate_partial_acc_and_jerk
//                                                 (RECURSIVE, for acc and jerk
//                   |                                         on a node)
//                   |
//                   |-------------------------------------------
//                   |                                           |
//                   v                                           v
//
//     perturber_acc_and_jerk_on_leaf       tree_walk_for_partial_acc_and_jerk_on_leaf
//       (if perturber list used)
//                                                               |
//                   |-------------------------------------------
//                   |                                           |
//                   v                                           v
//
// 	accumulate_acc_and_jerk             tree_walk_for_partial_acc_and_jerk_on_leaf
//                                            (RECURSIVE, if source node is not
//                                                        a leaf)
//


//=============================================================================
//  functions for calculating acceleration & jerk, for given nodes this and b
//=============================================================================

//-----------------------------------------------------------------------------
// accumulate_acc_and_jerk:  Accumulate acc, jerk, and potential
//           from particle b.  Note that the relative position and
//           velocity of b relative to the calling particle are
//           calculated by the caller and passed as arguments.
//           The only reason to pass b here is to access its mass.
//
//	     Note: this function is only referenced from three places:
//
//		flat_calculate_acc_and_jerk()			--> direct
//		perturber_acc_and_jerk_on_leaf()		--> indirect
//		tree_walk_for_partial_acc_and_jerk_on_leaf()	--> indirect
//
//-----------------------------------------------------------------------------

inline void accumulate_acc_and_jerk(hdyn* b,
				    vec& d_pos,
				    vec& d_vel,
				    real eps2,
				    vec& a,
				    vec& j,
				    real& p,
				    real& distance_squared)
{
    distance_squared = d_pos * d_pos;
    real r2inv = 1.0 / (distance_squared + eps2);
    real a3 = -3.0 * (d_pos * d_vel) * r2inv;
    real mrinv = b->get_mass() * sqrt(r2inv);
    p -= mrinv;
    real mr3inv = mrinv * r2inv;
    a += mr3inv * d_pos;
    j += mr3inv * (d_vel + a3 * d_pos);
}


//-----------------------------------------------------------------------------
// flat_calculate_acc_and_jerk: calculate acc and jerk on top node 'this'
// due to all other top-level nodes in the system except 'mask', always in
// the point-mass approximation, optionally determining a perturber list.
//
// (Overloaded function)
//-----------------------------------------------------------------------------

int hdyn::flat_calculate_acc_and_jerk(hdyn * b,    	// root node
				      bool make_perturber_list)
{

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "    flat_calculate_acc_and_jerk for "
	     << format_label() << endl;
    }
#endif

    acc = jerk = 0;
    pot = 0;

    // Initialize both nn and coll:

    nn = NULL;
    d_nn_sq = VERY_LARGE_NUMBER;
    coll = NULL;
    d_coll_sq = VERY_LARGE_NUMBER;

    real distance_squared;

    // ***********************************************************************
    // ** Note: the perturber list code here does not (yet) allow for "low" **
    // ** perturber lists in the event that the top-level list overflows.   **
    // ** See hdyn_grape6.C.  (Steve, 3/03)                                 **
    // ***********************************************************************

    int n_top = 0;

#ifdef USEMPI
    n_top = my_start_count;
#endif

    // Note special care in interpreting radius.  This code is duplicated
    // in function get_sum_of_radii (kira_encounter.C).

    real rad = get_radius();
    bool is_bh = false;
    if (rad < 0) {
	is_bh = true;
	rad = -rad;
    }

    // The use of get_pred_pos() and get_pred_vel() below will trigger
    // calls to predict_loworder() even though no prediction is
    // necessary here, as that was done previously.  Most of the time
    // the function will return without making a prediction, but the 
    // cumulative effect of all these calls within the inner loop can be
    // substantial.  Could replace by member functions get_nopred_xxx(),
    // which simply return the vector without prediction, but simpler just
    // to  access bi->pred_pos and bi->pred_vel directly.   Steve (1/05)


    // *** The following loop is a candidate for threading (Steve, 5/05). ***

#ifdef USEMPI

    // wwvv mpidoit is for determining if the calculation has to be skipped 
    // (another process will take care of this calculation) or not
    // wwvv could have used daughter_counter or n_top for the same
    // purpose, but again, I want to avoid double uses of the same variable

//    int mpidoit;
//    mpidoit = mpi_myrank;

    // wwvv
    // in the MPI case, we calculate the acc and jerk only for a part of 
    // the daughters in each process. Later on, we add things up..

#endif
#ifdef USEMPI
    int daughter_counter=0;
    for ( hdyn* bi = my_start_daughter; 
	  daughter_counter++ < my_daughter_count;
	  bi = bi->get_younger_sister())
#else
    for_all_daughters(hdyn, b, bi) 
#endif
    {
	n_top++;

#ifdef USEMPI

	// Check if we have to calculate for this daughter...

//	if (mpidoit-- != 0) continue;
//	mpidoit = mpi_nprocs-1;
	hdyn_ev_count++;
#endif
	if (bi != this) {

	    // Determine the sum of the radii for use in finding the
	    // coll.  Note that we do NOT check the story for black hole
	    // information -- now the default for get_sum_of_radii().

	    vec d_pos = bi->pred_pos - pred_pos;
	    vec d_vel = bi->pred_vel - pred_vel;
	    accumulate_acc_and_jerk(bi,				// (inlined)
				    d_pos, d_vel,
				    eps2, acc, jerk, pot, distance_squared);

	    // Note from Steve (1/05):  Must be very careful in this loop
	    // to avoid extra function calls or memory accesses that may
	    // cause cache misses.  Best to confine calculations to data
	    // close to the basic dynamical quantities, and to inline any
	    // functions used.  Specifically, the use of black hole tidal
	    // radii below can cause significant performance degradation
	    // if not handled carefully.  In the absence of a GRAPE, this
	    // is the innermost loop that dominates the total cost of the
	    // integration.

	    // Determine the sum of the radii for use in finding coll.
	    // Don't look in the starbase or the log story for flags.

	    real radi = bi->radius;
	    bool bi_is_bh = false;
	    if (radi < 0) {
		bi_is_bh = true;
		radi = -radi;
	    }

	    real sum_of_radii
		= compute_sum_of_radii(this, rad, is_bh,
				       bi, radi, bi_is_bh);	// (inlined)

	    // The nn and coll passed to update_nn_coll here are the actual
	    // pointers, so this update changes the data in b.  

	    update_nn_coll(this, 1,		// (1 = ID)	// (inlined)
			   distance_squared, bi, d_nn_sq, nn,
			   sum_of_radii, d_coll_sq, coll);

	    // Note: we do *not* take the slowdown factor into account
	    //	     when computing the perturber list.

	    // See equivalent code for use with GRAPE in
	    // hdyn_grape[46].C/get_neighbors_and_adjust_h2.

	    if (make_perturber_list
		&& is_perturber(this, bi->mass,			// (inlined)
				distance_squared,
				perturbation_radius_factor)) {
		if (n_perturbers < MAX_PERTURBERS) {

#ifdef USEMPI
#ifdef MPIDEBUG
		    cerr << ":"<<mpi_myrank<<": TD:"<<__FILE__<<":"
			 <<__LINE__<<" perturber modified"<<endl;
#endif
#endif

		    perturber_list[n_perturbers] = bi;
#ifdef USEMPI
		    // to keep track of the order in which the
		    // perturbers would have been found in the serial
		    // version:

		    seq[n_perturbers] = n_top;
#endif
#if 0
		    cerr << "added " << bi->format_label();
		    cerr << " to perturber list of "
			 << format_label()
			 << endl;
#endif
		}
		n_perturbers++; // note: n_pertubers can get bigger than
		                // MAX_PERTURBERS, but only MAX_PERTURBERS
				// are saved. wwvv
	    }
	}
    }

#ifdef USEMPI
    n_top = num_daughters;
#endif
    return n_top;
}

void hdyn::perturber_acc_and_jerk_on_leaf(vec &a,
					  vec &j,
					  real &p,
					  real &p_d_nn_sq,
					  hdyn *&p_nn,
					  hdyn *pnode,
					  hdyn *step_node)
{
    // Calculate acc and jerk on leaf 'this' due to the perturber
    // list associated with node pnode.

    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "        perturber_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
    }
    dbg_message("perturber_acc_and_jerk_on_leaf", this);
#endif

    a = j = 0.0;
    p = 0;

    // Do *not* set p_d_nn_sq = VERY_LARGE_NUMBER here!
    //						(Steve, 11/00)

    d_coll_sq = VERY_LARGE_NUMBER;

    // Note from Steve (3/03): pnode may have "low-only" perturbers available.
    // Use this list if necessary (better than an O(N) calculation).

    if (!pnode) return;

    int np = 0;
    if (pnode->valid_perturbers)
	np = pnode->n_perturbers;
    else if (pnode->valid_perturbers_low)
	np = pnode->n_perturbers_low;

    if (np <= 0) return;
    
    // Determine absolute position and velocity of 'this', using explicit code.
    // Removed four "hdyn_something_relative_to_root" references, which are
    // quite inefficient (Steve, 8/20/98).

    // vec d_pos = -hdyn_something_relative_to_root(this,
    // 						    &hdyn::get_pred_pos);
    // vec d_vel = -hdyn_something_relative_to_root(this,
    // 						    &hdyn::get_pred_vel);

    vec d_pos = -pred_pos;
    vec d_vel = -pred_vel;

    hdyn* par = get_parent();
    if (par) {
 	hdyn* gpar = par->get_parent();
 	while (gpar) {
 	    d_pos -= par->pred_pos;
 	    d_vel -= par->pred_vel;
 	    par = gpar;
 	    gpar = gpar->get_parent();
 	}
    }

    // Loop over the perturber list.

    vec d_pos_p;
    vec d_vel_p;

    bool reset = false;

    for (int i = 0; i < np; i++) {

	hdyn *bb = pnode->perturber_list[i];

	if (!bb || !bb->is_valid()) {

	    // Should not occur, and the logic in the calling function
	    // calculate_partial_acc_and_jerk_on_leaf assumes that it
	    // does not occur (so we can't just quit here).  Use any
	    // valid perturbers on the list, but flag a problem and
	    // force recomputation of the list.

	    if (!reset) {

		cerr << endl
		     << "warning: perturber_acc_and_jerk_on_leaf:"
		     << endl
		     << "         NULL or invalid perturber "
		     << bb << " for node "
		     << format_label() << " at time " << time
		     << endl
		     << "         forcing recomputation of perturber list"
		     << endl;

		pnode->valid_perturbers = false;
		reset = true;
	    }

	    // Go on to the next perturber on the list.

	    continue;
	}

	// Seems to be necessary the first time through -- perturbers are
	// apparently not (necessarily) predicted elsewhere (Steve 8/13/98):

	if (!elder_sister)
	    predict_loworder_all(bb, time);

	// Determine absolute position and velocity of perturber, using
	// explicit code rather than hdyn_something_relative_to_root.

// 	if (bb->is_top_level_node()) {		// (time saver)
// 	    d_pos_p = d_pos + bb->pred_pos;
// 	    d_vel_p = d_vel + bb->pred_vel;
//  	} else {
//  	    d_pos_p = d_pos +
//  		hdyn_something_relative_to_root(bb, &hdyn::get_pred_pos);
//  	    d_vel_p = d_vel +
//  		hdyn_something_relative_to_root(bb, &hdyn::get_pred_vel);
//  	}

 	d_pos_p = d_pos + bb->pred_pos;
 	d_vel_p = d_vel + bb->pred_vel;

 	hdyn* par = bb->get_parent();
 	if (par) {
 	    hdyn* gpar = par->get_parent();
 	    while (gpar) {
 		d_pos_p += par->pred_pos;
		d_vel_p += par->pred_vel;
 		par = gpar;
 		gpar = gpar->get_parent();
 	    }
 	}

	real d2_bb;
	
	accumulate_acc_and_jerk(bb, d_pos_p, d_vel_p,		// (inlined)
				eps2, a, j, p, d2_bb);
	
	// Note: p_d_nn_sq and p_nn come from the argument list,
	// and may be copies of d_nn_sq and nn, or the real thing.
	// However, d_coll_sq and coll are real, so this update
	// actually changes the coll data.

	real sum_of_radii = get_sum_of_radii(this, bb);

	update_nn_coll(this, 2,					// (inlined)
		       d2_bb, bb, p_d_nn_sq, p_nn,
		       sum_of_radii, d_coll_sq, coll);

    }

    step_node->inc_indirect_force(np);				// bookkeeping
}


//-----------------------------------------------------------------------------
// tree_walk_for_partial_acc_and_jerk_on_leaf:  Accumulate acc and jerk
//      while traversing the portion of the tree under source node b,
//	excluding everything under node mask (inclusive).  Vectors offset_pos
//	and offset_vel MUST be the position and velocity of b relative to
//	the field point (i.e. the point where acc and jerk are wanted).
//
// 	On return, a, j, and p are the acc, jerk, and pot at the field
//	point, updated to include the effect of b (modulo the value of
//	point_mass_flag), p_d_nn_sq is the minimum distance squared
//	from any component of b to the field point, and p_nn is an
//	updated pointer to the the nearest neighbor of the field point.
//
//	In the present code, the field point is the location of 'this'.
//	It is ASSUMED that 'mask' is 'this' or an ancestor of 'this'.
//-----------------------------------------------------------------------------

// This routine calculates the force by simple recursion.  If a leaf is
// reached, calculate the force.  Otherwise descend the tree.
//
// NOTE: RECURSIVE LOOPS are quite inefficient on some (most?) systems.
// However, by construction, this function is only used in circumstances
// where a system-wide O(N) traversal of the system is NOT needed, so the
// inefficiency is tolerable.  All functions that require O(N) operations
// are coded as explicit loops (see calculate_acc_and_jerk_on_top_level_node).

void hdyn::tree_walk_for_partial_acc_and_jerk_on_leaf(hdyn *b,
						      hdyn *mask,
						      vec& offset_pos,
						      vec& offset_vel,
						      vec& a,
						      vec& j,
						      real& p,
						      real& p_d_nn_sq,
						      hdyn * &p_nn,
						      bool point_mass_flag,
						      hdyn *step_node)
{
    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "        tree_walk_for_partial_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
	cerr << "        b = " << b->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }
    dbg_message("tree_walk_for_partial_acc_and_jerk_on_leaf", this, b);
#endif

    if (b == mask)
	return;

    // Stop descent of the tree at an unperturbed CM, *except* when
    // that CM is the parent of 'this' node (to get the correct 2-body
    // internal force).

    if (b->is_leaf()
	|| (b->is_top_level_node() && point_mass_flag)
	|| (b->kep != NULL && b != parent)
	) {

	if (b != this) {
	    real d2;

	    accumulate_acc_and_jerk(b, offset_pos, offset_vel,	  // (inlined)
				    eps2, a, j, p, d2);

	    step_node->inc_indirect_force();

	    // Note: p_d_nn_sq and p_nn come from the argument list,
	    // and may be copies of d_nn_sq and nn, or the real thing.
	    // However, d_coll_sq and coll are real, so this update
	    // actually changes the coll data.

	    real sum_of_radii = get_sum_of_radii(this, b);

	    update_nn_coll(this, 3,
			   d2, b, p_d_nn_sq, p_nn,
			   sum_of_radii, d_coll_sq, coll);

	} else {

	    // Probably can't occur...

	    cerr << "In tree_walk_for_partial_acc_and_jerk_on_leaf:\n";
	    cerr << "b = "; b->pretty_print_node(cerr);
	    cerr << "  mask = "; mask->pretty_print_node(cerr);

	}

    } else {

	for_all_daughters(hdyn, b, d) {		       // recursive loop

	    vec ppos = offset_pos + d->get_pred_pos();
	    vec pvel = offset_vel + d->get_pred_vel();

	    tree_walk_for_partial_acc_and_jerk_on_leaf(d, mask,
						       ppos, pvel,
						       a, j, p,
						       p_d_nn_sq, p_nn,
						       point_mass_flag,
						       step_node);
	}
    }
}

//=============================================================================
//  functions for determining which nodes should exert forces on each other
//=============================================================================

//-----------------------------------------------------------------------------
// calculate_partial_acc_and_jerk_on_leaf:  Calculate the force on `this' from
//      all particles under `top', excluding everything under node mask.
//      This routine first calculates the position of `top' relative to `this',
//      then calls tree_walk_for_partial_acc_and_jerk_on_leaf.
//-----------------------------------------------------------------------------

void hdyn::calculate_partial_acc_and_jerk_on_leaf(hdyn * top,
						  hdyn * common,
						  hdyn * mask,
						  vec& a,
						  vec& j,
						  real & p,
						  real & p_d_nn_sq,
						  hdyn * &p_nn,
						  bool point_mass_flag,
						  hdyn * pnode,
						  hdyn * step_node)
{
    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "      calculate_partial_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
	cerr << "      top = " << top->format_label();
	cerr << ", common = " << common->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }
    dbg_message("calculate_partial_acc_and_jerk_on_leaf", this);
    dbg_message("                 from particle", top);
#endif

    a = j = 0.0;
    p = 0;

    // Loop over the appropriate perturber list for external perturbations,
    // and traverse the tree below top, masked by mask, for internal
    // perturbations due to other members of the parent clump.

    if (pnode)
	perturber_acc_and_jerk_on_leaf(a, j, p,
				       p_d_nn_sq, p_nn,
				       pnode, step_node);

    if (top == mask) return;

    // Determine absolute position and velocity of "this":

    vec d_pos = 0;
    vec d_vel = 0;

    hdyn *b;
    for (b = this; b != common; b = b->get_parent()) {
	d_pos -= b->get_pred_pos();
	d_vel -= b->get_pred_vel();
    }

    for (b = top; b != common; b = b->get_parent()) {
	d_pos += b->get_pred_pos();
	d_vel += b->get_pred_vel();
    }

    tree_walk_for_partial_acc_and_jerk_on_leaf(top, mask, d_pos, d_vel,
					       a, j, p,
					       p_d_nn_sq, p_nn,
					       point_mass_flag,
					       step_node);
}

//-----------------------------------------------------------------------------
// calculate_partial_acc_and_jerk:  Calculate the force on `this' from
//	all particles under `top', excluding everything under node mask.
//      This routine recursively calls itself and calculates the weighted
//      average of all forces on its daughters.  For a leaf, it calls
//      calculate_partial_acc_and_jerk_on_leaf.  The recursion is OK in
//      this case because it is used only to descend within nodes, not to
//      traverse the entire system.
//-----------------------------------------------------------------------------

void hdyn::calculate_partial_acc_and_jerk(hdyn * top,
					  hdyn * common,
					  hdyn * mask,
					  vec& a,
					  vec& j,
					  real & p,
					  real & p_d_nn_sq,
					  hdyn * &p_nn,
					  bool point_mass_flag,
					  hdyn* pnode,
					  hdyn* step_node)
{
    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.  If this is a node, the nn returned
    // is the closer of its component nns.

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "    calculate_partial_acc_and_jerk for "
	     << format_label() << endl;
	cerr << "    top = " << top->format_label();
	cerr << ", common = " << common->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }
    dbg_message("calculate_partial_acc_and_jerk", this);
#endif

    // Operations of this function:
    //
    // (1) Calculate acc and jerk on a top-level node due to all
    //     other top-level nodes, in the point-mass approximation
    //     (this is what the GRAPE hardware does...).
    // (2) Calculate masked acc and jerk on a node due to the portion
    //     of the tree below top.
    // (3) Calculate acc and jerk on "this" node due to the particles
    //     in its perturber list, recursively looping over daughters.
    // (4) Calculate acc and jerk on "this" node due to the particles
    //     in its perturber list, in the point-mass approximation.
    // (5) Determine the perturber list of a binary during its
    //     center-of-mass step.

    a = j = 0;
    p = 0;
    real mtot = 0;

    // Don't resolve an unperturbed binary (Steve, 8/20/98).

    if (is_leaf() || point_mass_flag || get_oldest_daughter()->kep) {

	calculate_partial_acc_and_jerk_on_leaf(top, common, mask,
					       a, j, p,
					       p_d_nn_sq, p_nn,
					       point_mass_flag,
					       pnode, step_node);
    } else {

	vec a_daughter, a_save;
	vec j_daughter;
	real p_daughter;
	real m_daughter;

	for_all_daughters(hdyn, this, bd) {

	    // Note from Steve, 7/24/98.  Because of the way nns and
	    // colls are passed around, if we compute the acc and jerk
	    // on a leaf in order to obtain the acc and jerk on its parent
	    // node, the nn information of the leaf will be preserved
	    // because p_nn here is actually the node nn, but the coll
	    // information will be overwritten.  If the coll of bd is the
	    // other binary component (as is often the case), an error
	    // will be made.

	    // The best solution would be to treat colls in exactly the
	    // same way as nn pointers.  For now, however, we just save
	    // and restore the coll data of any node referenced here.

	    hdyn* save_coll = bd->get_coll();		// save the
	    real save_d_coll_sq = bd->get_d_coll_sq();	// coll data

	    m_daughter = bd->get_mass();
	    bd->calculate_partial_acc_and_jerk(top, common, mask,
					       a_daughter,
					       j_daughter, p_daughter,
					       p_d_nn_sq, p_nn,
					       point_mass_flag,
					       pnode, step_node);
	    a += m_daughter * a_daughter;
	    j += m_daughter * j_daughter;
	    p += m_daughter * p_daughter;
	    mtot += m_daughter;

	    if (bd->get_elder_sister() == NULL)		// in case we want
		a_save = a_daughter;			// the perturbation

	    bd->set_coll(save_coll);			// restore the
	    bd->set_d_coll_sq(save_d_coll_sq);		// coll data

	}
	real minv = 1 / mtot;
	a *= minv;
	j *= minv;
	p *= minv;

	// Special case:  It is helpful to know the perturbation on a
	// top-level node whose perturber list has overflowed (or which
	// is being calculated exactly).  In that case, the following
	// should all be true:
	//
	//	top = common = root
	//	mask = this
	//	this is top-level node
	//	point_mass_flag = false
	//	pnode = NULL
	//
	//	valid_perturbers = false
	//	n_perturbers <= 0
	//

	// In this case, we should set perturbation_squared before
	// leaving the function.  Note that perturbation_squared is
	// attached to the daughter nodes, but the perturber list,
	// valid_perturbers, and n_perturbers are properties of the
	// parent.

	if (!point_mass_flag
	    && !pnode
	    && top->is_root()
	    && common == top
	    && mask == this
	    && !valid_perturbers) {

	    // (We probably don't need to check all these!)

	    // Need to set a perturbation in this case.  In the stored
	    // accelerations, "save" refers to the older daughter,
	    // "daughter" to the younger daughter.

	    hdyn* od = get_oldest_daughter();
	    hdyn* yd = get_oldest_daughter()->get_younger_sister();

	    // Perturbation_squared does *not* contain the slowdown factor.

	    real r2 = square(od->pos - yd->pos);
	    od->perturbation_squared =
			square(a_save - a_daughter)
			        * square(r2/mass);
 	}
    }
}


// check_add_perturber:  check if p is a perturber of 'this' and add it
//			 to the perturber list if necessary.

void hdyn::check_add_perturber(hdyn* p, vec& this_pos)
{
    // Changed pred to nopred... (Steve, 4/05)

    vec ppos = hdyn_something_relative_to_root(p, &hdyn::get_nopred_pos);

     bool print = false;

     // print = get_top_level_node()->n_leaves() >= 3;

     if (print) {
 	cerr << "check_add_perturber:  checking " << p->format_label()
 	     << " at time " << system_time;
 	cerr << " for node " << format_label() << endl;
 	PRL(mass);
 	PRL(square(this_pos - ppos));
 	PRL(perturbation_radius_factor);
    }

    // Note 1: we do *not* take the slowdown factor into account
    //	       when computing the perturber list.

    // *********************************************************************
    // * Note 2: the perturber list code here does not yet allow for "low" *
    // * perturber lists in the event that the top-level list overflows.   *
    // * See hdyn_grape6.C.  (Steve, 3/03)                                 *
    // *********************************************************************

    if (is_perturber(this, p->mass,
		     square(this_pos - ppos),
		     perturbation_radius_factor)) {	  // (inlined)

	// If p is a compound system, include all components, without
	// checking if they satisfy the perturbation criteria individually.

	if (RESOLVE_UNPERTURBED_PERTURBERS) {
	    for_all_leaves(hdyn, p, pp) {
		if (n_perturbers < MAX_PERTURBERS)
		    perturber_list[n_perturbers] = pp;
		n_perturbers++;
	    }
	} else {
	    for_all_leaves_or_unperturbed(hdyn, p, pp) {
		if (n_perturbers < MAX_PERTURBERS)
		    perturber_list[n_perturbers] = pp;
		n_perturbers++;
	    }
	}
//	if (print) cerr << "...accepted" << endl;
    } // else
//	if (print) cerr << "...rejected" << endl;
}

void hdyn::create_low_level_perturber_list(hdyn* pnode)
{

    // Create a perturber list attached to this node, based on the list found
    // in pnode.  Extra possibility (Steve, 3/03): pnode may not have a complete
    // perturber list, if the "low" option is set, but the list there is still
    // usable for construction of this low-level list.

//   if (system_time >= 44.15328 && system_time <= 44.1533) {
//     cerr << "create_low_level_perturber_list for " << format_label()
// 	 << " at time " << system_time << endl;
//     if (pnode) {
//       PRC(pnode); PRC(pnode->format_label());
//       PRL(pnode->valid_perturbers);
//     } else
//       cerr << "pnode = NULL" << endl;
//   }

    valid_perturbers = false;
    if (!pnode) return;

    int np = 0;
    if (pnode->valid_perturbers)
	np = pnode->n_perturbers;
    else if (pnode->valid_perturbers_low)
	np = pnode->n_perturbers_low;

    if (np < 0) return;

    // Looks as if we have enough info to build the list.

    new_perturber_list();

    perturbation_radius_factor
	= define_perturbation_radius_factor(this, gamma23);

    // Changed pred to nopred... (Steve 4/05).

    vec this_pos = hdyn_something_relative_to_root(this,
						   &hdyn::get_nopred_pos);
    // First add sisters, aunts, etc. up to pnode...

    hdyn* p = this;
    while (p != pnode) {
	check_add_perturber(p->get_binary_sister(), this_pos);
	p = p->get_parent();
    }

    // ...then accept any node on the pnode list that satisfies the
    // inner node perturbation criterion.

    for (int i = 0; i < np; i++)
	check_add_perturber(pnode->perturber_list[i], this_pos);

    valid_perturbers = true;

    if (n_perturbers > MAX_PERTURBERS)			// can't happen?

	remove_perturber_list();

    else {

// 	if (get_top_level_node()->n_leaves() >= 4) {
// 	    cerr << ">>>> this->"; PRL(format_label());
// 	    cerr << ">>>> "; PRL(pnode->format_label());
// 	    cerr << ">>>> "; PRL(pnode->n_perturbers);
// 	    cerr << ">>>> this->"; PRL(n_perturbers);
// 	    print_perturber_list();
// 	}

    }
}

void hdyn::create_low_level_perturber_lists(bool only_if_null) // default = true
{
  // Possible extra option might be to allow each node to construct its list
  // from a specified pnode (e.g. the top-level), in case we don't think the
  // parent list is good enough...

    if (ALLOW_LOW_LEVEL_PERTURBERS && is_parent()) {

        // Create low-level perturber lists for all nodes below this node.
        // Create lists for daughters first, then recursively apply this
        // function to the daughter nodes.

	if (valid_perturbers || valid_perturbers_low) {

	    for_all_daughters(hdyn, this, bb)
		if (bb->is_parent()) {
		    if (!only_if_null
			|| !(bb->valid_perturbers || bb->valid_perturbers_low))
			bb->create_low_level_perturber_list(this);
		    bb->create_low_level_perturber_lists(only_if_null);
		}

	} else {

	    // End the recursion and clear the perturber lists of this
	    // node and all descendants.

	    for_all_nodes(hdyn, this, bb)
		if (bb->is_parent())
		    bb->remove_perturber_list();

	}
    }
}

//-----------------------------------------------------------------------------
// calculate_acc_and_jerk_on_low_level_node:  Calculate the acc (etc) on
//         one component of a binary node in the following steps:
//
// First, calculate the perturbation on the relative motion as the
// the difference between the accelerations of the node and that of
// its sister, multiplied by the fractional mass of the sister.
// These accelerations do not include the mutual interaction between
// the node and its sister.
//
// Second, calculate the acceleration due to the sister.
//
// Third, calculate the total acceleration, given by the sum of the
// perturbation acc and the acc from the sister.
//
// Note that all three components are calculated by one call to
// calculate_partial_acc_and_jerk.
//-----------------------------------------------------------------------------

void hdyn::calculate_acc_and_jerk_on_low_level_node()
{
#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  calculate_acc_and_jerk_on_low_level_node for "
	     << format_label() << endl;
    }
    dbg_message("calculate_acc_and_jerk_on_low_level_node", this);

    if (parent->get_oldest_daughter()->get_younger_sister()
				     ->get_younger_sister() != NULL)
	err_exit("calculate_acc_and_jerk_on_low_level_node: Not binary node");

    hdyn *sister = get_binary_sister();		// too elaborate!
#else
    hdyn *sister = (hdyn*)younger_sister;	// assume binary tree
    if (!sister) {
	sister = (hdyn*)elder_sister;
	if (!sister)
	    return;
//	    err_exit("calculate_acc_and_jerk_on_low_level_node: no sister!");
    }
#endif

    hdyn *root = get_root();

    // New formulation of perturber lists introduced by Steve 8/98:

    // Perturber list may not necessarily be associated with the
    // top-level node.  Any CM node above this node may have a
    // valid perturber list.  Use the lowest-level one, but make
    // sure that sisters and aunts are properly included.

    // Note from Steve (3/03): Don't worry about "low-only" perturber
    // lists in the perturbation calculation.  See note in function
    // perturber_acc_and_jerk_on_leaf() above.

    // Now pass a pointer to the node containing the perturber list
    // (NULL ==> no valid list) instead of bool "use_perturber_list"
    // flag.  Also, since the low-level perturber list includes all
    // members of the parent "clump" except those below the parent
    // node, we need only descend the tree below pnode to pick up the
    // remaining perturbations.

    hdyn* top_level = get_top_level_node();

    hdyn* pnode = find_perturber_node();    // returns first valid node found

    if (ALLOW_LOW_LEVEL_PERTURBERS && pnode && pnode != get_parent()) {

	// Must be a gap somewhere in the pnode chain.  New formulation
	// expects that the chain is intact, so fix that now.  Start by
        // finding the highest-lying gap in the chain (maybe at the top).

	hdyn *tmp = pnode, *start = pnode;

	while (tmp != root) {
	    if (!tmp->valid_perturbers) start = tmp;	// topmost gap
	    tmp = tmp->get_parent();
	}

	if (start != top_level) start = start->get_parent();

	// Note that start can be any node, including the top-level node,
	// and start may or may not itself have valid perturbers (if the
	// top-level node is invalid).

	start->create_low_level_perturber_lists(false);	// all
	pnode = find_perturber_node();			// may now be NULL
    }

//     if (system_time >= 44.15328 && system_time <= 44.1533) {
//       PRL(pnode);
//       for (hdyn *bb = this; bb != root; bb = bb->get_parent()) {
// 	PRC(bb->format_label()); PRL(bb->valid_perturbers);
//       }
//     }

    hdyn *top = pnode;
    int np = -1;

    // From Steve (3/03): If no valid perturber list is found, setting
    // top = root here will force a O(N) calculation of the force using
    // the front end.  We'd prefer not to do this...  In case of overflow,
    // the GRAPE version of the code now creates a "low-only" perturber
    // list, specifically for the computation of low-level perturbations.
    // If the low-level list hasn't yet been constructed using that list,
    // use the top-level list here.  Perturber_acc_and_jerk() will do the
    // right thing if sent a pnode with only low-only data available.

    // A "valid" pnode is one containing either a full perturber list or
    // a valid low-level list.  The two are distinguished by the value of
    // n_perturbers, which is MAX_PERTURBERS + 1 in the latter case.

    bool partial_only = false;

    if (!pnode) {

	// No pnode is available.  Use the top-level node if it at least
	// contains a valid low-level perturber list.

	if (top_level->valid_perturbers_low) {

	    // Use the partial top-level list.

	    top = pnode = top_level;
	    np = pnode->n_perturbers_low;
	    partial_only = true;

	} else {

	    // No list of any sort is available.  At present, this will
	    // force an O(N) calculation on the front end.  Better to
	    // use GRAPE if possible...

	    top = root;
	}

    } else {

	// Use pnode, but remember that it may only contain usable
        // low-level data.

	np = pnode->n_perturbers;
    }

    // The value of np is only used as a boolean flag below, but it
    // actually has the following meaning.  If np < 0 (-1), no perturber
    // list is available and we will have to use the entire system in
    // computing the perturbation.  If np = 0, the list is valid but
    // empty.  If np > 0, the list is valid, although it may only be
    // partial (partial_only = true).  If np > MAX_PERTURBERS, the list
    // is partial.

    if (np > MAX_PERTURBERS) partial_only = true;

    bool grape_pert = false;
    bool full_pert = false;

    // Handle (and flag) problematical cases when they occur.

    if (top == root				// O(N) calculation
	|| partial_only				// incomplete list
	|| np > MAX_PERTURBERS - 10) {		// long list

#if 0
	cerr << endl << "***** perturbers: ";
	if (top == root) cerr << "top = root, ";
        PRC(partial_only); PRL(np);
	if (pnode) PRL(pnode->format_label());
#endif

	// In these cases, we should compute the perturbation due
	// to the entire system, and use the GRAPE if possible.
	// Using top = root will do what is needed on the front-end
	// if no GRAPE is available.  However, if we have a GRAPE,
	// then we should compute the perturbation due to the other
	// components of the clump (top = top_level), then use GRAPE
	// for the rest of the system (new function, Steve 3/05).

	if (top == root) full_pert = true;
	if (has_grape6()) {
	    top = top_level;
	    pnode = NULL;
	    full_pert = true;
	    grape_pert = true;
	}
    }

    if (!kep) {

	// Bookkeeping (if this isn't part of an unperturbed step):

	if (pnode && !full_pert)
	    kc->pert_with_list += pnode->n_perturbers;
	else {
	    kc->pert_without_list++;

//  	    cerr << "no perturber node for " << format_label()
//  		 << " at time " << system_time << endl;
//  	    pp3(get_top_level_node());

	}
    }
    nn = NULL;
    d_nn_sq = VERY_LARGE_NUMBER;
    sister->set_nn(NULL);
    sister->set_d_nn_sq(VERY_LARGE_NUMBER);

    vec apert1, jpert1;
    vec apert2, jpert2;
    real p_dummy;

    if (grape_pert) {

	// Calculate perturbations to this and sister simultaneously
	// on GRAPE.  Don't update nn, etc and use point-mass
	// approximation for all forces.  (Assume for now that the
	// sister or a component will be the nn/coll.)  Probably
	// should expand this function to get nn and coll properly...
	// Note that this call violates the convention that all calls
	// that might use GRAPE be switched in kira_ev.C.
       
	grape6_calculate_perturbation(get_parent(),
				      apert1, apert2, jpert1, jpert2);


	if (get_parent() != top_level) {

	    int p = cerr.precision(HIGH_PRECISION);

	    //cerr << endl
	    //	 <<  "computing perturbation on "
	    //	 << get_parent()->format_label()
	    //	 << " using GRAPE; np = " << np << endl << flush;
	    // cerr << "this = " << format_label() << ", ";
	    // cerr << "sister = " << sister->format_label() << ", ";
	    // PRL(system_time);

	    // Shortened version:

	    cerr << endl
	    	 <<  "GRAPE perturbation on "
	    	 << get_parent()->format_label()
	    	 << " at t = " << system_time << endl << flush;

	    cerr.precision(p);

#if 0

	    // Check: do the O(N) calculation on the front-end too...

	    vec apert1_FE, jpert1_FE;
	    vec apert2_FE, jpert2_FE;
	    real d_nn_sq_dummy;
	    hdyn *nn_dummy;

	    calculate_partial_acc_and_jerk(root, root, top_level,
					   apert1_FE, jpert1_FE, p_dummy,
					   d_nn_sq_dummy, nn_dummy,
					   USE_POINT_MASS,    // explicit loop
					   NULL,	      // no list
					   this);	      // node to charge

	    sister
	      ->calculate_partial_acc_and_jerk(root, root, top_level,
					       apert2_FE, jpert2_FE, p_dummy,
					       d_nn_sq_dummy, nn_dummy,
					       USE_POINT_MASS, // explicit loop
					       NULL,	       // no list
					       this);	       // node to charge

	    cerr << "test 1:" << endl;
	    vec da1 = apert1_FE - apert1;
	    PRC(apert1); PRL(da1/abs(apert1));
	    vec da2 = apert2_FE - apert2;
	    PRC(apert2); PRL(da2/abs(apert2));
	    PRL(apert2 - apert1);
	    PRL(apert2_FE - apert1_FE);

#endif

	}

	// The above calculation *excludes* the perturbation due to
	// other components of the top-level node.  Include them here,
	// if necessary.

	// (Note also that we don't resolve neighboring binaries.)

	if (get_parent() != top_level) {

	    vec apert1_local, jpert1_local;
	    vec apert2_local, jpert2_local;
	    real d_nn_sq_dummy;
	    hdyn *nn_dummy;

	    // (The following calls just repeat earlier code.)

	    calculate_partial_acc_and_jerk(top_level, top_level, get_parent(),
					   apert1_local, jpert1_local, p_dummy,
					   d_nn_sq_dummy, nn_dummy,
					   !USE_POINT_MASS,
					   NULL,
					   this);
	    apert1 += apert1_local;
	    jpert1 += jpert1_local;

	    sister
		->calculate_partial_acc_and_jerk(top_level, top_level,
						     get_parent(),
						 apert2_local, jpert2_local,
						     p_dummy,
						 d_nn_sq_dummy, nn_dummy,
						 !USE_POINT_MASS,
						 NULL,
						 this);
	    apert2 += apert2_local;
	    jpert2 += jpert2_local;

#if 0
	    cerr << "adding perturbation due to aunt and cousins"
		 << endl << flush;

	    PRL(apert1);
	    PRL(apert2);

#endif

#if 0

	    // Alternate test: replacing top_level in the previous
	    // check by get_parent() should also include the
	    // perturbation due to the clump.  Ought to work, but
	    // seems not to...  However, the values just computed seem
	    // to be the right ones!

	    vec apert1_test, jpert1_test;
	    vec apert2_test, jpert2_test;

	    calculate_partial_acc_and_jerk(root, root, get_parent(),
					   apert1_test, jpert1_test,
					       p_dummy,
					   d_nn_sq_dummy, nn_dummy,
					   USE_POINT_MASS,    // explicit loop
					   NULL,	      // no list
					   this);	      // node to charge

	    sister
	      ->calculate_partial_acc_and_jerk(root, root, get_parent(),
					       apert2_test, jpert2_test,
					           p_dummy,
					       d_nn_sq_dummy, nn_dummy,
					       USE_POINT_MASS,// explicit loop
					       NULL,	      // no list
					       this);	      // node to charge

	    cerr << "test 2:" << endl;
	    vec da1 = apert1_test - apert1;
	    PRL(apert1); PRC(apert1_test); PRL(da1/abs(apert1));
	    vec da2 = apert2_test - apert2;
	    PRL(apert2); PRC(apert2_test); PRL(da2/abs(apert2));
	    PRL(apert2 - apert1);
	    PRL(apert2_test - apert1_test);




	    hdyn *aunt = get_parent()->get_binary_sister();

	    if (aunt) {
	      cerr << "relative to parent: " << endl;
	      PRL(aunt->format_label());
	      vec thispos = pred_pos;
	      vec sispos = sister->pred_pos;
	      vec auntpos = aunt->get_nopred_pos()-get_parent()->pred_pos;
	      PRL(thispos);
	      PRL(sispos);
	      PRL(auntpos);
	      PRL(thispos-auntpos);
	      PRL(abs(thispos-auntpos));
	      vec tmp = thispos-auntpos;
	      vec aprt1 = -aunt->mass*tmp/pow(square(tmp), 1.5);
	      PRL(aprt1);
	      PRL(sispos-auntpos);
	      PRL(abs(sispos-auntpos));
	      tmp = sispos-auntpos;
	      vec aprt2 = -aunt->mass*tmp/pow(square(tmp), 1.5);
	      PRL(aprt2);
	      tmp = thispos-sispos;
	      vec a12 = -sister->mass*tmp/pow(square(tmp), 1.5);
	      PRL(a12);
	    }


#endif

	}

    } else if (np != 0) {

	// (The following functions will have no effect if np = 0.)

	// Compute acceleration and jerk on this component due to the rest
        // of system, using the perturber list in pnode and including the
        // effect of the clump under the top-level node.

	calculate_partial_acc_and_jerk(top, top, get_parent(),
				       apert1, jpert1, p_dummy,
				       d_nn_sq, nn,
				       !USE_POINT_MASS,	    // explicit loop
				       pnode,		    // with list
				       this);		    // node to charge

	// Acceleration and jerk on other component due to rest of system:

	sister
	    ->calculate_partial_acc_and_jerk(top, top, get_parent(),
					     apert2, jpert2, p_dummy,
					     sister->d_nn_sq, sister->nn,
					     !USE_POINT_MASS,  // explicit loop
					     pnode,	       // with list
					     this);	       // node to charge

	// Note:  The first two calls to calculate_partial_acc_and_jerk pass
	// the d_nn_sq and nn from the hdyn, so the hdyn data are actually
	// updated by update_nn_coll.  The third call (below) passes local
	// variables, so it has no direct effect on d_nn_sq and nn.  (Test
	// afterwards if d_nn_sq and nn must be updated.)  The coll data
	// are always updated by these calls.

	// cerr << "Computed perturbation on " << get_parent()->format_label()
	//      << " using front end; np = " << np << endl << flush;
     }

    // Relative acceleration and jerk due to other (sister) component:

    vec a_2b, j_2b;
    real p_2b;
    hdyn *p_nn_sister;
    real d_min_sister = VERY_LARGE_NUMBER;

    vec d_pos, d_vel;

    if (!oldest_daughter && !sister->oldest_daughter) {

	// Make a concession to efficiency in the most common case...
	// *** May want to reconsider efficiency in other cases too. ***

	a_2b = j_2b = 0;
	p_2b = 0;
	d_pos = sister->pred_pos - pred_pos;
	d_vel = sister->pred_vel - pred_vel;
	accumulate_acc_and_jerk(sister, d_pos, d_vel, eps2,
				a_2b, j_2b, p_2b, d_min_sister);
	p_nn_sister = sister;

    } else

	calculate_partial_acc_and_jerk(get_parent(), get_parent(), this,
				       a_2b, j_2b, p_2b,
				       d_min_sister, p_nn_sister,
				       !USE_POINT_MASS,
				       NULL,		// no perturbers
				       this);		// node to charge

#if 0
    if (time >= 0.952164 && time <= 0.952273) {
	PRL(format_label());
	PRL(a_2b);
	PRL(apert1);
	PRL(apert2);
	PRL(d_pos);
	PRL(d_vel);
	PRL(pos);
	PRL(sister->pos);
	PRL(pred_pos);
	PRL(sister->pred_pos);
	PRL(sqrt(d_nn_sq));
	PRL(sqrt(d_min_sister));
    }
#endif

    if (np != 0) {

	real m2 = sister->get_mass();
	real pscale = m2 / (m2 + mass);

	// Note that perturbation_squared does *not* contain a slowdown factor.

	perturbation_squared = square(pscale * (apert1 - apert2))
	    				/ square(a_2b);

	a_2b += pscale * (apert1 - apert2) * get_kappa();
	j_2b += pscale * (jpert1 - jpert2);

    } else

	perturbation_squared = 0;

    set_acc_and_jerk_and_pot(a_2b, j_2b, p_2b);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(system_time) && T_DEBUG_LEVEL > 0 && name_is("23")) {
	cerr << "DEBUG: calculate_acc_and_jerk_on_low_level_node" << endl;
	PRI(4); PRC(p_2b); PRL(a_2b);
	PRI(4); PRL(apert1);
	PRI(4); PRL(apert2);
    }
#endif

#if 0
    if (is_low_level_node()) {
	PRC(time); PRL(format_label());
	PRL(abs(apert1 - apert2)/abs(jpert1 - jpert2));
    }
#endif

#if 0
    if (slow) {
	PRL(pred_pos);
	PRL(get_binary_sister()->pred_pos);
	PRL(a_2b);

	real m2 = get_binary_sister()->mass;
	vec sep = pred_pos * (1 + mass/m2);
	real r2 = sep*sep;
	vec a2 = -m2*sep / (r2*sqrt(r2));
	PRL(a2);

	PRL(pscale * (apert1 - apert2) * get_kappa());
    }
#endif

    if (d_min_sister < d_nn_sq) {

	nn = p_nn_sister;
	d_nn_sq = d_min_sister;

#if 0
	if (get_real_system_time() > 13.62265) {
	    cerr << "update_nn_coll(" << 999 << "): ";
	    PRL(format_label());
	    PRI(4); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
	}
#endif
    }

    // if (name_is("xxx")) {
    //     print_nn(this, 2);
    //     print_coll(this, 2);
    // }

    // The above procedure will fail to correctly determine the sister's
    // nearest neighbor.  Correct this here.  Note that, in this case,
    // the sister's nn may be a node (this), not a leaf; however, the
    // stored distance will accurately reflect the actual distance to the
    // nn leaf.  Too inconvenient and inefficient to routinely descend the
    // tree below this to (re)locate the nn leaf here -- however, we MUST
    // take this "feature" into account in new_sister_node (hdyn_tree.C)
    // when assessing the need for tree reconfiguration.

    if (d_min_sister < sister->get_d_nn_sq()) {

        sister->set_nn(this);
        sister->set_d_nn_sq(d_min_sister);

	// (OK to have sister->nn be a CM node.)
    }

    if ((nn == NULL) || (get_binary_sister()->get_nn() == NULL)) {
	cerr << "calculate_acc_and_jerk_on_low_level_node:  nn = NULL:\n";
	PRL(top_level->n_perturbers);
	pp3(this, cerr);
	pp3(top_level, cerr);
    }

    valid_perturbers = false;

//     if (system_time >= 44.15328 && system_time <= 44.1533) {
//       PRC(format_label()); PRL(pnode);
//       if (pnode) PRL(pnode->format_label());
//     }

    if (ALLOW_LOW_LEVEL_PERTURBERS && pnode) {

	// Create/revise the low-level perturber list(s).

	// Could revise list as we compute the force, but the distances used
	// should be from the center of mass, not from one of the components.
	// This may be a little less efficient, but it is much easier to code,
	// so keep it for now.  (Steve, 8/98)

	// Note that create_low_level_perturber_list() knows what to do if
	// only "low-only" perturber data are found in pnode.

	// Also continue down the hierarchy, creating lists for nodes whose
	// current lists are invalid.

        get_parent()->create_low_level_perturber_lists(false);

    }
}

local inline int n_leaves_or_unperturbed(hdyn* b)
{
    // Like n_leaves, but counting fully unperturbed center of mass nodes
    // as single particles.

    hdyn* od = b->get_oldest_daughter();
    if (od == NULL || (od->get_kepler() && od->get_fully_unperturbed())) {
	return 1;
    } else {
	int n = 0;
	for_all_daughters(hdyn, b, bb)
	    n += n_leaves_or_unperturbed(bb);
	return n;
    }
}

// expand_nodes:  expand center-of-mass nodes on a perturber list into
//		  single stars (i.e. leaves).
//
//		  NOTE from Steve, 8/20/98:  we (optionally) no longer
//		  expand fully unperturbed binaries.  This leaves the system
//		  in a slightly dangerous state, as it is possible that
//		  an unperturbed binary center of mass perturbing another
//		  binary may vanish before it can be resolved back into
//		  its components.  This is unlikely to occur dynamically,
//		  as the binary will first become perturbed and the
//		  components restored on the next center of mass time step.
//		  However, it can (and does!) happen.  To avoid this, all
//		  unperturbed binaries are expanded in integrate_unperturbed
//		  if the unperturbed motion is terminated.
//
//		  Nodes may be deleted as a result of stellar evolution (see
//		  merge_nodes), but the components are merged into the center
//		  of mass in that case.  (Nodes are also deleted as they
//		  escape, but this affects all particles, not just binaries.)
//
//		  Binaries undergoing partial unperturbed motion (pericenter
//		  reflection) are always expanded.

local inline bool expand_nodes(int &n, hdyn ** list, bool debug = false)
{
     if (debug)
	 cerr << endl << "in expand_nodes, n = " << n << endl;

    for (int i = 0; i < n; i++) {

 	if (debug)
 	    cerr << i << ".  " << list[i] << "  "
 		 << list[i]->format_label() << endl;

	if (!list[i]->is_leaf()) {

	    if (debug)
		cerr << "not a leaf" << endl;

	    int nl;

	    if (RESOLVE_UNPERTURBED_PERTURBERS)
		nl = list[i]->n_leaves();
	    else
		nl = n_leaves_or_unperturbed(list[i]);		// <-- new

	    if (debug) {
		pp3(list[i]->get_top_level_node());
		PRL(nl);
	    }

	    if (n + nl - 1 > MAX_PERTURBERS)
		return false;

	    // Make room for the expansion.

	    int j;
	    for (j = n - 1; j > i; j--)
		list[j + nl - 1] = list[j];

	    // Add the new nodes to the list.

	    j = 0;
	    hdyn *tmp = list[i];

	    if (debug)
		PRL(tmp->format_label());

	    if (RESOLVE_UNPERTURBED_PERTURBERS) {
		for_all_leaves(hdyn, tmp, b) {
		    list[i + j] = b;
		    j++;
		}
	    } else {

		for_all_leaves_or_unperturbed(hdyn, tmp, b) {	// <-- new
		    list[i + j] = b;
		    if (debug)
			cerr << "added " << b->format_label() << endl;
		    j++;
		}
	    }

	    n += nl - 1;
	    i += nl - 1;

	    if (debug && j != nl)
		cerr << "expand_nodes -- oops... "
		     << j << " != " << nl << endl;
	}
    }

    return true;
}

//-----------------------------------------------------------------------------
//  NOTE: calculate_acc_and_jerk_on_top_level_node is now split into three
//  pieces to allow use of GRAPE hardware.
//-----------------------------------------------------------------------------

void hdyn::calculate_acc_and_jerk_on_top_level_node(bool exact)
{
#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  calculate_acc_and_jerk_on_top_level_node for "
	     << format_label() << endl;
    }
#endif

    top_level_node_prologue_for_force_calculation(exact);

    if (!exact) {
	top_level_node_real_force_calculation();
	top_level_node_epilogue_force_calculation();
    }

}


void hdyn::top_level_node_prologue_for_force_calculation(bool exact)
{
#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  top_level_node_prologue_for_force_calculation for "
	     << format_label() << endl;
    }
#endif

    hdyn *root = get_root();
    d_coll_sq = VERY_LARGE_NUMBER;
    coll = NULL;

    if (exact) {

	// Perform the entire O(N) force calculation on the front end...

	n_perturbers = 0;
	valid_perturbers = false;

	nn = NULL;
	d_nn_sq = VERY_LARGE_NUMBER;

	clear_interaction();
	calculate_partial_acc_and_jerk(root, root, this,
				       acc, jerk, pot, d_nn_sq, nn,
				       !USE_POINT_MASS,
				       NULL,		// no perturbers
				       this);		// node to charge

    } else if (is_parent()) {

	// Set up computation of perturber list.

        // new_perturber_list();	// don't do this here -- only when
					// we are sure we will recompute the
					// perturbers

	perturbation_radius_factor
		= define_perturbation_radius_factor(this, gamma23);

    }
}

int hdyn::top_level_node_real_force_calculation()
{
    // *************************************************************
    // **** This function is NEVER called if GRAPE is available ****
    // **** -- it is replaced by grape_calculate_acc_and_jerk.  ****
    // *************************************************************

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  top_level_node_real_force_calculation for "
	     << format_label() << endl;
    }
#endif

    // Special treatment of traversal of top level only (for
    // efficiency, and for GRAPE implementation).

    // Calculate forces first in the POINT-MASS approximation,
    // without using any perturber list, and construct a new list,
    // in the case of CM nodes.

    // Should NOT be called in exact case.

    hdyn * root = get_root();

    if (is_parent()) {

	if (get_oldest_daughter()->slow)
	    clear_perturbers_slow_perturbed(this);

	new_perturber_list();		// effectively choose pert_freq = 1

        n_perturbers = 0;
        valid_perturbers = true;
    }

    clear_interaction();
#ifdef USEMPI
    double t0=MPI_Wtime();
#endif
    int n_top = flat_calculate_acc_and_jerk(root, is_parent());
#ifdef USEMPI
#if 0
    double t1=MPI_Wtime()-t0;
    double *timers=new double[mpi_nprocs];
    MPI_Gather(&t1,1,MPI_DOUBLE,timers,1,MPI_DOUBLE,0,mpi_communicator);
    if (mpi_myrank == 0)
    {
      cerr << "flat_calculate_acc_and_jerk: ";
      for (int i=0; i<mpi_nprocs; i++)
	cerr << timers[i] << " ";
      cerr << endl;
    }
    delete[]timers;
#endif
#endif

    if (nn == NULL) {
        // pretty_print_node(cerr); cerr << " nn NULL after flat " << endl;
    }

    return n_top;
}

void hdyn::top_level_node_epilogue_force_calculation()
{
    static const char *func = "top_level_node_epilogue_force_calculation";

#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  " << func << " for " << format_label() << endl;
    }
#endif

#if 0
    if (time > 13.62265) {
	cerr << func << "(1): " << endl;
	PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
    }
#endif

    // Should NOT be called in exact case.

    // Take care of the case where the nearest neighbor is
    // a complex node, since otherwise nn may disappear.
    // (Jan 21, 1998, JM and SPZ)

    // Exclude the case of nn = this (possible when GRAPE is used.)
    // (Sep 9, 1998, SLWM)

    if (nn && nn != this) {
	while (nn->is_parent()) {
	    // cerr << "nn correction for "; pretty_print_node(cerr);
	    // cerr << "and  "; nn->pretty_print_node(cerr);

	    nn = nn->get_oldest_daughter();

	    // cerr << " as   ";   nn->pretty_print_node(cerr);
	    // cerr << endl;
	}
    }

#if 0
    if (time > 13.62265) {
	cerr << func << "(2): " << endl;
	PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
    }
#endif

    hdyn* od = get_oldest_daughter();
    if (!od) return;

    // Apply center-of-mass correction to binary node.

    hdyn * root = get_root();
	
    // if (is_parent()) {
    //     print_label(cerr); cerr << " nn before epilogue ";
    //     n->print_label(cerr); cerr << " ";
    //     RL(d_nn_sq);
    // }

    // Perturber list was computed in "real_force_calculation" or by GRAPE.

    // If the perturber list has overflowed and no perturber list is available,
    // we will have to use the entire tree (below).  We'd prefer to avoid this
    // if possible, so the GRAPE version of the code provides a "low-only" list
    // -- a subset of the complete list, for the MAX_PERTURBERS closest nodes.
    // It is intended for use in computing low-level perturbers, but it is as
    // large as possible, so it can be used here as an approximate means of
    // avoiding a potentially expensive calculation on the front end.  It is
    // constructed so that it should not overflow when CM nodes are expanded
    // into components, so the "low" list should always be usable.  Note that
    // calculate_partial_acc_and_jerk() knows what to do if sent a node with
    // only "low" data available.
    //							(Steve, 3/03)

    if (n_perturbers > MAX_PERTURBERS || !valid_perturbers) {
	valid_perturbers = false;
	kc->perturber_overflow++;
    }

    // Use whatever data we have to try to avoid an expensive calculation.
    // Note that the decision to use the low-level perturber list for the
    // top-level node as an efficiency measure (with errors) is completely
    // independent of the use of the list for low-level nodes.
	
    if (valid_perturbers && n_perturbers > 0
	|| valid_perturbers_low && n_perturbers_low > 0) {
	
	// Expand a copy of the perturber list to see in advance if the list
	// will overflow.  Shrink the perturber list until the expansion is
	// successful.

	hdynptr plist[MAX_PERTURBERS];
	int np = n_perturbers, nlist;
	if (!valid_perturbers) np = n_perturbers_low;
	nlist = np;
	for (int i = 0; i < np; i++)
	    plist[i] = perturber_list[i];

	real dmax2 = 0;
	real d2[MAX_PERTURBERS];

	while (!expand_nodes(nlist, plist)) {

	    // Expanded list is too long.  Reduce the length and retry.

	    if (dmax2 == 0) {

		kc->perturber_overflow++;

		for (int i = 0; i < np; i++) {
		    d2[i] = square(perturber_list[i]->get_pred_pos()
				   - pred_pos);
		    if (d2[i] > dmax2) dmax2 = d2[i];
		}
	    }

	    // Shrink the list.

	    dmax2 *= 0.95;

	    int nskip = 0;
	    for (int i = 0; i < np; i++) {
		if (nskip > 0) {
		    perturber_list[i-nskip] = perturber_list[i];
		    d2[i-nskip] = d2[i];
		}
		if (d2[i] > dmax2) nskip++;
	    }
	    np -= nskip;

	    valid_perturbers = false;
	    valid_perturbers_low = true;	// by construction

	    // Make a new copy and continue.

	    nlist = np;
	    for (int i = 0; i < np; i++)
		plist[i] = perturber_list[i];
	}

	// Original list is perturber_list, of length np;
	// expanded list is plist, of length nlist.

	// Use the (partial) perturber list to correct the center-of-mass force.
	// First, subtract out the center-of-mass contribution due to perturbers
	// (note that, at this point, the perturber list contains only CMs):
	
	vec a_cm, a_p, j_cm, j_p;
	real p_p;
	calculate_partial_acc_and_jerk(this, this, this,
				       a_cm, j_cm, p_p, d_nn_sq, nn,
				       USE_POINT_MASS,	    	    // explicit
				       this,	   		    // loop
				       this);

#if 0
	if (time > 13.62265) {
	    cerr << func << "(3a): " << endl;
	    PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
	}
#endif

	// (Note: 'this' is OK here because this is a top-level node.)

	// acc -= a_cm;
	// jerk -= j_cm;

	pot -= p_p;

	// Replace nodes by components on the perturber list and
	// add in contribution due to components.

	for (int i = 0; i < nlist; i++)
	    perturber_list[i] = plist[i];

	if (valid_perturbers)
	    n_perturbers = nlist;
	else if (valid_perturbers_low) {
	    n_perturbers = MAX_PERTURBERS + 1;
	    n_perturbers_low = nlist;
	}

	// Include component forces in the CM force (explicit loop).
	
	calculate_partial_acc_and_jerk(this, this, this,
				       a_p, j_p, p_p, d_nn_sq, nn,
				       !USE_POINT_MASS,
				       this,
				       this);

#if 0
	if (time > 13.62265) {
	    cerr << "top_level_node_epilogue_force_calculation(3b): "
		 << endl;
	    PRC(format_label()); PRC(nn->format_label());
	    PRL(sqrt(d_nn_sq));
	}
#endif

	// Include the "slow" term in correcting acc and jerk.  The
	// back reaction is handled in correct_acc_and_jerk.

	// In the case of slow binary motion, on exit from this function,
	// acc and jerk remain computed in the point-mass approximation;
	// the perturbative corrections are stored separately.

	// This may not be quite right if one of the perturbers is
	// itself a slow binary.  However, the procedure followed in
	// correct_and_update() should properly treat the dominant
	// term due to the internal motion of this node.


// 	if (time >= xreal(2296, 3651856000000000000)) {
// 	    cerr << endl << "in top_level_node_epilogue_force_calculation"
// 		 << endl;
// 	    int p = cerr.precision(HIGH_PRECISION);
// 	    PRL(format_label());
// 	    pp3(get_top_level_node());
// 	    cerr << endl;
//	    cerr.precision(p);
// 	}


	if (!od->slow) {

	    // Normal correction in non-slow case.

	    acc += a_p - a_cm;
	    jerk += j_p - j_cm;

	} else {

	    // Keep the CM approximation and save the perturbation
	    // for use in correct_and_update().

	    od->slow->set_acc_p(a_p - a_cm);
	    od->slow->set_jerk_p(j_p - j_cm);

	    // The new perturber list is valid and intact.  Update
	    // the slow_perturbed lists of the top_level_nodes of all
	    // perturbers.

	    for (int j = 0; j < n_perturbers; j++) {
		hdyn *pert_top = perturber_list[j]->get_top_level_node();
#if 0
		cerr << "adding " << format_label();
		cerr << " to slow_perturbed list of "
		     << pert_top->format_label()
		     << endl;
#endif
		pert_top->add_slow_perturbed(this, diag->slow_perturbed);
	    }
	}

	pot += p_p;
	

	// if (nn == NULL) {
	//     pretty_print_node(cerr);
	//     cerr << " nn NULL after add back "<<endl;
	// }

#if 0
	if (time > 13.62265) {
	    cerr << "top_level_node_epilogue_force_calculation(4): " << endl;
	    PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
	}
#endif

    }

    // Recompute the force exactly if the perturber list has overflowed
    // (shouldn't happen, now, at least in GRAPE-6 code...).

    if (!valid_perturbers && !valid_perturbers_low) {

#if 0
	cerr << func << ": calculate_partial_acc_and_jerk for "
	     << format_label() << " at time " << system_time << endl;
#endif

	n_perturbers = -1;
	calculate_partial_acc_and_jerk(root, root, this,
				       acc, jerk, pot, d_nn_sq, nn,
				       !USE_POINT_MASS,
				       NULL,		// no perturbers
				       this);
    }

#if 0
    if (time > 13.62265) {
	cerr << "top_level_node_epilogue_force_calculation(5): " << endl;
	PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
    }
#endif

    // NOTE:  On exit, if valid_perturbers is true, and
    //        RESOLVE_UNPERTURBED_PERTURBERS is true, then
    //	      perturber_list contains *leaves only.*
    //
    //	      If valid_perturbers is false, then the force on 'this'
    //	      node has been computed exactly; no correction is needed.

    // if (nn == NULL) {
    //     pretty_print_node(cerr);
    //     cerr << " nn NULL after corrections "<<endl;
    // }

    // if (is_parent()) {
    //     print_label(cerr); cerr << " nn after epilogue ";
    //     nn->print_label(cerr); cerr << " ";
    //     PRL(d_nn_sq);
    //     PRC(valid_perturbers); PRL(n_perturbers);
    // }

    // Finally, if the perturber list is valid, propagate the information
    // to lower levels.

    if (ALLOW_LOW_LEVEL_PERTURBERS
	&& (valid_perturbers || valid_perturbers_low))
	create_low_level_perturber_lists(false);	// false ==> all
}


//-----------------------------------------------------------------------------
// calculate_acc_and_jerk:  Generic routine to calculate the acceleration
// of a node.  Call appropriate functions depending on whether or not the
// node is at the top level of the tree.
//-----------------------------------------------------------------------------

void hdyn::calculate_acc_and_jerk(bool exact)
{
#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "calculate_acc_and_jerk for "
	     << format_label() << endl;
    }
    dbg_message("calculate_acc_and_jerk", this);

    hdyn *sister = get_binary_sister();
#else
    hdyn *sister = (hdyn*)younger_sister;		// assume binary tree
    if (!sister) {
	sister = (hdyn*)elder_sister;
	if (!sister)
	    return;
//	    err_exit("calculate_acc_and_jerk_on_low_level_node: no sister!");
    }
#endif

//   if (system_time >= 44.15328 && system_time <= 44.1533) {
//     cerr << "calculate_acc_and_jerk for " << format_label() << ",  ";
//     PRL(exact);    
//   }

    d_coll_sq = VERY_LARGE_NUMBER;
    coll = NULL;
    if (is_low_level_node())
	sister->set_d_coll_sq(VERY_LARGE_NUMBER);

    if (is_top_level_node())
	calculate_acc_and_jerk_on_top_level_node(exact);
    else
	calculate_acc_and_jerk_on_low_level_node();

    // NOTE: On return, we must correct the acc and jerk on:
    //
    //       (1) a top-level C.M. node due to other nodes not on its
    //           perturber list, but on whose perturber lists it lies.
    // and
    //       (2) a single top-level node due to any C.M.s on whose perturber
    //           lists it lies.

    if (nn == NULL) {
         // pretty_print_node(cerr); cerr << " nn NULL after calc " << endl;
    }
}


//=========================================================================
//
// correct_acc_and_jerk:  By construction, the acc and jerk on a top-level
//			  node bi_top will be computed incorrectly if it is a
//			  perturber of another top-level node bj that is not
//			  itself a perturber of bi_top.  (Note that if bi_top
//			  is a leaf, then this latter criterion is necessarily
//			  satisfied.)  The reason that the force calculation
//			  requires correction is that it is written to enable
//			  efficient computation of the top-level interactions
//			  with the GRAPE.  The calculation is first performed
//			  in the point-mass approximation.  Then the acc and
//			  jerk on bi_top due to its perturbers are corrected
//			  to resolve bi_top into its components.  However, if
//			  bi_top is a perturber of bj and bj is not a perturber
//			  of bi_top, no correction has yet been applied.
//			  Apply it here.
//
//=========================================================================

// apply_correction: bi is a top-level node on the integration list,
//		     whose interaction with bj has not been correctly
//		     calculated.  Include the missing (tidal) terms here.

local inline void apply_correction(hdyn * bj, hdyn * bi)
{
    vec a_c = 0, j_c = 0;
    real p_c = 0;
    vec a_p = 0, j_p = 0;
    real p_p = 0;
    real dum_d_sq = VERY_LARGE_NUMBER;
    hdyn *dum_nn = NULL;

#ifdef KIRA_DEBUG
    if (bj->get_kira_diag()->correct) {
	cerr << "correcting " << bi->format_label();
	cerr << " for " << bj->format_label() << " at "
	     << bi->get_system_time() << endl;
    }
#endif

    // cerr << "before correction" << endl;
    // PRL(bi->get_acc());
    // PRL(bj->get_acc());

    bi->calculate_partial_acc_and_jerk(bj, bi->get_parent(), bi,
				       a_c, j_c, p_c, dum_d_sq, dum_nn,
				       USE_POINT_MASS,
				       NULL,		// no perturbers
				       bi);
    bi->calculate_partial_acc_and_jerk(bj, bi->get_parent(), bi,
				       a_p, j_p, p_p, dum_d_sq, dum_nn,
				       !USE_POINT_MASS,
				       NULL,		// no perturbers
				       bi);

    // Separate out the "perturbative" components of acc and jerk if bj is
    // a slow binary.  See also top_level_node_epilogue_force_calculation().

    kira_diag *kd = bi->get_kira_diag();

    hdyn *od = bj->get_oldest_daughter();
    if (od && od->get_slow()) {

	// Save the perturbative acc and jerk for use in the corrector.
	// Don't update the center-of-mass quantities yet.  Because the
	// correction is different for different slowdown factors, it
	// is necessary to save the perturbations from different slow
	// binaries separately, hence the unpleasant bookkeeping...

	// *** Need to monitor these lists to ensure that they don't
	// *** get out of control!

	// Procedure:	(1) See if bj is on the list of slow binaries
	//		    associated with node bi, and check that
	//		    its kappa hasn't changed.
	//		(2) If it is, modify the entry.  If not, create
	//		    a new entry (but retain the center of mass
	//		    approximation for this step).
	//		(3) Also clean up the list by removing links to
	//		    invalid nodes.

	slow_perturbed *s = bi->find_slow_perturbed(bj);
	if (!s) {

	    // Shouldn't happen...

	    cerr << "apply_correction: adding " << bj->format_label();
	    cerr << " to slow_perturbed list of " << bi->format_label()
		 << endl;
	    bj->print_perturber_list(cerr);

	    s = bi->add_slow_perturbed(bj, kd->slow_perturbed);
	}

	if (s) {
	    s->set_acc_p(a_p - a_c);
	    s->set_jerk_p(j_p - j_c);
	}

	// Failure (!s) shouldn't occur here, but if it does, the default
	// is that we simply proceed in the center of mass approximation.

	bi->check_slow_perturbed(kd->slow_perturbed);	// clean up if necessary

#if 0
	cerr << "apply_correction: slow_perturbed treatment for bi = "
	     << bi->format_label();
	cerr << " and bj = " << bj->format_label() << endl;
	PRL(bi->count_slow_perturbed());
#endif

    } else {

	// Normal correction of acc and jerk.

	bi->inc_acc(a_p - a_c);
	bi->inc_jerk(j_p - j_c);

	bi->check_slow_perturbed(kd->slow_perturbed);	// clean up if necessary
    }

    bi->inc_pot(p_p - p_c);

    // Note: changes to indirect_force counter are handled at lower levels.

    // We should update nn etc here, since otherwise
    // perturbers of a perturbed binary see the binary
    // as the nearest neighbor, which is unsafe.
    // The following correction added Jan 21 1998.

    if (bi->get_d_nn_sq() > dum_d_sq) {

	bi->set_d_nn_sq(dum_d_sq);
	bi->set_nn(dum_nn);

#if 0
	if (bi->get_real_system_time() > 13.6235) {
	    cerr << "nn correction for ";   bi->pretty_print_node(cerr);
	    cerr << " and ";   bj->pretty_print_node(cerr);
	    cerr << " as  ";   bi->get_nn()->pretty_print_node(cerr);
	    cerr << "  new distance = " << bi->get_d_nn_sq() <<endl;
	}
#endif
    }

    bi->get_kira_counters()->force_correction++;

    // cerr << "after correction" << endl;
    // PRL(bi->get_acc());
    // PRL(bj->get_acc());
}


//-----------------------------------------------------------------------------
//
// need_correction:  return a pointer to a node whose force must be corrected.
//
//     bj: a perturbed top-level CM node
//     bi: one of bj's perturbers
//     t:  present system time
//
// The requirement is as follows:
//
//     a) the actual pointer returned is the top_level_node of bi (btop)
//     b) bi must be the left-most leaf of the tree under btop, to
//        guarantee uniqueness (i.e. that the correction is applied only once
//     c) btop must be in the present integration block
//     d) bj must not be on the perturber list of btop
//
//
// ***** OLD CODE -- goes with old version of correct_acc_and_jerk. *****
//
//-----------------------------------------------------------------------------

local inline hdyn *need_correction(hdyn * bj, hdyn * bi)
{
    hdyn *btop = bi->get_top_level_node();

    if (!btop->is_on_integration_list()) return NULL;

    // Now we know that btop is on the integration list.
    // If it is a single star, it needs correction.

    if (btop ->is_leaf()) return btop;

    hdyn *b_oldest = btop;
    hdyn *bret = NULL;

    // btop is the top-level node of a clump requiring correction.
    // Test if bi is the left-most ("oldest") leaf of the clump, in
    // order to guarantee that btop is corrected only once.

    // Find the left-most leaf.

    while (b_oldest->is_parent())
	b_oldest = b_oldest->get_oldest_daughter();

    // Note that bi can be a top-level node in the case
    // of correction from a node without valid perturber
    // list.

    if ( (b_oldest == bi) || (btop == bi)) {

	bret = btop;

	// No need to apply a correction if some component of bj is a
	// perturber of btop (interaction already properly calculated).

	// Top-level node btop is on the perturber list of node bj.
	// Return NULL iff bj is a perturber of btop.

	// Well, what if perturber list of btop is INVALID?  This in practice
	// *should not* occur unless previous force calculation is exact...
	// ...or if perturber-list overflow has occurred (Steve 7/98).

	if (btop->get_valid_perturbers() || btop->get_valid_perturbers_low()) {
	    int np = btop->get_n_perturbers();
	    if (!btop->get_valid_perturbers()) np = btop->get_n_perturbers_low();
	    for (int j = 0; j < np; j++)
		if (btop->get_perturber_list()[j]->get_top_level_node() == bj)
		    bret = NULL;
	}

	// Code modified by Steve 7/98

	else
	    bret = NULL;	// btop force should be correct in this case...

    }

    return bret;
}


//-----------------------------------------------------------------------------
//
// correct_acc_and_jerk:  OLD version.  Inefficient, but works...
//
// To avoid the use of inverse perturber lists, we identify nodes bi_top
// and bj by scanning across the entire top-level (bj), then checking for
// perturbers (bi).  Correction, if necessary, is applied to bi_top, the
// top-level node of bi.  Membership on the integration list is
// determined by checking the flag on_integration_list(), so there is no
// need for the list itself to be provided to this routine.
//
// Note from Steve, 9/98.  This is a very inefficient function, as it checks
// all the perturber lists of all nodes, even though in practice very few
// corrections are actually applied.  The new version below uses the list
// of perturbed binaries to speed things up.
//
//-----------------------------------------------------------------------------

// Static data:

static int work_size = 0;
static hdyn ** nodes = NULL;
static int nnodes = -1;

// Allow possibility of cleaning up if necessary:

void clean_up_hdyn_ev() {if (nodes) delete [] nodes;}

void correct_acc_and_jerk(hdyn * root,		// OLD!
			  bool& reset)
{

    if (nnodes == -1) reset = true;

    if (reset) {

	// initialize the array of node pointers

        // cerr << "correct_acc_and_jerk, reset\n";

	int n = 0;
	for_all_daughters(hdyn, root, bb)
	    if (bb->is_parent()) n++;
	if  (work_size < n) {
	    nodes = new hdynptr[n*2+10];	// DEC C++ doesn't like (hdyn*)
	    work_size = n*2+10;
	}
	n = 0;
	for_all_daughters(hdyn, root, bbb) {
	    if (bbb->is_parent()) {
		nodes[n] = bbb;
		n++ ;
	    }
	}
	nnodes = n;
    }

    reset = false;
    for (int j = 0; j < nnodes; j++) {

	hdyn * bj = nodes[j];

	// for_all_daughters(hdyn, root, bj) {

	if (!bj->is_valid()) {

	    // Invalid node ==> should not occur; better reconstruct
	    // the internal list if found.

	    cerr << "correct_acc_and_jerk: warning: invalid node at time "
		 << root->get_system_time() << endl;
	    reset = true;

	} else if (bj->is_parent()
		   && bj->get_oldest_daughter()->get_kepler() == NULL) {

	    hdyn *bi_top;

	    // bj is a perturbed top-level CM node.

	    if (bj->get_valid_perturbers()) {

		// Look in bj's perturber list for something to correct.

		for (int i = 0; i < bj->get_n_perturbers(); i++) {

		    hdyn *bi = bj->get_perturber_list()[i];

		    // Flag an invalid perturber...

		    if (!bi->is_valid()) {

			cerr << "warning: correct_acc_and_jerk: "
			     << "invalid perturber #" << i
			     << " (" << bi << ")"
			     << "         for " << bj->format_label()
			     << " at time " << root->get_system_time()
			     << endl
			     << "         forcing recomputation "
			     << " of perturber list"
			     << endl;

			// Overkill:
			//
			// pp3(root);
			// pp3(bi);
			// exit(1);

			bj->set_valid_perturbers(false);
			break;

		    } else {

			if ((bi_top = need_correction(bj, bi))
			    != NULL) {
			    apply_correction(bj, bi_top);
			}
		    }
		}

	    } else {

		// Look at all top-level nodes for something to correct.

		// Note (J.M. 96-Aug-5):
		// Actually, one needs to loop over particles in current
		// block only...  Not implemented yet.

		for_all_daughters(hdyn, root, bi) {
		    if (bi != bj) {
			if ((bi_top = need_correction(bj, bi))
			      != NULL) {
			    apply_correction(bj, bi_top);
			}
		    }
		}
	    }
	}
    }
}


//-----------------------------------------------------------------------------
//
// check_and_apply_correction:	helper function for use with new version
//				of correct_acc_and_jerk.
//
//-----------------------------------------------------------------------------

local inline void check_and_apply_correction(hdyn * bj, hdyn * bi)
{
    // On entry, bj is a perturbed binary.  Node bi is on bj's perturber
    // list, if it exists.  Otherwise, bi is a node in the current time
    // step block.  We want to check if the top-level node of bi (btop)
    // needs correction.  Criteria are:
    //
    //	    1. btop is on the integration list
    //	    2. btop is not perturbed by bj

    hdyn *btop = bi->get_top_level_node();

    if (btop->is_on_integration_list()) {

	// Node btop is on the integration list and effectively perturbs bj.

	// If btop is a single star or unperturbed, it needs correction.
	// If not, apply more checks.

	// Hmmm.  It is possible, although it is not supposed to happen,
	// that the individual components of an unperturbed binary or
	// multiple system may appear on a perturber list even when
	// RESOLVE_UNPERTURBED_PERTURBERS is set false.  We take the
	// defensive position here of checking the components of all
	// binaries, even unperturbed systems.

	if (btop->is_parent()

//	    && !btop->get_oldest_daughter()->get_kepler()) {  // check moved
							      // to (*) below
	    ) {

	    // Since perturbers are determined in the CM approximation
	    // and subsequently are expanded into component leaves, all
	    // leaves should be on the list if one is.  Test if bi is
	    // the left-most ("oldest") leaf of the clump, in order to
	    // guarantee that btop is corrected only once.

	    // Find the left-most leaf.

	    hdyn *b_oldest = btop;

	    while (b_oldest->is_parent())
		b_oldest = b_oldest->get_oldest_daughter();

	    // Note that bi can be a top-level node in the case of correction
	    // from a node without valid perturber list.

	    // cerr << "b_oldest = " << b_oldest->format_label() << endl;

	    if (b_oldest != bi && btop != bi) return;

	    if (!btop->get_oldest_daughter()->get_kepler()) {  //  (*) <---

		// Node btop is perturbed.  See if some component of bj
		// is a perturber of btop.

		// What if perturber list of btop is INVALID?  This in
		// practice *should not* occur unless previous force
		// calculation was exact, or if perturber-list overflow
		// has occurred (Steve 7/98).

		if (!btop->get_valid_perturbers()) return;

		for (int j = 0; j < btop->get_n_perturbers(); j++)
		    if (btop->get_perturber_list()[j]->get_top_level_node()
				== bj)
			return;
	    }
	}

	// Need to apply a correction to btop.

	apply_correction(bj, btop);

    }
}


//-----------------------------------------------------------------------------
//
// NEW version of correct_acc_and_jerk uses perturbed binary list to
// speed up identification of systems that need correction.
//
// Only calling function is calculate_acc_and_jerk_for_list() (kira_ev.C).
//
//-----------------------------------------------------------------------------

void correct_acc_and_jerk(hdyn** next_nodes,	// NEW
			  int n_next)
{
    if (n_next < 1) return;

    // Look on the list of top-level perturbed nodes for nodes
    // perturbed by a node in the current time step block.

    hdyn* root = next_nodes[0]->get_root();
    hdyn** perturbed_list = root->get_perturbed_list();

    if (!perturbed_list) {

	// No perturbed binary list.  Revert to the old code!

	bool reset = true;
	correct_acc_and_jerk(root, reset);
	return;
    }

    int n_perturbed = root->get_n_perturbed();

    for (int j = 0; j < n_perturbed; j++) {

	hdyn * bj = perturbed_list[j];  // bj is a perturbed top-level CM node

	if (bj->is_valid()			// should not be necessary...
	    && bj->get_oldest_daughter()	// also should not be necessary

	    && !bj->get_oldest_daughter()	// could be partially
	    	  ->get_kepler()) {		//     unperturbed motion

	    if (bj->get_valid_perturbers()) {

		// Look in bj's perturber list for something to correct.

		for (int i = 0; i < bj->get_n_perturbers(); i++) {

		    hdyn *bi = bj->get_perturber_list()[i];

		    // Flag an invalid perturber...

		    if (!bi->is_valid()) {

			cerr << "warning: correct_acc_and_jerk: "
			     << "invalid perturber #" << i
			     << " (" << bi << ")" << endl;
			cerr << "         for " << bj->format_label()
			     << " at time " << root->get_system_time()
			     << endl
			     << "         forcing recomputation "
			     << " of perturber list"
			     << endl;

			bj->set_valid_perturbers(false);
			break;

		    } else {

			// bi is a valid node on bj's perturber list.

			check_and_apply_correction(bj, bi);
		    }
		}

	    } else {

		// Look at all top-level nodes associated with the current
		// block for something to correct.

		for (int i = 0; i < n_next; i++)
		    if (next_nodes[i] != bj
			&& next_nodes[i]->is_top_level_node())
			check_and_apply_correction(bj, next_nodes[i]);

	    }
	}
    }
}


//=============================================================================
//  driver function for a one-time-step integration of this node
//=============================================================================

//-----------------------------------------------------------------------------
//  integrate_node -- advance this node from time to system_time.
//
//                    Note that the resulting step may or may not be equal
//                    to the node's timestep (new feature, Steve, 4/03).
//
//                    We don't predict here -- must be taken care of by
//                    the calling function.  In that case, the only use 
//                    of system_time is in correct_and_update(), in
//                    determining the actual step.
//
//		      Currently, this function is called only by
//		      synchronize_node().
//
//-----------------------------------------------------------------------------

void hdyn::integrate_node(bool exact,			// default = true
			  bool integrate_unperturbed,	// default = true
			  bool force_unperturbed_time)	// default = false

// From Steve: force_unperturbed_time can cause an unperturbed binary
// to be advanced to an undesirable phase.

{
#ifdef KIRA_DEBUG
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "integrate_node for " << format_label()
	     <<" at time " << time  + timestep << endl;
	// pretty_print_node(cerr);
    }
#endif

    if (time == system_time) return;

    bool local_debug = false;	// time > 0.614;

    if (!kep) {

	clear_interaction();
	calculate_acc_and_jerk(exact);
	if (exact) set_valid_perturbers(false);

	if (tidal_type && is_top_level_node())
	    get_external_acc(this, pred_pos, pred_vel,
			     pot, acc, jerk);

	if (local_debug) {PRC(1); PRL(timestep);}

	correct_and_update();		// sets time, timestep,...
					// note: no retry on error

	if (local_debug) {PRC(2); PRL(timestep);}

	store_old_force();

	// Note that old_acc = acc at the end of a step.

    } else if (integrate_unperturbed) {

	if (eps2 == 0) {

	    bool reinitialize;
	    integrate_unperturbed_motion(reinitialize,
					 force_unperturbed_time);

	    if (reinitialize)
		cerr << endl
		     << "integrate_node: received reinitialization flag "
		     << "during synchronization..." << endl;

	}
	else
	    err_exit(
	        "integrate_node: unperturbed binary with non-zero softening");
    }
}

void synchronize_tree(hdyn * b)		// update all top-level nodes
{					// better to use GRAPE if we can

    // No prediction is performed here -- up to the calling function.
    // Also, the step will take us to system_time, which also should
    // be set on entry.

    // cerr << "synchronize_tree for " << b->format_label()
    //      << " at time " << b->get_system_time() << endl;

    // Added by Steve (11/04) to avoid losing perturber list
    // information when a multiple clump is synced.  Not clear if
    // exact = true was ever necessary -- watch this!!  See also
    // contrary note in synchronize_node() -- to be resolved...

    bool exact = false;

    if (!b->is_root())
	b->synchronize_node(exact);

    for_all_daughters(hdyn, b, bb)
	synchronize_tree(bb);
}

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  src/dyn/evolve/hdyn_ev.C
//  +---------------+                     +------------------------------//
//=========================
