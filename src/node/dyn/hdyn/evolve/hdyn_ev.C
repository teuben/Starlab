
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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
#include "hdyn_inline.C"

#define INITIAL_STEP_LIMIT	0.0625	// Arbitrary limit on the first step
#define USE_POINT_MASS		true

// General debugging flag:

static bool dbg = false;


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
    hdyn *sister = bi->get_binary_sister();

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
	bi->set_jerk(vector(0,0,0));
	sister->set_jerk(vector(0,0,0));
    }

    sister->set_pot(-factor * bi->get_pot());
    sister->store_old_force();

    sister->set_perturbation_squared(bi->get_perturbation_squared());
}

//-----------------------------------------------------------------------------
// synchronize_node -- force the time of a node to a specified value,
//                     by integrating it forward in time to system_time.
//-----------------------------------------------------------------------------

void hdyn::synchronize_node()
{
    bool do_diag = (diag->ev && diag->check_diag(this));

    if (do_diag)
	cerr << "Enter synchronize_node for " << format_label() << endl;

    if (time == system_time) return;

    if (is_low_level_node()
	&& (get_parent()->get_oldest_daughter() != this)) {
	if (do_diag)
	    cerr << "Low-level node: "<<format_label()<<endl<<flush;
	get_elder_sister()->synchronize_node();
	update_binary_sister(get_elder_sister());
	return;
    }

    real old_timestep = timestep;
    timestep = system_time - time;

    if (timestep < 0.0) {

	if (do_diag)
	    cerr << "synchronize_node: negative timestep..." << system_time
		 << " " << time << flush << endl;
	put_node(cerr, *this, options->print_xreal);
	exit(-1);
    }

    if (do_diag) {
	cerr << "Synchronize node " << format_label()
	     << " at time " << system_time << endl
	     <<" ...before integrate_node " << endl << flush;
	if (is_low_level_node()) pp3(get_parent());
	cerr << flush;
    }

    // cerr << "Energy before "; print_recalculated_energies(get_root(), 1);

    // NOTE: the "false" here means that we do NOT attempt to
    //	     synchronize unperturbed binaries...

    integrate_node(get_root(), false, true);

    // cerr << "Energy after "; print_recalculated_energies(get_root(), 1);

    if (do_diag)
	cerr <<"After integrate_node "<<endl;

    // NOTE: Timestep adjustment below is not clean...

//    cerr << format_label() << " "; PRC(old_timestep);

    while (fmod(time, old_timestep) != 0.0)
	old_timestep *= 0.5;
    timestep = old_timestep;

//    PRL(timestep);

    if (is_low_level_node())
	update_binary_sister(this);

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
//               32             \   | jerk |       | acc |    /         ~
//
//  The first expression gives the time scale on which the acceleration is
//  changing significantly.
//  The second expression gives an estimate for the free fall time
//  in the case of near-zero initial velocity.
//-----------------------------------------------------------------------------

void hdyn::set_first_timestep(real additional_step_limit) // default = 0
{
    real eta_init = eta / 32;	// initial accuracy parameter

    real j2 = jerk * jerk;
    real a2 = acc * acc;

    real dt_adot = eta_init/4;	// time step based on rate of change of acc
    real dt_ff   = eta_init/4;	// time step based on free-fall time
    real dt;			// minimum of the two

    if (j2 > 0)
	dt_adot = eta_init * sqrt(a2 / j2);

    if (pot < 0 && a2 > 0)
	dt_ff = eta_init * sqrt(-pot / a2);

    dt = min(dt_adot, dt_ff);

    // Apply an upper limit to the time step:

    dt = min(dt, eta_init/4);	// eta_init/4 is arbitrary limit...

    real true_limit = step_limit;
    if (additional_step_limit > 0) true_limit = min(true_limit,
						    additional_step_limit);

    timestep = adjust_number_to_power(dt, true_limit);

    if (time != 0)
	while (fmod(time, timestep) != 0)
	    timestep *= 0.5;

    xreal tnext = time + timestep;

    if (timestep == 0 || time == tnext) {

	cerr << endl << "set_first_timestep: dt = 0 at time "
	     << get_system_time() << endl;

	PRC(dt_adot); PRL(dt_ff);
	PRC(a2); PRC(j2); PRL(pot);
	pp3(this);

	cerr << endl << endl << "System dump:" << endl << endl;

	pp3(get_root());
	exit(0);
    }
}


//-----------------------------------------------------------------------------
// new_timestep:  Calculate the next timestep following the Aarseth
//                formula (Aarseth 1985), including necessary adjustments
//                for the hierarchical timestep scheme.
//-----------------------------------------------------------------------------

static int count = 0;

local inline real new_timestep(vector& at3,		// 3rd order term
			       vector& bt2,		// 2nd order term
			       vector& jerk,		// 1st order term
			       vector& acc,		// 0th order term
			       real timestep,		// old timestep
			       real time,		// present time
			       real correction_factor,
			       hdyn * b,
			       real pert_sq)
{
    if (timestep == 0)
	cerr << "new_timestep: old timestep = 0" << endl;

    // N.B. "timestep" may really be dtau (if b->get_kappa() > 1).

    real newstep, altstep;
    bool keplstep = (pert_sq >= 0 && pert_sq < 0.0001
		     && b->is_low_level_node()
		     && b->get_kira_options()->allow_keplstep);

    real dist, dtff2, dtv2;

    bool timestep_check = (b->get_kira_diag()->timestep_check
			   && b->get_kira_diag()->check_diag(b));

    if (keplstep) {

	// Use a simple "kepler" time step criterion if b is a binary
	// and the perturbation is fairly small.

	hdyn* s = b->get_younger_sister();

	if (s) {				// excessively cautious?

	    dist = abs(b->get_pos() - s->get_pos());
	    dtff2 = dist*dist*dist / (2*b->get_parent()->get_mass());
	    dtv2  = square(b->get_pos()) / square(b->get_vel());

	    newstep =
		altstep = 0.5 * correction_factor	// 0.5 is empirical
		    	      * b->get_eta()
			      * sqrt(min(dtff2, dtv2));

	    // PRC(dist); PRC(square(b->get_vel())); PRC(dtff2); PRL(dtv2);

	} else

	    keplstep = false;

    }

    real a2, j2, k2, l2;
    real tmp1, tmp2, tmp3, tmp4;

    if (!keplstep || timestep_check) {

	// Simple criterion:
	// real newstep = b->get_eta() * sqrt(sqrt(a2 / k2));

	// Use the Aarseth criterion here.

	real dt2 = timestep * timestep;
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
#endif

	}

	newstep = aarsethstep * correction_factor;
#if 0
	PRC(timestep); PRL(newstep);
#endif
    }

    if (keplstep && timestep_check) {

	if (newstep < 1.e-13 || altstep < 1.e-13) {

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl << "in new_timestep:" << endl;

	    PRC(b->format_label()); PRC(dist);
	    PRC(b->get_posvel()); PRL(pert_sq);
	    PRC(newstep); PRC(altstep); PRL(altstep/newstep);

	    if (count > 0 || altstep/newstep > 10) {
		PRL(b->format_label());
		PRC(b->get_system_time()); xprint(b->get_system_time());
		PRC(b->get_time()); xprint(b->get_time());
		PRC(b->get_t_pred()); xprint(b->get_t_pred());
		PRL(timestep);
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
		if (altstep/newstep > 10) count++;
		if (count > 100) exit(1);
	    }
	    cerr.precision(p);
	}

	newstep = altstep;	// comment out to retain Aarseth step

    }

    // The natural time step is newstep.  Force it into an appropriate
    // power-of-two block.  A halving of timestep is always OK.  To
    // preserve the synchronization, a doubling is OK only if the
    // current time could be reached by an integral number of doubled
    // time steps (see use of fmod below).

    if (newstep < timestep)

	return 0.5 * timestep;

    else if (newstep < 2 * timestep)

	return timestep;

    else if (fmod(time, 2 * timestep * b->get_kappa()) != 0.0
	     || 2 * timestep * b->get_kappa() > b->get_step_limit())

	return timestep;

    else {

	// Added by Steve 7/28/98 to deal with pathological ICs from SPZ...
	// Do not double if this is a strongly perturbed center of mass node...

	if (b->is_leaf()
	    || b->get_oldest_daughter()->get_perturbation_squared() < 1)
	    return 2 * timestep;
	else
	    return timestep;
    }
}

//-----------------------------------------------------------------------------
// update:  Update the time and the timestep.
//-----------------------------------------------------------------------------

void hdyn::update(vector& bt2, vector& at3)    // pass arguments to
						      // avoid recomputation
{
    time += timestep;
    real dt = timestep;

    if (slow) {
	dt = slow->get_dtau();
	slow->inc_tau(dt);
    }

    // vector at3 = 2 * (old_acc - acc) + dt * (old_jerk + jerk);
    // vector bt2 = -3 * (old_acc - acc) - dt * (2 * old_jerk + jerk);

    // Reminder:	at3  =  a''' dt^3 / 6		(a' = j)
    //			bt2  =  a''  dt^2 / 2

#if 0
    if (time > 0) {
	cerr << endl;
	PRI(4); PRC(index); PRL(dt);
	PRI(4); PRL(old_acc);
	PRI(4); PRL(acc);
	PRI(4); PRL(old_jerk);
	PRI(4); PRL(jerk);
	PRI(4); PRL((acc-old_acc)/dt);
	PRI(4); PRL(2*bt2/(dt*dt));
	PRI(4); PRL((jerk-old_jerk)/dt);
	PRI(4); PRL(6*at3/(dt*dt*dt));
    }
#endif

    // Define correction factor for use in new_timestep.

    real correction_factor = 1;

    if (is_low_level_node()) {

	// Correction factor reduces the timestep in a close binary,
	// in order to improve energy errors.

	real m = mass/mbar;
	real pot_sq = m*m*d_min_sq/(pos*pos);

	// Steve 8/98:  1. 0.125 here is conservative -- expect *cumulative*
	//		   energy error to be ~fourth order.
	//		2. Could rewrite to replace pow() if necessary...
	//		   *** see kira_approx.C ***

	// Hmmm...  This pow() seems to fall victim to the strange math.h
	// bug in Red Hat Linux 5.  In that case, use Steve's approximate
	// version instead.

	// correction_factor = pow(pot_sq, -0.1);

	if (pot_sq > 1)
	    // correction_factor = pow(pot_sq, -0.125);
	    correction_factor = pow_approx(pot_sq);	// -0.125 is built in!

	// correction_factor = pow(pos*pos/d_min_sq, 0.025);

	// May be desirable to place limits on the correction...

	if (correction_factor > 1.0) correction_factor = 1.0;
	if (correction_factor < 0.2) correction_factor = 0.2;

#if 0
	if (perturbation_squared < 0.0001 && pot_sq > 1
	    && timestep < 1.e-11) {
	    PRC(format_label()); PRL(correction_factor);
	}
#endif

    }

    real new_dt = new_timestep(at3, bt2, jerk, acc, dt, time,
			       correction_factor, this, perturbation_squared);

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

    if (is_low_level_node()) {
	prev_posvel = posvel;
	posvel = pos * vel;
    }

    // Finally, store a'' for use in prediction of (top-level, for now) nodes.

    if (is_top_level_node())
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

local inline void get_derivs(vector& acc, vector& jerk,
			     vector& old_acc, vector& old_jerk,
			     real dt, vector& bt2, vector& at3)
{
    bt2 = -3 * (old_acc - acc) - dt * (2 * old_jerk + jerk);
    at3 =  2 * (old_acc - acc) + dt * (old_jerk + jerk);
}

local inline void update_derivs_from_perturbed(vector& acc_p,
					       vector& jerk_p,
					       vector& old_acc_p,
					       vector& old_jerk_p,
					       real kdt,
					       real ki2,	// not used
					       real ki3,	// not used
					       vector& acc,
					       vector& jerk,
					       vector& bt2,
					       vector& at3,
					       bool print = false)
{
    // Return here to neglect all perturbed terms and compute slow binary
    // CM motion in the center of mass approximation.

    // if (CM_ONLY) return;

    vector bt2_p, at3_p;
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
			       vector &acc, vector &jerk,
			       vector &bt2, vector &at3,
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

    vector acc_p = s->get_acc_p();
    vector jerk_p = s->get_jerk_p();

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

    vector old_acc_p = s->get_old_acc_p();
    vector old_jerk_p = s->get_old_jerk_p();

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

    real dt = timestep;
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


	    vector a_2b, j_2b;
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

	    vector dx = pred_pos - get_younger_sister()->pred_pos;
	    real r2 = dx*dx;
	    vector a = -get_younger_sister()->mass * dx / pow(r2, 1.5);

	    PRL(a);

	    // pp3(get_parent());

	    sb->slow->set_acc_p(acc-a_2b);
	    sb->slow->set_jerk_p(jerk-j_2b);

	    acc = a_2b;
	    jerk = j_2b;
	}
    }

    vector bt2, at3;
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

		    vector acc_p = s->get_acc_p();	// as above, better to
		    vector jerk_p = s->get_jerk_p();	// have _dyn_ as friend?

		    // Scale acc_p and jerk_p (see comment above).

		    acc_p *= ki2;
		    jerk_p *= ki3;

		    // Propogate the changes back to the class.

		    s->set_acc_p(acc_p);
		    s->set_jerk_p(jerk_p);

		    if (corr) {

			// Apply the correction if the old_ quantities are set.

			vector old_acc_p = s->get_old_acc_p();

			if (square(old_acc_p) > 0) {

			    vector old_jerk_p = s->get_old_jerk_p();
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

    vector new_pos = pred_pos + (0.05 * at3 + ONE12 * bt2) * dt * dt;
    vector new_vel = pred_vel + (0.25 * at3 + ONE3 * bt2) * dt;

#if defined(STARLAB_HAS_GRAPE4) || defined(STARLAB_HAS_GRAPE4)

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

#if defined(STARLAB_HAS_GRAPE4)
		PRL(get_grape_chip(this));	// direct access to data
						// in hdyn_grape4.C
#endif

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
	    sprintf(tmp, "runaway in correct at time %f", time);
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

#endif

    pos = new_pos;
    vel = new_vel;

    //--------------------------------------------------------
    // Call update here to avoid recomputation of bt2 and at3.
    //--------------------------------------------------------

    update(bt2, at3);

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
				    vector& d_pos,
				    vector& d_vel,
				    real eps2,
				    vector& a,
				    vector& j,
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

void hdyn::flat_calculate_acc_and_jerk(hdyn * b,    	// root node
				       bool make_perturber_list)
{
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "    flat_calculate_acc_and_jerk for "
	     << format_label() << endl;
    }

    acc = jerk = 0;
    pot = 0;

    // Initialize both nn and coll:

    nn = NULL;
    d_nn_sq = VERY_LARGE_NUMBER;
    coll = NULL;
    d_coll_sq = VERY_LARGE_NUMBER;

    real distance_squared;

    for_all_daughters(hdyn, b, bi)	
	if (bi != this) {

	    vector d_pos = bi->get_pred_pos() - pred_pos;
	    vector d_vel = bi->get_pred_vel() - pred_vel;
	    accumulate_acc_and_jerk(bi,				// (inlined)
				    d_pos, d_vel,
				    eps2, acc, jerk, pot, distance_squared);

	    // Note that the nn and coll passed to update_nn_coll here
	    // are the actual pointers, so this update changes the
	    // data in b.

	    real sum_of_radii = get_sum_of_radii(this, bi);
	    update_nn_coll(this, 1,		// (1 = ID)	// (inlined)
			   distance_squared, bi, d_nn_sq, nn,
			   sum_of_radii, d_coll_sq, coll);

	    // Note: we do *not* take the slowdown factor into account
	    //	     when computing the perturber list.

	    // See equivalent code for use with GRAPE-4 in
	    // hdyn_grape4.C/get_neighbors_and_adjust_h2.

	    if (make_perturber_list
		&& is_perturber(this, bi->mass,			// (inlined)
				distance_squared,
				perturbation_radius_factor)) {
		if (n_perturbers < MAX_PERTURBERS) {
		    perturber_list[n_perturbers] = bi;
#if 0
		    cerr << "added " << bi->format_label();
		    cerr << " to perturber list of "
			 << format_label()
			 << endl;
#endif
		}
		n_perturbers++;
	    }
	}
}

void hdyn::perturber_acc_and_jerk_on_leaf(vector &a,
					  vector &j,
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

    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "        perturber_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
    }

    dbg_message("perturber_acc_and_jerk_on_leaf", this);

    a = j = 0.0;
    p = 0;

    // Do *not* set p_d_nn_sq = VERY_LARGE_NUMBER here!
    //						(Steve, 11/00)

    d_coll_sq = VERY_LARGE_NUMBER;

    if (!pnode->valid_perturbers || pnode->n_perturbers <= 0)
	return;

    // Added &hdyn:: to four "hdyn_something_relative_to_root" (Steve, 6/27/97).

    // Removed all four "hdyn_something_relative_to_root" references, which are
    // quite inefficient (Steve, 8/20/98).

    // Determine absolute position and velocity of 'this', using explicit code.

    vector d_pos = -pred_pos;
    vector d_vel = -pred_vel;

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

    // vector d_pos = -hdyn_something_relative_to_root(this,
    // 						    &hdyn::get_pred_pos);
    // vector d_vel = -hdyn_something_relative_to_root(this,
    // 						    &hdyn::get_pred_vel);

    // Loop over the perturber list.

    vector d_pos_p;
    vector d_vel_p;

    bool reset = false;

    for (int i = 0; i < pnode->n_perturbers; i++) {

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

    step_node->inc_indirect_force(pnode->n_perturbers);		// bookkeeping
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
						      vector& offset_pos,
						      vector& offset_vel,
						      vector& a,
						      vector& j,
						      real& p,
						      real& p_d_nn_sq,
						      hdyn * &p_nn,
						      bool point_mass_flag,
						      hdyn *step_node)
{
    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.

    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "        tree_walk_for_partial_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
	cerr << "        b = " << b->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }

    dbg_message("tree_walk_for_partial_acc_and_jerk_on_leaf", this, b);

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

	    vector ppos = offset_pos + d->get_pred_pos();
	    vector pvel = offset_vel + d->get_pred_vel();

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
						  vector& a,
						  vector& j,
						  real & p,
						  real & p_d_nn_sq,
						  hdyn * &p_nn,
						  bool point_mass_flag,
						  hdyn * pnode,
						  hdyn * step_node)
{
    // The input arguments p_d_nn_sq and p_nn may be the actual
    // d_nn_sq and nn, or copies.

    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "      calculate_partial_acc_and_jerk_on_leaf for "
	     << format_label() << endl;
	cerr << "      top = " << top->format_label();
	cerr << ", common = " << common->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }

    dbg_message("calculate_partial_acc_and_jerk_on_leaf", this);
    dbg_message("                 from particle", top);

    a = j = 0.0;
    p = 0;

    // Determine absolute position and velocity of "this":

    vector d_pos = 0;
    vector d_vel = 0;

    // Loop over the appropriate perturber list for external perturbations,
    // and traverse the tree below top, masked by mask, for internal
    // perturbations due to other members of the parent clump.

    if (pnode)
	perturber_acc_and_jerk_on_leaf(a, j, p,
				       p_d_nn_sq, p_nn,
				       pnode, step_node);

    if (top == mask) return;

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
					  vector& a,
					  vector& j,
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

    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "    calculate_partial_acc_and_jerk for "
	     << format_label() << endl;
	cerr << "    top = " << top->format_label();
	cerr << ", common = " << common->format_label();
	cerr << ", mask = " << mask->format_label() << endl;
    }

    dbg_message("calculate_partial_acc_and_jerk", this);

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

	vector a_daughter, a_save;
	vector j_daughter;
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
	// top-level node whose neighbor list has overflowed (or which
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

	    // "Daughter" quantities refer to the younger daughter on
	    // leaving the loop.  The vector a_2b is the two-body
	    // relative acceleration of the components (computed here
	    // in the point-mass approximation).

	    hdyn* od = get_oldest_daughter();
	    hdyn* yd = get_oldest_daughter()->get_younger_sister();

	    // Perturbation_squared does *not* contain the slowdown factor.

	    real pscale = m_daughter / (m_daughter + mass);
	    real r2 = square(od->pos - yd->pos);
	    od->perturbation_squared =
			square(pscale * (a_save - a_daughter))
			        * square(r2/mass);
 	}
    }
}


// check_add_perturber:  check if p is a perturber of 'this' and add it
//			 to the perturber list if necessary.

void hdyn::check_add_perturber(hdyn* p, vector& this_pos)
{
    vector ppos = hdyn_something_relative_to_root(p, &hdyn::get_pred_pos);

//     bool print = get_top_level_node()->n_leaves() >= 4;
//     if (print) {
// 	cerr << "check_add_perturber:  checking " << p->format_label()
// 	     << " at time " << system_time;
// 	cerr << " for node " << format_label() << endl;
// 	PRL(mass);
// 	PRL(square(this_pos - ppos));
// 	PRL(perturbation_radius_factor);
//    }

    // Note: we do *not* take the slowdown factor into account
    //	     when computing the perturber list.

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
    } //else
//	if (print) cerr << "...rejected" << endl;
}

void hdyn::create_low_level_perturber_list(hdyn* pnode)
{
    valid_perturbers = true;
    if (perturber_list == NULL)
	perturber_list = new hdyn *[MAX_PERTURBERS];

    perturbation_radius_factor
	= define_perturbation_radius_factor(this, gamma23);

    vector this_pos = hdyn_something_relative_to_root(this,
						      &hdyn::get_pred_pos);
    n_perturbers = 0;

    // First add sisters, aunts, etc. up to pnode...

    hdyn* p = this;
    while (p != pnode) {
	check_add_perturber(p->get_binary_sister(), this_pos);
	p = p->get_parent();
    }

    // ...then accept any node on the pnode list that satisfies the
    // inner node perturbation criterion.

    for (int i = 0; i < pnode->n_perturbers; i++)
	check_add_perturber(pnode->perturber_list[i], this_pos);

    if (n_perturbers > MAX_PERTURBERS) {

	valid_perturbers = false;
	delete [] perturber_list;

    } else {

// 	if (get_top_level_node()->n_leaves() >= 4) {
// 	    cerr << ">>>> this->"; PRL(format_label());
// 	    cerr << ">>>> "; PRL(pnode->format_label());
// 	    cerr << ">>>> "; PRL(pnode->n_perturbers);
// 	    cerr << ">>>> this->"; PRL(n_perturbers);
// 	    print_perturber_list();
// 	}

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
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  calculate_acc_and_jerk_on_low_level_node for "
	     << format_label() << endl;
    }

    dbg_message("calculate_acc_and_jerk_on_low_level_node", this);

    hdyn *root = get_root();

    if (parent->get_oldest_daughter()->get_younger_sister()
				     ->get_younger_sister() != NULL)
	err_exit("calculate_acc_and_jerk_on_low_level_node: Not binary node");

    hdyn *sister = get_binary_sister();

    real mtot = 0;

    vector apert1, jpert1;
    vector apert2, jpert2;
    real p_dummy;
    hdyn *top;

    // New formulation of perturber lists introduced by Steve 8/98:

    // Perturber list may not necessarily be associated with the
    // top-level node.  Any CM node above this node may have a
    // valid perturber list.  Use the lowest-level one.

    // Now pass a pointer to the node containing the perturber list
    // (NULL ==> no valid list) instead of bool "use_perturber_list"
    // flag.  Also, since the low-level perturber list includes all
    // members of the parent "clump" except those below the parent
    // node, we need only descend the tree below pnode to pick up the
    // remaining perturbations.

    hdyn* top_level = get_top_level_node();
    hdyn* pnode = find_perturber_node();

    if (pnode) {
	top = pnode;
	kc->pert_step++;
	kc->pert_with_list += pnode->n_perturbers;
    } else {
	top = root;
	kc->pert_without_list++;
    }

    // Acceleration and jerk on this component due to rest of system:

    nn = NULL;
    d_nn_sq = VERY_LARGE_NUMBER;

    calculate_partial_acc_and_jerk(top, top, get_parent(),
				   apert1, jpert1, p_dummy,
				   d_nn_sq, nn,
				   !USE_POINT_MASS,		// explicit loop
				   pnode,			// with list
				   this);			// node to charge

    // Acceleration and jerk on other component due to rest of system:

    sister->set_nn(NULL);
    sister->set_d_nn_sq(VERY_LARGE_NUMBER);

    sister->calculate_partial_acc_and_jerk(top, top, get_parent(),
					   apert2, jpert2, p_dummy,
					   sister->d_nn_sq, sister->nn,
					   !USE_POINT_MASS,	// explicit loop
					   pnode,		// with list
					   this);		// node to charge

    // Note:  The first two calls to calculate_partial_acc_and_jerk pass
    // the d_nn_sq and nn from the hdyn, so the hdyn data are actually
    // updated by update_nn_coll.  The third call (below) passes local
    // variables, so it has no direct effect on d_nn_sq and nn.  (Test
    // afterwards if d_nn_sq and nn must be updated.)  The coll data
    // are always updated by these calls.

    // Relative acceleration and jerk due to other (sister) component:

    vector a_2b, j_2b;
    real p_2b;

    hdyn *p_nn_sister;
    real d_min_sister = VERY_LARGE_NUMBER;
    calculate_partial_acc_and_jerk(get_parent(), get_parent(), this,
				   a_2b, j_2b, p_2b,
				   d_min_sister, p_nn_sister,
				   !USE_POINT_MASS,
				   NULL,			// no perturbers
				   this);			// node to charge

    real m2 = sister->get_mass();
    real pscale = m2 / (m2 + mass);

    set_acc_and_jerk_and_pot(a_2b + pscale * (apert1 - apert2) * get_kappa(),
			     j_2b + pscale * (jpert1 - jpert2), p_2b);

    // Note that perturbation_squared does *not* contain the slowdown factor.

    perturbation_squared = square(pscale * (apert1 - apert2)) / square(a_2b);

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
	vector sep = pred_pos * (1 + mass/m2);
	real r2 = sep*sep;
	vector a2 = -m2*sep / (r2*sqrt(r2));
	PRL(a2);

	PRL(pscale * (apert1 - apert2) * get_kappa());
    }
#endif

    if (d_nn_sq > d_min_sister) {

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

    if (ALLOW_LOW_LEVEL_PERTURBERS && pnode) {

	// Create/revise the low-level perturber list(s).

	// Could revise list as we compute the force, but the distances used
	// should be from the center of mass, not from one of the components.
	// This may be a little less efficient, but it is much easier to code,
	// so keep it for now.  (Steve, 8/98)

	if (is_parent())
	    create_low_level_perturber_list(pnode);

	if (get_binary_sister()->is_parent())
	    get_binary_sister()->
		create_low_level_perturber_list(pnode);
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
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  calculate_acc_and_jerk_on_top_level_node for "
	     << format_label() << endl;
    }

    top_level_node_prologue_for_force_calculation(exact);

    if (!exact) {
	top_level_node_real_force_calculation();
	top_level_node_epilogue_force_calculation();
    }

}


void hdyn::top_level_node_prologue_for_force_calculation(bool exact)
{
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  top_level_node_prologue_for_force_calculation for "
	     << format_label() << endl;
    }

    hdyn *root = get_root();
    d_coll_sq = VERY_LARGE_NUMBER;
    coll = NULL;

    if (exact) {

	// Perform the entire force calculation.

	n_perturbers = 0;
	valid_perturbers = false;

	nn = NULL;
	d_nn_sq = VERY_LARGE_NUMBER;

	calculate_partial_acc_and_jerk(root, root, this,
				       acc, jerk, pot, d_nn_sq, nn,
				       !USE_POINT_MASS,
				       NULL,		// no perturbers
				       this);		// node to charge
    } else if (is_parent()) {

	// Set up computation of perturber list.

	if (perturber_list == NULL)
	    perturber_list = new hdyn *[MAX_PERTURBERS];

	perturbation_radius_factor
		= define_perturbation_radius_factor(this, gamma23);
    }
}

void hdyn::top_level_node_real_force_calculation()
{
    // ***************************************************************
    // **** This function is NEVER called if GRAPE is available   ****
    // **** -- it is replaced by grape4/6_calculate_acc_and_jerk. ****
    // ***************************************************************

    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  top_level_node_real_force_calculation for "
	     << format_label() << endl;
    }

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

        n_perturbers = 0;
        valid_perturbers = true;
    }

    flat_calculate_acc_and_jerk(root, is_parent());

    if (nn == NULL) {
        pretty_print_node(cerr); cerr << " nn NULL after flat " << endl;
    }
}

void hdyn::top_level_node_epilogue_force_calculation()
{
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "  top_level_node_epilogue_force_calculation for "
	     << format_label() << endl;
    }

#if 0
    if (time > 13.62265) {
	cerr << "top_level_node_epilogue_force_calculation(1): " << endl;
	PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
    }
#endif

    // Should NOT be called in exact case.

    // Take care of the case where the nearest neighbor is
    // a complex node, since otherwise nn may disappear.
    // (Jan 21, 1998, JM and SPZ)

    // Exclude the case of nn = this (possible when GRAPE is used.)
    // (Sep 9, 1998, SLWM)

    if (nn != this) {
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
	cerr << "top_level_node_epilogue_force_calculation(2): " << endl;
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

    if (n_perturbers > MAX_PERTURBERS) {
	
	// Perturber list has overflowed.  Use the entire tree (below).

	valid_perturbers = false;
	kc->perturber_overflow++;
	
    } else if (n_perturbers > 0) {
	
	// Use the perturber list to correct the center-of-mass force.
	// First, subtract out center-of-mass contribution:
	
	vector a_cm, a_p, j_cm, j_p;
	real p_p;
	calculate_partial_acc_and_jerk(this, this, this,
				       a_cm, j_cm, p_p, d_nn_sq, nn,
				       USE_POINT_MASS,	    	    // explicit
				       this,	   		    // loop
				       this);

#if 0
	if (time > 13.62265) {
	    cerr << "top_level_node_epilogue_force_calculation(3a): " << endl;
	    PRC(format_label()); PRC(nn->format_label()); PRL(sqrt(d_nn_sq));
	}
#endif

	// (Note: 'this' is OK here because this is a top-level node.)

	// acc -= a_cm;
	// jerk -= j_cm;

	pot -= p_p;

	// Replace nodes by components on the perturber list and
	// add in contribution due to components.

	bool debug = false;

	if (!expand_nodes(n_perturbers, perturber_list, debug)) {
	
	    // Perturber list has overflowed.  Use the entire tree (below).
	    // Note: may cause problems if this is a slow binary...
	
	    valid_perturbers = false;
	    kc->perturber_overflow++;
	
	} else {

	    // (Explicit loop)
	
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




// 	    if (time >= xreal(2296, 3651856000000000000)) {
// 		cerr << endl << "in top_level_node_epilogue_force_calculation"
// 		     << endl;
// 		int p = cerr.precision(HIGH_PRECISION);
// 		PRL(format_label());
// 		pp3(get_top_level_node());
// 		cerr << endl;
//		cerr.precision(p);
// 	    }




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
	
	}

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

    // Recompute the force exactly if the perturber list has overflowed.

    if (!valid_perturbers) {
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

}


//-----------------------------------------------------------------------------
// calculate_acc_and_jerk:  Generic routine to calculate the acceleration
// of a node.  Call appropriate functions depending on whether or not the
// node is at the top level of the tree.
//-----------------------------------------------------------------------------

void hdyn::calculate_acc_and_jerk(bool exact)
{
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "calculate_acc_and_jerk for "
	     << format_label() << endl;
    }

    dbg_message("calculate_acc_and_jerk", this);

    d_coll_sq = VERY_LARGE_NUMBER;
    coll = NULL;
    if (is_low_level_node())
	get_binary_sister()->set_d_coll_sq(VERY_LARGE_NUMBER);

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
         pretty_print_node(cerr); cerr << " nn NULL after calc " << endl;
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
    vector a_c = 0, j_c = 0;
    real p_c = 0;
    vector a_p = 0, j_p = 0;
    real p_p = 0;
    real dum_d_sq = VERY_LARGE_NUMBER;
    hdyn *dum_nn = NULL;

    if (bj->get_kira_diag()->correct) {
	cerr << "correcting " << bi->format_label();
	cerr << " for " << bj->format_label() << " at "
	     << bi->get_system_time() << endl;
    }

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
	// ...or if neighbor-list overflow has occurred (Steve 7/98).

	if (btop->get_valid_perturbers()) {
	    for (int j = 0; j < btop->get_n_perturbers(); j++)
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
    //	    1. btop is on the integration list (time = system time)
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
		// calculation was exact, or if neighbor-list overflow
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
//  integrate_node -- advance this node by one time step
//-----------------------------------------------------------------------------

void hdyn::integrate_node(hdyn * root,
			  bool integrate_unperturbed,	// default = true
			  bool force_unperturbed_time)	// default = false

// From Steve: force_unperturbed_time can cause an unperturbed binary
// to be advanced to an undesirable phase.

{
    if (diag->ev_function_id && diag->check_diag(this)) {
	cerr << "integrate_node for " << format_label()
	     <<" at time " << time  + timestep << endl;
	// pretty_print_node(cerr);
    }

    if (kep == NULL) {

	clear_interaction();
	calculate_acc_and_jerk(true);
	set_valid_perturbers(false);

	if (tidal_type && is_top_level_node()) {
	    real dpot;
	    vector dacc, djerk;
	    get_external_acc(this, pred_pos, pred_vel,
			     dpot, dacc, djerk);
	    inc_pot(dpot);
	    inc_acc(dacc);
	    inc_jerk(djerk);
	}

	correct_and_update();			// note: no retry on error
	// update();

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

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  src/dyn/evolve/hdyn_ev.C
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
