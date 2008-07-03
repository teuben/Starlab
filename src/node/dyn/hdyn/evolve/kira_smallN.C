
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                                     _\|/_
//=======================================================//              /|\ ~


//// Self-contained few-body integrator, using a fourth-order Hermite
//// scheme incorporating a modified unperturbed treatment of close
//// approaches and time symmetrization.  The program reads a snapshot
//// from standard input and optionally writes snapshot(s) to standard
//// output.  Optional periodic log output is sent to standard error.
//// This version runs as a standalone program, or incorporated into kira
//// to handle multiple systems containing perturbed hard binaries.
////
//// Usage: kira_smallN [OPTIONS] < input > ouptut
////
//// Options:
////         -a    set accuracy parameter [0.03]
////         -d    set log output interval [0 --> no output]
////         -D    set snap output interval [0 --> no output]
////         -E    set energy output interval [0 --> no output]
////         -g    set unperturbed limit [1.e-5]
////         -n    set number of symmetrization iterations [1]
////         -r    set termination radius [infinite]
////         -t    set termination time [200]
////
//// Authors: Fan-Chi Lin and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// Externally visible library functions:
//
//	void set_smallN_eta()
//	void set_smallN_gamma()
//	void set_smallN_niter()
//	void set_smallN_dtcrit()
//	void set_smallN_rcrit()
//
//	real get_smallN_eta()
//	real get_smallN_gamma()
//	int  get_smallN_niter()
//	real get_smallN_dtcrit()
//	real get_smallN_rcrit()
//
//	real smallN_evolve()
//	real integrate_multiple()
//
// As currently envisaged, on entry from kira we will have an isolated
// (that is, relatively unperturbed) few-body (N <~ 5) system in which
// a hard binary has just become perturbed by a neighbor.  A typical
// scenario might be that a stable triple has been perturbed into an
// unstable orbit, and now the outer component has become a perturber
// of the inner pair.  Fow now, we
//
//	* integrate the multiple in isolation,
//	* stop when an unperturbed hard binary is again identified,
//	  using the same criteria as used in kira,
//	* reinsert the multiple into the main integration after
//	  zero elapased system time,
//	* keep track of any resulting tidal error.
//
// Later, we will
//
//	* include external perturbations, probably using a multipole
//	  approximation rather than summing over a perturber list,
//	* treat the multiple as an object in the main kira integration,
//	  with its own time step and internal structure that is visible
//	  to the rest of the system.  Since predictions always go in the
//	  same (forward) direction in time, predicting the internal
//	  structure will simply cause the multiple to be integrated to
//	  the specified time.  Correction will have no additional effect.
//	  Need to check that the normally flat but not binary tree used
//	  here will work as is with the hdyn_ev functions, or if additional
//	  functions need be written.
//	* largely elimate the tidal error by these actions.
//
// The internal interaction currently uses scaled times and coordinates
// to avoid precision problems with small time steps and to facilitate
// visualization.  However, this may not actually be necessary so long as
// the system time is always offset so that all times are on the order of
// the time steps used.  Need to experiment with this.  For compatibility
// with the rest of kira, do not use system_time anywhere here, except
// temporarily modified for output purposes.
//
//						Fan-Chi Lin
//						Steve McMillan (10/04)

#include  "hdyn.h"
#include  "hdyn_inline.C"
#include "kira_defaults.h"



#ifndef TOOLBOX

// Global stuff...

// *** Too many global variables! ***

// Integration parameters are global to the library function, settable
// from outside via accessor functions.

static real eta = 0.03;			// empirical
static real gamma2_unpert = 1.e-10;	// conservative, given first-order
					// correction now employed

// Number of iterations to symmetrize(0 == > explicit, > 0 ==> implicit).
// Default of 1 seems to work well.

static int  n_iter = 1;

// These parameters simply limit the test for unperturbed motion.
// Default values are valid only for standard units.

static real dt_crit = 1.e-4;
static real r2_crit = 1.e-6;

static real smallN_energy_scale = 0;
static real smallN_time_scale = 0;

void set_smallN_eta(real a)	{eta = a;}
void set_smallN_gamma(real g)	{gamma2_unpert = g*g;}
void set_smallN_niter(int n)	{n_iter = n;}
void set_smallN_dtcrit(real dt)	{dt_crit = dt;}
void set_smallN_rcrit(real r)	{r2_crit = r*r;}

real get_smallN_eta()		{return eta;}
real get_smallN_gamma()		{return sqrt(gamma2_unpert);}
int  get_smallN_niter()		{return n_iter;}
real get_smallN_dtcrit()	{return dt_crit;}
real get_smallN_rcrit()		{return sqrt(r2_crit);}

// Debugging:

static bool dbg_allow_binary = true;	// set false to disable binaries
static real rbin = 0;			// use to keep track of close pair

const bool VERBOSE = false;

// Global pointers to the closest pair and next nearest object.

static hdyn *bi_min = NULL, *bj_min = NULL, *bk_nn = NULL;
static real min_distance2 = VERY_LARGE_NUMBER, rij_2 = VERY_LARGE_NUMBER;



local inline real kepler_step_sq(real distance2, real distance3,
				 real mass, real vel2)
{
    real dtff2 = 0.5 * distance3 / mass;

    // When this function is called from integrate_multiple to set
    // scales, the velocity may be unrealistic, or even negative.
    // Small velocities are OK, but make sure negative ones don't
    // mess us up...

    real dtv2 = distance2 / abs(vel2);

    return Starlab::min(dtff2, dtv2);
}

local real kepler_step_sq(real distance, real mass, real vel2)
{
    return kepler_step_sq(pow(distance,2), pow(distance,3), mass, vel2);
}

local inline real get_pairwise_acc_and_jerk(hdyn *bi, hdyn *bj,
					    vec &force, vec &jerk)
{
    // Compute the force and jerk on node bi due to node bj by summing
    // all pairwise component forces.  Return the unscaled time step
    // appropriate to the minimum distance between any component of bi
    // and any component of bj.

    force = jerk = 0;
    real min_distance3;
    min_distance2 = VERY_LARGE_NUMBER;

    vec delx = 0, delv = 0;
    if (bi->is_parent()) {		// NOTE: assuming for now that our trees
	delx -= bi->get_nopred_pos();	//	 are just a single level deep
	delv -= bi->get_nopred_vel();
    }
    if (bj->is_parent()) {
	delx += bj->get_nopred_pos();
	delv += bj->get_nopred_vel();
    }

    for_all_leaves(hdyn, bi, bbi) {

	// Compute the force and jerk on leaf bbi due to all leaves under bj.

	vec iforce = 0, ijerk = 0;
	for_all_leaves(hdyn, bj, bbj) {

	    // Note use of nopred to return current value of pred
	    // quantities without prediction or setting t_pred.

	    vec dx = bbj->get_nopred_pos() - bbi->get_nopred_pos() + delx;
	    vec dv = bbj->get_nopred_vel() - bbi->get_nopred_vel() + delv;
	    real distance2 = dx*dx;
	    real distance  = sqrt(distance2);
	    real distance3 = distance2*distance;

	    iforce += bbj->get_mass() * dx / distance3;
	    ijerk  += bbj->get_mass() * (dv / distance3
					 - 3*dx*(dx*dv)/(distance3*distance2));

	    if (distance2 < min_distance2) {
		min_distance2 = distance2;
		min_distance3 = distance3;
	    }
	}

	force += bbi->get_mass() * iforce;
	jerk  += bbi->get_mass() * ijerk;
    }

    // Time step criterion is exactly the same as used in kira::kepler_step()
    // (see hdyn_ev).C, except it is squared and lacks a correction factor.

    real timestep2 = kepler_step_sq(min_distance2, min_distance3,
				    bj->get_mass() + bi->get_mass(),
				    square(bj->get_nopred_vel()
					     - bi->get_nopred_vel()));

    // Alternate squared time step criterion:

    // timestep2 = min_distance3 / (bj->get_mass() + bi->get_mass());

    return timestep2;
}



local inline real calculate_top_level_acc_and_jerk(hdyn *b)
{
    // Compute the acc and jerk on all top-level nodes.  We could use the
    // hdyn functions in hdyn_ev.C to do this more simply, but they will
    // double-count all pairwise interactions.  All nodes must be resolved
    // into components for purposes of computing the acc and jerk.  Return
    // the minimum time step associated with any top-level pair.

    for_all_daughters(hdyn, b, bi) {
	bi->set_acc(0);
	bi->set_jerk(0);
    }

    real min_timestep2 = VERY_LARGE_NUMBER;
    bi_min = bj_min = NULL;

    for_all_daughters(hdyn, b, bi) {
	for (hdyn *bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {

	    vec force, jerk;
	    real timestep2 = get_pairwise_acc_and_jerk(bi, bj, force, jerk);

#if 0
	    cerr << bi->format_label() << " ";
	    cerr << bj->format_label() << " ";
	    cerr << abs(acc) << " " << abs(jerk) << endl;
#endif

	    real mi = 1 / bi->get_mass();
	    bi->inc_acc(mi*force);
	    bi->inc_jerk(mi*jerk);

	    real mj = 1 / bj->get_mass();
	    bj->inc_acc(-mj*force);
	    bj->inc_jerk(-mj*jerk);

	    if (timestep2 < min_timestep2) {
		min_timestep2 = timestep2;
		bi_min = bi;
		bj_min = bj;
		rij_2 = min_distance2;
	    }
	}
    }

    real dt = eta*sqrt(min_timestep2);	// natural time step

#if 0

    // Experiment: force dt to a power of 2.

    real dt2 = 1;
    while (dt2 > dt) dt2 /= 2;
    dt = dt2;

#endif

    return dt;
}



//----------------------------------------------------------------------
//
// Advance the components of a binary to the specified time.

local inline void advance_components_to_time(hdyn *bi, real t)	  // unperturbed
{								  // "predictor"
    hdyn *od = bi->get_oldest_daughter();

    if (od) {	// redundant

	// Only way this can occur is for unperturbed motion.
	// Check and flag if no kepler found.

	kepler *k = od->get_kepler();

	if (!k)
	    warning("smallN_evolve: daughter node with no kepler.");

	else {

	    // Advance the components to time t.  We will not
	    // integrate the internal motion, so could just set
	    // pos = pred_pos (etc.) here (no corrector), but
	    // for now defer that until the corrector step.

	    k->transform_to_time(t);

	    real fac = od->get_mass()/bi->get_mass();
	    od->set_pred_pos(-fac*k->get_rel_pos());
	    od->set_pred_vel(-fac*k->get_rel_vel());

	    hdyn *yd = od->get_younger_sister();
	    yd->set_pred_pos((1-fac)*k->get_rel_pos());
	    yd->set_pred_vel((1-fac)*k->get_rel_vel());
	}
    }
}

local inline void update_components_from_pred(hdyn *bi)		  // unperturbed
{								  // "corrector"
    // Called after a step is completed.

    hdyn *od = bi->get_oldest_daughter();

    if (od) {	// redundant

	od->set_time(bi->get_time());
	od->set_pos(od->get_nopred_pos());
	od->set_vel(od->get_nopred_vel());
	hdyn *yd = od->get_younger_sister();
	yd->set_time(bi->get_time());
	yd->set_pos(yd->get_nopred_pos());
	yd->set_vel(yd->get_nopred_vel());
    }
}

local inline void correct_timestep_for_components(hdyn *bi, real t, real &dt)
{
    // Reduce the time step if necessary to end at the termination
    // time for unperturbed motion.

    hdyn *od = bi->get_oldest_daughter();

    if (od) {	// redundant

	// Termination time is stored in od->t_pred.

	real old_dt = dt;
	real t_term = od->get_t_pred();
	if (t + dt > t_term) dt = t_term - t;	// should check for dt small

	if (dt != old_dt) {
#if 0
	    cerr << "correct_dt for " << bi->format_label() << ":" << endl;
	    PRC(t); PRC(t_term); PRC(old_dt); PRL(dt);
#endif
	}
    }
}

local inline bool unperturbed_and_approaching(hdyn *b, hdyn *bi, hdyn *bj)
{
    // Settle for close and approaching for now.  Eventually
    // require unperturbed motion, but add that later.  Also only
    // allow one binary at a time, until the perturbation criterion
    // is in place.  OK to use pos and vel here because pred_pos
    // and pos are the same at the end of a step, which is when
    // this function is called.
 
    if (bi->is_parent() || bj->is_parent()) return false;
    if (square(bi->get_pos()-bj->get_pos()) > r2_crit) return false;
    if ((bi->get_pos()-bj->get_pos()) * (bi->get_vel()-bj->get_vel()) > 0)
      return false;

    // Particles are within r_crit and approaching.
    // Check the perturbation due to the rest of the system.

    bk_nn = NULL;
    real dis2, near_dis2 = VERY_LARGE_NUMBER;
    real min_dis2;
    for_all_leaves(hdyn, b, bk)
      {
	if( (bk != bi) && (bk != bj) )
	  { 
	    dis2 = square(bk->get_pos()-bi->get_pos());
	    if (dis2 < near_dis2)
	      {
		near_dis2 = dis2;
		bk_nn = bk;
	      }
	  }
      }

    // Particle bk_nn is the next nearest neighbor to bi.  (Criterion
    // is based on distance -- should eventually use potential to
    // include masses.)  Compute the perturbation.

    min_dis2 = square(bi->get_pos()-bj->get_pos());
    real mass = bi->get_mass() + bj->get_mass();
 
    if (gamma2_unpert / 0.25 <
	(pow(min_dis2, 3)/square(mass)) / 
	(pow(near_dis2, 3)/square(mass+bk_nn->get_mass())))
      return false;

    // Perturbation is small.  Estimate the time scale for bk_nn
    // to approach (bi,bj).

#if 0

    // Not completely clear what this is doing -- Steve.

    if (gamma2_unpert < 
	square(min_dis2/square(bi->get_vel()-bj->get_vel())/
	       near_dis2*square((bi->get_vel()*bi->get_mass()+
	       bj->get_vel()*bj->get_mass())/mass-bk_nn->get_vel())
	       )
       ) return false;
#endif

    return true;
}



local inline bool create_binary(hdyn *bi, hdyn *bj)
{
    // Combine two nodes, replacing the first by the center of mass of
    // the two, and moving both to lie below the new center of mass node.
    // On entry we have just completed a step, so pos and pred_pos should
    // be the same (vel similarly).

    if (!bi || !bj) return false;
    if (bi->get_pot() != 0) return false;

    // Construct a kepler structure describing the relative motion.
    // Code below follows dyn_to_kepler() in dyn/kepler/kepler.C.

    real time = bi->get_time();		// all times should be the same
    real mi = bi->get_mass();
    real mj = bj->get_mass();

    // As above, OK to use pos and vel here because pred_pos
    // and pos are the same at the end of a step, which is when
    // this function is called.

    kepler *k = new kepler;
    k->set_time(time);
    k->set_total_mass(mi+mj);
    k->set_rel_pos(bj->get_pos() - bi->get_pos() );
    k->set_rel_vel(bj->get_vel() - bi->get_vel() );
    k->initialize_from_pos_and_vel();

    // Set a time step.  Choose time to separation at the same radius.
    // Could use pred_advance_to_periastron(), but this may be quite
    // inaccurate.  Better to use mean motion.  (Steve, 4/99)

    real mean_anomaly = k->get_mean_anomaly();
    if (k->get_energy() < 0)
	mean_anomaly = sym_angle(mean_anomaly);	    // place in (-PI, PI]

    real peri_time;
    if (mean_anomaly >= 0) {	// should never happen:
	delete k;
	return false;
    } else
	peri_time = time - mean_anomaly / k->get_mean_motion();

    // Create a new center of mass node.

    hdyn *cm = new hdyn;
    cm->set_mass(mi+mj);
    cm->set_pos((mi*bi->get_pos()+mj*bj->get_pos())/cm->get_mass());
    cm->set_vel((mi*bi->get_vel()+mj*bj->get_vel())/cm->get_mass());
    cm->set_time(time);

    // Compute the tidal potential (components - CM) due to the rest
    // of the system.  Store it in bi->pot.

    real dphi_tidal = 0;
    for_all_daughters(hdyn, bi->get_parent(), bk)
	if (bk != bi && bk != bj) {
	    real dphi_k = -bi->get_mass()/abs(bi->get_pos()-bk->get_pos())
			  -bj->get_mass()/abs(bj->get_pos()-bk->get_pos())
			  +cm->get_mass()/abs(cm->get_pos()-bk->get_pos());
	    dphi_tidal += bk->get_mass()*dphi_k;
	}
    bi->set_pot(dphi_tidal);

    //-----------------------------------------------------------------
    // Debugging:

    if (!dbg_allow_binary) {

	// Don't actually create the binary.  Print info and return.

	bj->set_pot(dphi_tidal);
	PRL(bi->format_label());
	PRL(bj->format_label());
	PRL(dphi_tidal);
	print_binary_from_dyn_pair(bi, bj, 0, 0, true, true);
	return true;
    }

    //-----------------------------------------------------------------

    // Offset components to the center of mass frame.

    bi->inc_pos(-cm->get_pos());
    bi->inc_vel(-cm->get_vel());
    bj->inc_pos(-cm->get_pos());
    bj->inc_vel(-cm->get_vel());

    hdyn *p;
    if (VERBOSE) {
      cerr << endl << "too close: " << bi->format_label();
      cerr << " and " << bj->format_label() << endl;
      p = bi->get_parent();
      pp(p), cerr << endl;
    }

    // Restructure the tree.

    add_node_before(cm, bi);
    detach_node_from_general_tree(bi);
    detach_node_from_general_tree(bj);

    cm->set_oldest_daughter(bi);
    bi->set_parent(cm);
    bj->set_parent(cm);
    bi->set_elder_sister(NULL);
    bi->set_younger_sister(bj);
    bj->set_elder_sister(bi);
    bj->set_younger_sister(NULL);
    label_binary_node(cm);

    // Kira convention is that kepler is associated with components.
    // Really makes more sense to tie kepler to the CM, but for now
    // we won't do that here...  Criterion here is that the relevant
    // kepler is attached to (only) the first component.

    bi->set_kepler(k);
    bi->set_t_pred(2*peri_time-time);	// end of unperturbed segment

    // Eventually we will allow possibility of completely unperturbed
    // motion extending over several orbits, as in kira.

    if (VERBOSE) {
      cerr << "created new binary " << cm->format_label() << endl;
      pp(p), cerr << endl;
      // k->print_all();
      real dt_unpert = bi->get_t_pred() - time;
      PRC(dt_unpert); PRL(bi->get_t_pred());
    }

    // Make sure pos and pred_pos agree.

    cm->set_t_pred(time);
    for_all_nodes(hdyn, cm, bb) {
	bb->set_pred_pos(bb->get_pos());
	bb->set_pred_vel(bb->get_vel());
    }

    return true;
}



local inline bool terminate_binary(hdyn*& bi)
{
    if (!bi) return false;
    hdyn *od = bi->get_oldest_daughter();
    if (!od) return false;
    hdyn *yd = od->get_younger_sister();
    if (!yd) return false;

    kepler *k = od->get_kepler();
    if (!k) err_exit("smallN: binary with no kepler.");

    hdyn *p;
    if (VERBOSE) {
      cerr << endl << "terminating binary " << bi->format_label() << endl;
      // k->print_all();
      PRL(od->get_t_pred());
      p = bi->get_parent();
      pp(p), cerr << endl;
    }

    // Offset the components to include the center of mass pos and vel.

    od->inc_pos(bi->get_pos());
    od->inc_vel(bi->get_vel());
    yd = od->get_younger_sister();
    yd->inc_pos(bi->get_pos());
    yd->inc_vel(bi->get_vel());

    od->set_kepler(NULL);
    delete k;

    // Recompute the change in the tidal potential (components - CM) due
    // to the rest of the system.

    real dphi_tidal = 0;
    for_all_daughters(hdyn, bi->get_parent(), bk)
	if (bk != bi) {
	    real dphi_k = -od->get_mass()/abs(od->get_pos()-bk->get_pos())
			  -yd->get_mass()/abs(yd->get_pos()-bk->get_pos())
			  +bi->get_mass()/abs(bi->get_pos()-bk->get_pos());
	    dphi_tidal += bk->get_mass()*dphi_k;
	}

    if (VERBOSE) {
      PRC(dphi_tidal); PRC(od->get_pot());
    }
    dphi_tidal -= od->get_pot();
    if (VERBOSE)
      PRL(dphi_tidal);

    // Clear flags.

    od->set_pot(0);
    yd->set_pot(0);

    // Absorb dphi_total into the binary energy -- it represents the
    // omitted effect of the tidal field.  For now, simply rescale
    // the velocity.  We will investigate later how best to correct
    // the orbital elements.

    real vrel2 = square(yd->get_vel()-od->get_vel());

    // 0.5 * mu * vrel2 + phi = const
    // 0.5 * mu * d(vrel2) = -dphi
    // d(vrel2) = -2*dphi/mu
    // vrel2 *= 1 - dphi/(0.5*mu*vrel2)

    real mu = od->get_mass()*yd->get_mass()/bi->get_mass();
    real vfac = sqrt(1 - 2*dphi_tidal/(mu*vrel2));

    od->set_vel(vfac*od->get_vel());
    yd->set_vel(vfac*yd->get_vel());

    if (VERBOSE) {
      cerr << "corrected component velocities: "; PRL(vfac);
    }

    // Update the tree.

    add_node_before(od, bi);
    add_node_before(yd, bi);
    detach_node_from_general_tree(bi);

    // Delete the CM node (careful to detach the daughters to
    // avoid recursive deletion!).

    bi->set_oldest_daughter(NULL);
    delete(bi);

    if (VERBOSE) {
      pp(p), cerr << endl;
    }

    // Make sure pos and pred_pos agree.  Note that init_pred sets
    // t_pred = time as well as pred_pos = pos.

    od->init_pred();
    yd->init_pred();

    // Change bi to yd, so next node is correct in loop below.

    bi = yd;
    return true;
}



local inline real acc_and_jerk_and_get_dt(hdyn *b)
{
    // Compute acc and jerk, and return the top-level time step.

    real dt = calculate_top_level_acc_and_jerk(b);

    // Adjust dt to take unperturbed systems into account.

    real t = b->get_time();
    for_all_daughters(hdyn, b, bi)
	if (bi->is_parent()) correct_timestep_for_components(bi, t, dt);

    // Don't store the time steps yet, as the old time step may
    // be needed for the corrector.

    return dt;
}

local inline void set_all_timesteps(hdyn *b, real dt)
{
    // Update all time steps (root and top-level nodes).

    for_all_daughters(hdyn, b, bi) bi->set_timestep(dt);
    b->set_timestep(dt);
}

local void set_all_times(hdyn *b)
{
    // Update times of all nodes in the system from system_time.

    real sys_t = b->get_time();
    for_all_nodes(hdyn, b, bi) bi->set_time(sys_t);
}



// Routines transplanted from sdyn/evolve/sdyn_ev.C.  Naming scheme is
// that pred_pos here is new_pos there (similarly vel), old_acc here is
// acc there, acc here is new_acc there (similarly jerk).

local void correct_acc_and_jerk(hdyn *bi,
				const real new_dt, const real prev_new_dt)
{
    // Correct the values of acc and jerk from time prev_new_dt to new_dt.
    // We simply fit a polynomial to old_acc and old_jerk at time 0 and acc
    // and jerk at time prev_new_dt, then evaluate it at time new_dt.
    // Conventions and notation for time steps are taken from sdyn(3).

    if (new_dt == prev_new_dt) return;

    real dt_off = new_dt - 0.5 * prev_new_dt;
                                    // offset from midpoint of prev_new_dt step
    real theta = 0.5 * prev_new_dt;
    real tau = dt_off / theta;          // equals 1 if new_dt = prev_new_dt

    real inv_theta = 1 / theta;
    real tau2 = tau * tau;
    real tau3 = tau2 * tau;

    // Wouldn't need these if this were a member function...

    vec old_acc = bi->get_old_acc();
    vec acc = bi->get_acc();
    vec old_jerk = bi->get_old_jerk();
    vec jerk = bi->get_jerk();

    vec prev_acc = acc;
    vec prev_jerk = jerk;

    acc = 0.25 * (old_acc * (2 - 3 * tau + tau3)
		   + prev_acc * (2 + 3 * tau - tau3)
		   + old_jerk * theta * (1 - tau - tau2 + tau3)
		   + prev_jerk * theta * (-1 - tau + tau2 + tau3));

    jerk = 0.25 * (old_acc * inv_theta * 3 * (-1 + tau2)
		    + prev_acc * inv_theta * 3 * (1 - tau2)
		    + old_jerk * (-1 - 2*tau + 3*tau2)
		    + prev_jerk * (-1 + 2*tau + 3*tau2));

    // ...or these.

    bi->set_acc(acc);
    bi->set_jerk(jerk);
}

local void correct_pos_and_vel(hdyn *bi, const real new_dt)
{
    // Apply a corrector in the form presented by Hut et al. (1995).
    // The "pred" quantities are those at the end of the step.

    real new_dt2 = new_dt * new_dt;
    real new_dt3 = new_dt2 * new_dt;

    // Wouldn't need these if this were a member function...

    vec old_acc = bi->get_old_acc();
    vec acc = bi->get_acc();
    vec old_jerk = bi->get_old_jerk();
    vec jerk = bi->get_jerk();

    vec vel = bi->get_vel();
    vec pred_vel = bi->get_pred_vel();

    // Corrector as in Hut et al. (1995).

    pred_vel = vel + new_dt * (acc + old_acc)/2
		   - new_dt2 * (jerk - old_jerk)/12;
    bi->set_pred_vel(pred_vel);

    bi->set_pred_pos(bi->get_pos() + new_dt * (pred_vel + vel)/2 
		     - new_dt2 * (acc - old_acc)/12);
}



local void print_pred_energies(hdyn *b)
{
    // Compute and print energies based on "pred" quantities,
    // for debugging purposes.

    int n = 0;
    for_all_nodes(hdyn, b, bb) n++;
    vec *p = new vec[n];
    vec *v = new vec[n];

    n = 0;
    for_all_nodes(hdyn, b, bb) {
	p[n] = bb->get_pos();
	bb->set_pos(bb->get_nopred_pos());
	v[n] = bb->get_vel();
	bb->set_vel(bb->get_nopred_vel());
	n++;
    }

    int pr = cerr.precision(HIGH_PRECISION);
    cerr << b->get_time() << " (pred): ";
    cerr.precision(pr);
    print_recalculated_energies((dyn*)b, 1);

    n = 0;
    for_all_nodes(hdyn, b, bb) {
	bb->set_pos(p[n]);
	bb->set_vel(v[n]);
	n++;
    }
}



// Take a step to time t.  Return the actual system time at the end
// of the step, after symmetrization if specified.

local real take_a_step(hdyn *b,		// root node
		       real &dt)	// natural step at start/end
{
    // Predict to time t + dt.

    real t = b->get_time();

#if 0

    if (t > 121.934) {

	eta = 0.01;
	gamma2_unpert = 1.e-10;
	PRC(eta); PRL(gamma2_unpert);

 	print_pred_energies(b);
	int pr = cerr.precision(HIGH_PRECISION);
	cerr << "t = " << t << "  E = "
	     << get_total_energy((dyn*)node_with_name("3", b),
				 (dyn*)node_with_name("5", b)) << endl;
	print_binary_from_dyn_pair((dyn*)node_with_name("3", b),
				   (dyn*)node_with_name("5", b),
				   0, 0, true, true);
	cerr.precision(pr);
    }

#endif

    for_all_daughters(hdyn, b, bi) {
	bi->store_old_force();
	bi->predict_loworder(t+dt);
	if (bi->is_parent()) advance_components_to_time(bi, t+dt);
    }

    // Compute forces, etc. and correct, iterating if desired.
    // Note that dt is the natural time step associated with the
    // state of the system on entry and the natural time step of
    // the new system on exit.   We use "pred" quantities to
    // represent the current iterate of quantities the end of the
    // time step (even after correction; these values are called
    // "new" in the sdyn and sdyn3 versions of the symmetrization
    // scheme.)  The acc and jerk at the start of the step are
    // "old_acc" etc.  Those at the pred time are acc and jerk.

    real new_dt = dt;
    real end_point_dt = dt;

    // new_dt will be the actual step.
    // end_point_dt will be the natural time step at the end
    // (apparently don't need to recompute after final correction).

    for (int i = 0; i <= n_iter; i++) {

	real prev_new_dt = new_dt;
	end_point_dt = acc_and_jerk_and_get_dt(b);

	// Obtain the next iterate of the time step.

	if (i < n_iter) {
	    new_dt = 0.5 * (dt + end_point_dt);
	    new_dt = dt + 0.5 * (end_point_dt - dt) * (new_dt/prev_new_dt);
	}

	// Extrapolate acc and jerk to the end of the new step, and
	// apply the corrector for this iteration.

	for_all_daughters(hdyn, b, bi) {
	    if (new_dt != prev_new_dt)
		correct_acc_and_jerk(bi, new_dt, prev_new_dt);
	    correct_pos_and_vel(bi, new_dt);
	    if (bi->is_parent()) advance_components_to_time(bi, t+new_dt);
	}
     }

    // Complete the step.

    for_all_daughters(hdyn, b, bi) {
	bi->set_pos(bi->get_pred_pos());
	bi->set_vel(bi->get_pred_vel());
	if (bi->is_parent())
	    update_components_from_pred(bi);
    }

    dt = end_point_dt;
    return t + new_dt;
}



local void set_tag(hdyn *b, const char *s)
{
    putiq(b->get_log_story(), s, 1);
}

local bool check_tag(hdyn *b, const char *s)
{
    return (find_qmatch(b->get_log_story(), s) != NULL);
}

local void clear_tag(hdyn *b, const char *s)
{
    rmq(b->get_log_story(), s);
}

local void combine_nodes(hdyn *bi, hdyn *bj)
{
    // Merge bi and bj into bi, and call the node (bi, bj).

    real mi = bi->get_mass(), mj = bj->get_mass(), m = mi + mj;
    mi /= m;
    mj /= m;
    bi->set_mass(m);
    bi->set_pos(mi*bi->get_pos()+mj*bj->get_pos());
    bi->set_vel(mi*bi->get_vel()+mj*bj->get_vel());
    bi->set_name(construct_binary_label(bi, bj));
    detach_node_from_general_tree(bj);
    delete bj;
    set_tag(bi, "merged");
}

local bool search_for_multiples(hdyn* b,
				real energy_cutoff = 0)	// positive ==> bound!
{
    // Recursively search for bound binary pairs among top-level nodes.
    // After each pass, merge bound systems and repeat the search.

    int n = b->n_leaves();
    if (n <= 2) return false;

    hdynptr *imerge = new hdynptr[n];
    hdynptr *jmerge = new hdynptr[n];
    real *emerge = new real[n];
    int nmerge = 0;

    bool found = false;

    for_all_daughters(hdyn, b, bi) {
// 	PRL(bi);
// 	PRL(bi->format_label());
// 	PRL(bi->get_log_story());
// 	PRL(find_qmatch(bi->get_log_story(), "printed"));
	clear_tag(bi, "printed");
	bi->set_d_nn_sq(VERY_LARGE_NUMBER);
    }

    for_all_daughters(hdyn, b, bi)
	for (hdyn *bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {  // now bi precedes bj

	    vec dx = bj->get_pos() - bi->get_pos();
	    real dr2 = square(dx);

	    if (dr2 < bi->get_d_nn_sq()) {
		bi->set_d_nn_sq(dr2);
		bi->set_nn(bj);
	    }
	    if (dr2 < bj->get_d_nn_sq()) {
		bj->set_d_nn_sq(dr2);
		bj->set_nn(bi);
	    }

	    real M = bi->get_mass() + bj->get_mass();
	    vec dv = bj->get_vel() - bi->get_vel();
	    real E = 0.5*dv*dv - M / abs(dx);		  // energy/reduced mass

	    // Convention: energy_cutoff is in terms of E.

	    if (E < -energy_cutoff) {

		// Note: pair will be identified as (bi, bj).

		print_binary_from_dyn_pair(bi, bj,
					   0, vec(0), true,
					   true);
		found = true;
		set_tag(bi, "printed");
		set_tag(bj, "printed");

		imerge[nmerge] = bi;
		jmerge[nmerge] = bj;
		emerge[nmerge++] = E;

		// Until we figure out the problem with stories and
		// copies of the system, content ourselves here with
		// printing out the nearest-neighbor distance.

		vec cmpos = (bi->get_mass()*bi->get_pos()
		    		+ bj->get_mass()*bj->get_pos())
		    	    / (bi->get_mass()+bj->get_mass());
		real dnn2 = VERY_LARGE_NUMBER;
		hdyn *nn = NULL;

		for_all_daughters(hdyn, b, bk)
		    if (bk != bi && bk != bj) {
			real dr2 = square(bk->get_pos() - cmpos);
			if (dr2 < dnn2) {
			    dnn2 = dr2;
			    nn = bk;
			}
		    }

		if (nn)
		    cerr << "    nearest neighbor is " << nn->format_label()
			 << " at distance " << sqrt(dnn2) << endl;
		cerr << endl;
	    }
	}

    for_all_daughters(hdyn, b, bi)
	if (check_tag(bi, "merged") && !check_tag(bi, "printed")) {

	    // Want info on the nearest neighbor, even if unbound.

	    if (bi->get_nn())
		print_binary_from_dyn_pair(bi, bi->get_nn(),
					   0, vec(0), true,
					   true);
	}

#if 0

    // Merge bound pair(s) and repeat.

    PRL(nmerge);

    if (found) {

	// Take care to merge to most tightly bound pairs first.
	// Sort by energy and remove duplicates.

	// Use a simple insertion sort (small N!).

//	for (int i = 0; i < nmerge; i++) {PRC(i); PRL(emerge[i]);}

	for (int j = 1; j < nmerge; j++) {
	    hdyn *bi = imerge[j];
	    hdyn *bj = jmerge[j];
	    real e = emerge[j];
	    int i = j - 1;
	    while (i >= 0 && emerge[i] > e) {
		imerge[i+1] = imerge[i];
		jmerge[i+1] = jmerge[i];
		emerge[i+1] = emerge[i];
		i--;
	    }
	    imerge[i+1] = bi;
	    jmerge[i+1] = bj;
	    emerge[i+1] = e;
	}

//	PRL(nmerge);
//	for (int i = 0; i < nmerge; i++) {PRC(i); PRL(emerge[i]);}

	// Check for duplicates (keep first instance).  There should be
	// at most one instance of any node on either list.  Work backwards
	// through the list, eliminating pairs which share a component with
	// a more tightly bound pair.

	for (int i = nmerge-2; i >= 0; i--) {
	    int joff = 0;
	    for (int j = i+1; j < nmerge; j++) {
		if (joff > 0) {
		    imerge[j] = imerge[j+joff];
		    jmerge[j] = jmerge[j+joff];
		    emerge[j] = emerge[j+joff];
		}
		if (imerge[j] == imerge[i] || imerge[j] == jmerge[i]
		    || jmerge[j] == imerge[i] || jmerge[j] == jmerge[i]) {
		    joff++;
		    nmerge--;
		}
	    }
	}

	PRL(nmerge);
//	for (int i = 0; i < nmerge; i++) PRL(emerge[i]);

	// Do the mergers, in order.

	for (int i = 0; i < nmerge; i++) {
	    cerr << "merging " << imerge[i]->format_label();
	    cerr << " and " << jmerge[i]->format_label() << endl;
	    combine_nodes(imerge[i], jmerge[i]);
	}

//	search_for_multiples(b, energy_cutoff);
    }

#endif

    delete [] imerge;
    delete [] jmerge;
    delete [] emerge;

    return found;
}



local void copy_hdyn(hdyn *b, hdyn *c)
{
    // Copy hdyn contents of b to c.

    c->set_mass(b->get_mass());
    c->set_pos(b->get_pos());
    c->set_vel(b->get_vel());
    c->set_index(b->get_index());
    c->set_name(b->get_name());
}

local hdyn *copy_tree(hdyn *b)
{
    // Make a local hdyn copy of the (flat) tree under hdyn node b.

    cerr << "in copy_tree..." << endl << flush;
    hdyn *c = new hdyn();
    PRL(c);
    copy_hdyn(b, c);
    c->set_name("copytop");

    hdyn *pcc = NULL;
    for_all_daughters(hdyn, b, bb) {

	//PRL(bb->format_label());
	put_node(bb, cerr);

	// Pointers:

	hdyn *cc = new hdyn();
	cc->set_parent(c);
	if (!c->get_oldest_daughter()) c->set_oldest_daughter(cc);

	if (pcc) {
	    cc->set_elder_sister(pcc);
	    pcc->set_younger_sister(cc);
	}
	pcc = cc;
	PRL(cc);

	// Content:

	copy_hdyn(bb, cc);
	//put_node(cc, cerr);
    }

    return c;
}



static bool local_sys_stats = false;

local void sys_stats(hdyn *b)
{
    if (local_sys_stats) {

	// A local, minimal version of the sys_stats function.

	cerr << endl << "SmallN time = " << b->get_time();
	if (smallN_time_scale > 0)
	    cerr << " = " << (real)b->get_time()/smallN_time_scale
		 << " smallN units";

	cerr << endl;
	pp(b); cerr << endl;
	print_recalculated_energies((dyn*)b, 1);

	// Numbers:

	int N = b->n_leaves();
	int N_top_level = b->n_daughters();
	cerr << "    N = " << N << "  N_top_level = " << N_top_level << endl;

	// Masses and averages:

	real total_mass_nodes = 0;
	for_all_daughters(hdyn, b, bi)
	    total_mass_nodes += bi->get_mass();

	real total_mass_leaves = 0;
	real m_min = VERY_LARGE_NUMBER, m_max = -m_min;
	for_all_leaves(hdyn, b, bj) {
	    total_mass_leaves += bj->get_mass();
	    m_min = Starlab::min(m_min, bj->get_mass());
	    m_max = Starlab::max(m_max, bj->get_mass());
	}
	real m_av = total_mass_leaves / Starlab::max(1, N);

	cerr << "    total_mass = " << b->get_mass();
	cerr << "  nodes: " << total_mass_nodes;
	cerr << "  leaves: " << total_mass_leaves << endl;

	cerr << "    m_min = " << m_min << "  m_max = " << m_max
	     << "  m_av = " << m_av << endl;

	// Center of mass:

	vec com_pos = 0, com_vel = 0;
//	cerr << "to compute_com" << endl << flush;
	compute_com(b, com_pos, com_vel);
//	cerr << "back" << endl << flush;

	cerr << "    center of mass position = " << com_pos << endl
	     << "                   velocity = " << com_vel << endl;

	// Basic structure:

//	hdyn *c = copy_tree(b);

	// For unknown reasons, the combination of copying the system
	// and using search_for_multiples leads to problems with the
	// story mechanism.  Apparently related to memory management,
	// but not clear what the interaction is.  Need both copy_tree
	// *and* search_for_multiples for the failure to occur...

//	if (c) {
	    cerr << endl << "Bound pairs:" << endl;
	    if (!search_for_multiples(b, smallN_energy_scale))
	      cerr << "   (none)" << endl << endl;
//	    rmtree(c);
//	}

    } else

	// The usual version.

	sys_stats(b, 1, 1, true, false, 2, true, true, true);
}



local bool fully_unperturbed(hdyn *b)
{
    // Check if the binary (bi_min, bj_min) is perturbed by bk_nn.
    // Use the same criteria as in the rest of kira.

    if (!bk_nn) return false;

    // cerr << "checking full unperturbed" << endl;

    real mi = bi_min->get_mass();
    real mj = bj_min->get_mass();
    real mass = mi + mj;
    vec cm_pos = (mi*bi_min->get_pos() + mj*bj_min->get_pos()) / mass;
    vec cm_vel = (mi*bi_min->get_vel() + mj*bj_min->get_vel()) / mass;

    // The next nearest body is bk_nn.  If it is a component
    // of a close binary, then use the CM of that binary in
    // assessing the unperturbed status of (bi_min,bj_min).

    real m_nn  = bk_nn->get_mass();
    vec nn_pos = bk_nn->get_pos();
    vec nn_vel = bk_nn->get_vel();

    hdyn *nnn = NULL;
    real r2_nnn = VERY_LARGE_NUMBER;

    // Find the nearest neighbor of bk_nn (if any).

    for_all_daughters(hdyn, b, bb)
	if (bb != bi_min && bb != bj_min && bb != bk_nn) {
	    real r2 = square(bb->get_pos() - nn_pos);
	    if (r2 < r2_nnn) {
		r2_nnn = r2;
		nnn = bb;
	    }
	}

    if (nnn) {

	// Check the relative energy of (bk_nn, nnn).  If less than
	// smallN_energy_scale, replace nn_* by center of mass quantities.

	real M = bk_nn->get_mass() + nnn->get_mass();
	vec dx = bk_nn->get_pos() - nnn->get_pos();
	vec dv = bk_nn->get_vel() - nnn->get_vel();
	if (0.5*dv*dv - M / abs(dx) < -smallN_energy_scale) {
	    m_nn = M;
	    nn_pos = (bk_nn->get_mass()*bk_nn->get_pos()
		        + nnn->get_mass()*nnn->get_pos()) / M;
	    nn_vel = (bk_nn->get_mass()*bk_nn->get_vel()
		        + nnn->get_mass()*nnn->get_vel()) / M;
	}
    }

    vec dx = nn_pos - cm_pos;
    vec dv = nn_vel - cm_vel;

    // Is the neighbor moving away from the binary and perturbing
    // it negligibly?

    if (dx*dv > 0) {

	// Apply the perturbation test basically as in kira (with a
	// safety factor).  Use the binary semi-major axis as a
 	// representative binary separation, and compute the squared
	// perturbation using just the nearest neighbor (presumably
	// the dominant component).

//	real rbin = abs(bi_min->get_pos() - bj_min->get_pos());
	real rbin = 2*get_semi_major_axis(bi_min, bj_min);
	real rnn  = abs(dx);
	real perturbation_squared = square(2*(m_nn/mass)*pow(rbin/rnn, 3));

//	PRC(rbin); PRC(rnn); PRL(perturbation_squared);

	if (perturbation_squared
	      < b->get_kira_options()->full_merge_tolerance) return true;
    }

    return false;
}



// Evolve the system to time t_end, using the input data and settings
// without modification.  Stop if
//
//	(1) t >= t_end,
//	(2) any particle gets too far (break_r2) from the origin, or
//	(3) (optionally) if the closest binary is unperturbed.
//
// The return value is 1, 2, or 3 for these three cases.

#define NCHECK 100

int smallN_evolve(hdyn *b,
		  real t_end,		// default = VERY_LARGE_NUMBER
		  real break_r2,	// default = VERY_LARGE_NUMBER
		  real end_on_unpert,	// default = false
		  real dt_log,		// default = 0
		  real dt_energy,	// default = 0
		  real dt_snap)		// default = 0
{
    b->init_pred();
    for_all_daughters(hdyn, b, bi) bi->init_pred();

    int n_steps = 0;

    real dt = acc_and_jerk_and_get_dt(b);	// initialization
    set_all_timesteps(b, dt);
    set_all_times(b);

    // Note that we correctly maintain all times and time steps, but
    // we do not modify system_time.  It is not strictly necessary to
    // update all times, as times and time steps are shared, but it is
    // convenient to keep the system coherent for output purposes.

    real t_log = b->get_time();
    real t_energy = b->get_time();
    real t_snap = b->get_time();
    int  n_snap = 0;

    bk_nn = NULL;

    if (dt_energy > 0) {
	int p = cerr.precision(HIGH_PRECISION);
	cerr << b->get_time() << " (" << n_steps << "): ";
	cerr.precision(p);
	print_recalculated_energies((dyn*)b, 0);
	if (dt_energy > 0) t_energy += dt_energy;
    } else
	initialize_print_energies(b);

    // Optional output:

    if (dt_snap > 0 || dt_log > 0) {
	real t_sys = b->get_system_time();
	b->set_system_time(b->get_time());
	if (dt_log > 0) {
	    sys_stats(b);
	    t_log += dt_log;
	}
	if (dt_snap > 0) {
	    cerr << "Snapshot " << n_snap << " output at time "
		 << b->get_time() << " (" << n_steps << ")" << endl;
	    put_node(b);
	    t_snap += dt_snap;
	}
	b->set_system_time(t_sys);
    }

    // Termination criteria:
    //
    //     (1) t >= t_end  (checked at start of the while loop).
    //     (2) optional: r > break_r (checked after every NCHECK steps)
    //     (3) optional: inner binary is unperturbed (checked after each step)

    while (b->get_time() < t_end) {

	// Take a step.  Don't set time in advance, as the symmetrization
	// process will determine the actual time step.  Thus, during the
	// entire step, all times refer to the time at the *start* of the
	// step.  The time step dt on entry is the natural time step for
	// the system at the current time.  On return it is the new natural
	// time step for the system at the end of the step.  The return
	// value of the function is the new system time.

//      PRC(b->get_time()); PRC(dt); cerr << "a" << endl;
// 	int p = cerr.precision(18);
// 	for_all_daughters(hdyn, b, bb)
// 	  cerr << " " << bb->format_label() << " " << bb->get_pos() << endl;
// 	cerr.precision(p);

	b->set_time(take_a_step(b, dt));

	// Update all times and time steps.

	set_all_timesteps(b, dt);
	set_all_times(b);
	n_steps++;

// 	PRC(b->get_time()); cerr << "b" << endl;
// 	p = cerr.precision(18);
// 	for_all_daughters(hdyn, b, bb)
// 	  cerr << " " << bb->format_label() << " " << bb->get_pos() << endl;
// 	cerr.precision(p);

#if 0
	cerr << b->get_time() << " (" << n_steps << "): ";
	print_recalculated_energies((dyn*)b, 1);
#endif

	// The step is over and all times have been advanced.
	// Time step dt was set by bi_min and bj_min during the last acc
	// and jerk calculation.

	// Debugging (rbin can only be set later in this function):

	if (rbin > 0 && abs(bi_min->get_pos()-bj_min->get_pos()) > rbin) {

	    // Print debugging info on the binary we didn't actually create.

	    hdyn * od = bi_min;
	    hdyn * yd = bj_min;
	    PRL(od->format_label());
	    PRL(yd->format_label());
	    real mbin = od->get_mass() + yd->get_mass();
    	    vec cmpos = (od->get_mass()*od->get_pos()
			 + yd->get_mass()*yd->get_pos())/mbin;
	    real dphi_tidal = 0;
	    for_all_daughters(hdyn, od->get_parent(), bk)
		if (bk != od && bk != yd) {
		    real dphi_k =
			-od->get_mass()/abs(od->get_pos()-bk->get_pos())
			-yd->get_mass()/abs(yd->get_pos()-bk->get_pos())
			+mbin/abs(cmpos-bk->get_pos());
		    dphi_tidal += bk->get_mass()*dphi_k;
		}
	    PRC(dphi_tidal); PRC(od->get_pot());
	    dphi_tidal -= od->get_pot();
	    PRL(dphi_tidal);
	    print_binary_from_dyn_pair(od, yd, 0, 0, true, true);

	    rbin = 0;
	}

	// Check second (size) termination criterion.

	if (n_steps%NCHECK == 0) {
	    for_all_daughters(hdyn, b, bi)
	        if (square(bi->get_pos()) > break_r2)
		  return 2;
	}

	// Check for the start/end of unperturbed motion, and also the
	// third (unperturbed) termination condition.  Use various
	// thresholds to avoid this check at the end of every step:
	//
	//	time step dt < dt_crit
	//	bi_min and bj_min are leaves
	//	bi_min and bj_min are within distance r_crit
	//	bi_min and bj_min are approaching
	//	others...
	//
	// Once the check fails, a more clever search would defer
	// further checks until the components had approached each
	// other significantly, but this would entail significant
	// bookeeping.  For now, just live with the fact that we carry
	// out this check at roughly half the total number of staps,
	// as the initial hard binary will likely dominate the time
	// step.

	bool tree_changed = false;

	if (dt < dt_crit
	    && unperturbed_and_approaching(b, bi_min, bj_min)) {

	    // We have a candidate for unperturbed motion.  Check if
	    // it satisfies the unperturbed termination condition.

	    if (end_on_unpert && fully_unperturbed(b)) return 3;

	    // Keep going -- make a new binary.

	    tree_changed = create_binary(bi_min, bj_min);

	    // PRL(bk_nn->format_label());

	    // Debugging:

	    if (!dbg_allow_binary && tree_changed) {

		// Didn't actually create a binary.  Save the separation
		// for later use.

		rbin = abs(bi_min->get_pos()-bj_min->get_pos());
		tree_changed = false;
	    }

	    if (tree_changed && end_on_unpert) {

		// See if the new binary triggers the termination criterion.
		// Nearest neighbor is bk_nn.

		hdyn *cm = bi_min->get_parent();
		kepler *k = bi_min->get_kepler();
		if (VERBOSE) {
		  PRC(cm); PRL(k);
		}

		if (k->get_energy() < 0) {

		    // Note that estimated_perturbation is a kira
		    // function (in hdyn_inline.C).

		    real gamma = estimated_perturbation(bi_min->get_parent(),
							bk_nn);
		    if (VERBOSE) PRL(gamma);

		    // Conservatively advance the perturbation to apastron.

		    gamma *= pow(k->get_apastron()/k->get_separation(), 3);
		    if (VERBOSE) PRL(gamma);
		    vec dr = hdyn_something_relative_to_root(cm, &hdyn::get_pos)
		       - hdyn_something_relative_to_root(bk_nn, &hdyn::get_pos);
		    if (VERBOSE) PRL(dr);
		    vec dv = hdyn_something_relative_to_root(cm, &hdyn::get_vel)
		       - hdyn_something_relative_to_root(bk_nn, &hdyn::get_vel);
		    if (VERBOSE) {
		      PRC(dv); PRL(dr*dv);
		    }

		    // Exit if the binary would be unperturbed and
		    // approaching by (default) kira standards.

		    if (gamma < DEFAULT_FULL_MERGE_TOLERANCE && dr*dv > 0)
			return 3;
		}
	    }
	}

	// Check for termination of unperturbed motion.

	for_all_daughters(hdyn, b, bi) {
	    hdyn *od = NULL;
	    if (od = bi->get_oldest_daughter())
		if (od->get_t_pred() <= b->get_time()) {
// 		    int p = cerr.precision(HIGH_PRECISION);
// 		    cerr << b->get_time() << " (xxxxx): ";
// 		    cerr.precision(p);
// 		    print_recalculated_energies((dyn*)b, 1);
		    tree_changed |= terminate_binary(bi);
		}
	}

	if (tree_changed) {

	    // Recompute accs, jerks, and the time step.

	    dt = acc_and_jerk_and_get_dt(b);

	    if (VERBOSE) {
	      cerr << endl << "new structure: "; pp(b);
	      cerr << ",  ";
	      PRL(dt);
	    }

	    set_all_timesteps(b, dt);
// 	    int p = cerr.precision(HIGH_PRECISION);
// 	    cerr << b->get_time() << " (yyyyy): ";
// 	    cerr.precision(p);
// 	    print_recalculated_energies((dyn*)b, 1);
	}

	// Optional output:

	if (dt_log > 0 && b->get_time() >= t_log) {
	    real t_sys = b->get_system_time();
	    b->set_system_time(b->get_time());
	    sys_stats(b);
	    t_log += dt_log;
	    b->set_system_time(t_sys);
	}

	if (dt_energy > 0 && b->get_time() >= t_energy) {
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << b->get_time() << " (" << n_steps << "): ";
	    cerr.precision(p);
	    print_recalculated_energies((dyn*)b, 1);
	    t_energy += dt_energy;
	}

	if (dt_snap > 0 && b->get_time() >= t_snap) {
	    real t_sys = b->get_system_time();
	    b->set_system_time(b->get_time());
	    cerr << "Snapshot " << ++n_snap << " output at time "
		 << b->get_time() << " (" << n_steps << ")" << endl;
	    put_node(b);
	    t_snap += dt_snap;
	    b->set_system_time(t_sys);
	}
    }

    // PRL(n_steps);
    return 1;
}



local void old_restore_binary_tree(hdyn *b)	// to be deleted
{
    // Reinstate the binary tree structure below node b.
    //
    // Tree structure is dictated by the decomposition leading to
    // the end of integrate_multiple():
    //
    //	    (((bi_min, bj_min), bk_nn [+ partner, if any]) rest...)
    //
    // Build the new subtree from the bottom up based on this.
    // All positions, velocities, accs and jerks are correct.
    // Time steps are also OK (if conservative).

    // Strategy: order the members in an array, starting with bi_min,
    // bj_min, and bk_nn.  Each subsequent entry is the nearest
    // neighbor of the previous entry, drawn from the remaining stars.

    cerr << endl << "in restore_binary_tree()" << endl;

    b->flatten_node();			// just in case... maybe better
					// to require no unpert motion
    int n = b->n_daughters();

#if 0

    // Special care because energy functions expect low-level
    // nodes to be binary nodes.

    hdyn *tmp1 = b->get_parent();
    hdyn *tmp2 = b->get_root();
    b->set_parent(NULL);
    b->set_root(b);
    pp(b);
    print_recalculated_energies((dyn*)b);
    b->set_parent(tmp1);
    b->set_root(tmp2);
    pp3(b);

#endif

    hdynptr *list = new hdynptr[n];

    // bi_min and bj_min are computed at every step; bk_nn is only
    // computed if a check is made for unperturbed motion.  It should
    // be NULL or usable.

    list[0] = bi_min;
    list[1] = bj_min;

    int i = 2;
    if (bk_nn)
      list[2] = bk_nn;
    else
      i = 1;

    while (i < n-1) {

	// Add the remaining components (potential better than distance?).

	vec last_pos = list[i]->get_pos();
	real rnn2 = VERY_LARGE_NUMBER;
	hdyn *nn = NULL;
	for_all_daughters(hdyn, b, bb) {

	    bool counted = false;
	    for (int j = 0; j <= i; j++)
		if (bb == list[j]) {
		    counted = true;
		    break;
		}

	    if (!counted) {

		real r2 = square(bb->get_pos() - last_pos);
		if (r2 < rnn2) {
		    rnn2 = r2;
		    nn = bb;
		}
	    }
	}

	if (nn) list[++i] = nn;
	else break;
    }

    if (i != n-1)
	cerr << "warning: restore_binary_tree: i = " << i
	     << " != n = " << n << endl;

    PR(n); cerr << ": ";
    for (i = 0; i < n; i++) {
	cerr << list[i]->format_label() << " ";
//	PRL(list[i]->get_pos());
    }
    cerr << endl;

    pp(b);

    // Now combine the elements of the list.  Each element becomes the
    // binary sister of the current node bcurr, and bcurr becomes the
    // center of mass node.

    hdyn *bcurr = list[0];				// = bi_min
    cerr << "start with " << bcurr->format_label() << endl;
    for (i = 1; i < n; i++) {

	// Combine list[i] with bcurr.  Code follows that in
        // create_binary_from_toplevel_nodes(), except that we
	// don't insist that the nodes be top-level.

	hdyn *bi = bcurr;
	hdyn *bj = list[i];

	detach_node_from_general_tree(bi);
	bi->set_younger_sister(NULL);
	bi->set_elder_sister(NULL);

	bcurr = new hdyn();				// new CM
	insert_node_into_binary_tree(bi, bj, bcurr);	// --> (bi, bj)

	label_binary_node(bi->get_parent());

	// Above operations only handle the tree structure.
	// Must also reorganize pos, vel, acc, and jerk, and reset
	// all pred_ and old_ quantities.  Handle these using
	// existing tool in hdyn_tt.C.

	bcurr->setup_binary_node();

	// Set pred times also.  This probably should be done in
	// setup_binary_node(), as the operations there only make
	// sense if the pred times are the same in any case.

	bcurr->set_t_pred(b->get_time());
	bi->set_t_pred(b->get_time());
	bj->set_t_pred(b->get_time());

	cerr << "added " << list[i]->format_label() << ": ";
	pp(bcurr); cerr << endl;
    }

    // Note: the final bcurr is the entire clump, but the dynamical
    // data it contains should all be zero (within rounding error).
    // Check this, and then connect bcurr's daughters directly to b.

    PRL(bcurr->get_pos());
    PRL(bcurr->get_vel());
    PRL(bcurr->get_acc());
    PRL(bcurr->get_jerk());

    hdyn *od = bcurr->get_oldest_daughter();
    b->set_oldest_daughter(od);
    od->set_parent(b);
    od->get_younger_sister()->set_parent(b);

    bcurr->set_oldest_daughter(NULL);
    delete bcurr;    

//    pp3(b);
//    print_recalculated_energies((dyn*)b);
    cerr << "leaving restore_binary_tree()" << endl << endl;
}

local void restore_binary_tree(hdyn *b)
{
    // Reinstate the binary tree structure below node b.
    //
    // Closest pair is known from the decomposition leading to
    // the end of integrate_multiple():
    //
    //	    (((bi_min, bj_min), bk_nn [+ partner, if any]) rest...)
    //
    // However, build the new subtree from scratch.  All positions,
    // velocities, accs and jerks are correct.  Time steps are also OK
    // (if somewhat conservative).

    // Strategy: starting from a flat tree, combine the closest
    // (largest potential) pair.  Then repeat, considering only
    // top-level noes, until the binary tree is complete..

    // cerr << endl << "in restore_binary_tree()" << endl;

    b->flatten_node();			// just in case... maybe better
					// to require no unpert motion

    hdyn *bcm = NULL;

    // Find and pair the particles with the most negative pairwise
    // potential (per unit reduced mass).

    while (1) {
	real potmin = 0;
	hdyn *bimin = NULL;		// local, note
	hdyn *bjmin = NULL;
	for_all_daughters(hdyn, b, bi) {
	    real mi = bi->get_mass();
	    vec posi = bi->get_pos();
	    for (hdyn *bj = bi->get_younger_sister(); bj != NULL;
		 bj = bj->get_younger_sister()) {
		real pot = -(mi+bj->get_mass())/abs(posi - bj->get_pos());
		if (pot < potmin) {
		    potmin = pot;
		    bimin = bi;
		    bjmin = bj;
		}
	    }
	}

	// cerr << "combining ";
	// PRC(bimin->format_label()); PRL(bjmin->format_label());

	// Combine nodes bimin and bjmin under cm bcm.  Code follows
	// that in create_binary_from_toplevel_nodes(), except that we
	// don't insist that the nodes be top-level.

	detach_node_from_general_tree(bimin);
	bimin->set_younger_sister(NULL);
	bimin->set_elder_sister(NULL);

	bcm = new hdyn();				  // new CM
	insert_node_into_binary_tree(bimin, bjmin, bcm);  // --> (i,j)

	// Above operations only handle the tree structure.  Must also
	// reorganize pos, vel, acc, and jerk, and reset all pred_ and
	// old_ quantities.  Handle these using the existing tool in
	// hdyn_tt.C.

	label_binary_node(bcm);
	bcm->setup_binary_node();

	// Set pred times also.  This probably should be done in
	// setup_binary_node(), as the operations there only make
	// sense if the pred times are the same in any case.

	bcm->set_t_pred(b->get_time());
	bimin->set_t_pred(b->get_time());
	bjmin->set_t_pred(b->get_time());

	// cerr << "under "; PRL(bcm->format_label());

	if (b->n_daughters() <= 1) break;
    }

    // Note: the final bcm is the entire clump, but the dynamical
    // data it contains should all be zero (within rounding error).
    // Check this, and then connect bcm's daughters directly to b.

//     PRL(bcm->get_pos());
//     PRL(bcm->get_vel());
//     PRL(bcm->get_acc());
//     PRL(bcm->get_jerk());

    hdyn *od = bcm->get_oldest_daughter();
    b->set_oldest_daughter(od);
    od->set_parent(b);
    od->get_younger_sister()->set_parent(b);

    bcm->set_oldest_daughter(NULL);
    delete bcm;    

    // cerr << "leaving restore_binary_tree()" << endl << endl;
}



local kepler get_multiple_params(hdyn *b,
				 hdynptr &inner)
{
    // Clump b is a binary tree describing a multiple system.  We
    // expect that it contains at least one very close pair of leaves,
    // already in the form of two binary sisters.  Find the closest
    // such pair, to set time and length scales.  Use both minimum
    // distance and maximum potential criteria in the search.  Return
    // a kepler structure describing the closest pair and set inner
    //pointing to the inner CM.

    hdyn *bimin, *bjmin;
    real rmin = VERY_LARGE_NUMBER;
    hdyn *bipot, *bjpot;
    real pmax = 0;

    for_all_leaves(hdyn, b, bi)
      if (bi->get_elder_sister() == NULL) {	// older sister
	hdyn *bj = bi->get_younger_sister();	// younger sister
	    if (bj && bj->is_leaf()) {
		real r = abs(bi->get_pos() - bj->get_pos());
		if (r < rmin) {
		    bimin = bi;
		    bjmin = bj;
		    rmin = r;
		}
		real pot = bi->get_mass()*bj->get_mass()/r;
		if (pot > pmax) {
		    bipot = bi;
		    bjpot = bj;
		    pmax = pot;
		}
	    }
	}

    // PRC(rmin); PRC(bimin->format_label()); PRL(bjmin->format_label());
    // PRC(pmax); PRC(bipot->format_label()); PRL(bjpot->format_label());

    if (bimin != bipot)
	cerr << "get_multiple_params: closest and max. potential pair "
	     << "are different.  Using the latter." << endl;

    // Use bipot and bjpot as the closest pair.

    bi_min = bipot;
    bj_min = bjpot;
    inner = bipot->get_parent();

    // Determine their binary parameters.

    real tsave = b->get_system_time();
    b->set_system_time(0);
    kepler k; 
    initialize_kepler_from_dyn_pair(k, bipot, bjpot, false);
    // k.print_all();
    b->set_system_time(tsave);

    return k;
}

local int correct_external_nodes(hdyn *b, hdynptr nlist[])
{
    // Correct all nodes potentially affected by the multiple
    // integration -- those having an nn, coll, or perturber which is
    // a node (not a leaf, since those are preserved) under b.  Most
    // likely these are already pointers to leaves, but check and
    // correct that here (i.e.  replaced pointers to nodes by pointers
    // to leaves before proceeding).  By construction, the inner
    // binary is now perturbed, so it should also appear as components
    // on all lists.  Thus it is sufficiently unlikely that any
    // centers of mass appear that we should print a warning message
    // flag when we encounter one.  Also flag nodes with pointers to
    // the innermost components of the clump, as these will soon
    // become pointers to the newly unperturbed center of mass
    // (implemented at the end) -- returned in nlist.

    int ilist = 0;
    for_all_nodes(hdyn, b->get_root(), bb) {

	if (!is_descendent_of(bb, b, 1)) {	// bb is not part of the clump
	    int add = 0;

	    hdyn *nn = bb->get_nn();
	    if (nn && nn->is_valid()
		&& nn->is_parent()
		&& is_descendent_of(nn, b, 0)) {

		// Replace nn by a component.

		real minr2 = VERY_LARGE_NUMBER;
		hdyn *nncpt = NULL;
		vec bbpos = hdyn_something_relative_to_root(bb, &hdyn::get_pos);
		for_all_leaves(hdyn, nn, cpt) {
		    real r2 =
			square(bbpos -
			       hdyn_something_relative_to_root(cpt,
							       &hdyn::get_pos));
		    if (r2 < minr2) {
			minr2 = r2;
			nncpt = cpt;
		    }
		}

// 		cerr << "correct_external_nodes: replacing nn "
// 		     << nn->format_label();
// 		cerr << " by component " << nncpt->format_label();
// 		cerr << " for node " << bb->format_label() << endl;

		bb->set_nn(nncpt);
		nn = nncpt;
	    }

	    // Replace coll by nn if a change is needed.

	    hdyn *coll = bb->get_coll();
	    if (coll && coll != nn && coll->is_valid()
		&& coll->is_parent()
		&& is_descendent_of(coll, b, 0)) bb->set_coll(nn);

	    if (bb->is_parent() && bb->get_n_perturbers() > 0) {

		bool flag = false;
		int np = bb->get_n_perturbers();
		if (np > MAX_PERTURBERS) np = MAX_PERTURBERS;

		for (int i = 0; i < np; i++) {
		    hdyn *pert = bb->get_perturber_list()[i];

		    if (pert && pert->is_valid()) {

			if (is_descendent_of(pert, b, 1)) flag = true;

			if (pert->is_parent()
			    && is_descendent_of(pert, b, 0)) {

			    // Split pert into components (if
			    // possible: if overflow would occur,
			    // don't expand the list -- fix this
			    // soon!).

// 			    cerr << "integrate_multiple: expanding perturber "
// 				 << pert->format_label() << " into components";
// 			    cerr << " for node " << bb->format_label() << endl;

			    int np1 = bb->get_n_perturbers();   // np unchanged
			    int j = 0;
			    for_all_leaves(hdyn, pert, pp) {
				if (j == 0)
				    bb->get_perturber_list()[i] = pp;
				else
				    if (np1 < MAX_PERTURBERS)
				        bb->get_perturber_list()[np1++] = pp;
			    }
			    bb->set_n_perturbers(np1);
			}
		    }
		}
		if (flag) nlist[ilist++] = bb;
	    }
	}
    }
    return ilist;
}

local void flatten_hierarchy(hdyn *b)
{
    // b is the top-level node of the clump in question.  Begin by
    // flattening it to a single-level tree.

    // By construction, there should not be any kepler structures attached
    // to any of the nodes under b.  Make sure this is the case.

    for_all_nodes(hdyn, b, bb)
	if (bb->get_kepler())
	    if (bb->get_elder_sister() == NULL) {
		delete bb->get_kepler();
		bb->set_kepler(NULL);
		bb->get_younger_sister()->set_kepler(NULL);
	    }

    int n_mult = b->flatten_node();
//    cerr << "integrate_multiple: "; PRC(n_mult); PRL(b->n_daughters());
//    cerr << "flattened tree to a single level" << endl;
}

local real get_diameter(hdyn *b)
{
    // Return the diameter of the *flat* tree under b.

    real diameter2 = 0;
    for_all_daughters(hdyn, b, bi)
	for (hdyn *bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister()) {
	    real rij2 = abs(bi->get_pos() - bj->get_pos());
	    if (rij2 > diameter2) diameter2 = rij2;
	}
    return sqrt(diameter2);
}

local real perturber_pot(hdyn *b)
{
    // Return the potential on CM b due to its perturbers.

    hdyn *od = b->get_oldest_daughter();
    if (!od) return 0;

    int np = b->get_n_perturbers();
    if (np <= 0) return 0;

    int n = 0;
    for_all_leaves(hdyn, b, bb) n++;

    real masslist[n];
    vec poslist[n];
    n = 0;
    for_all_leaves(hdyn, b, bb) {
	masslist[n] = bb->get_mass();
	poslist[n] = hdyn_something_relative_to_root(bb, &hdyn::get_pos);
	n++;
    }

    // PRL(np);
    if (np > MAX_PERTURBERS) np = MAX_PERTURBERS;

    real phi = 0;
    for (int k = 0; k < np; k++) {
	hdyn *p = b->get_perturber_list()[k];
	if (p && p->is_valid()) {
	    vec ppos = hdyn_something_relative_to_root(p, &hdyn::get_pos);
	    real dphi = 0;
	    for (int i = 0; i < n; i++)
		dphi -= masslist[i]/abs(ppos-poslist[i]);
	    phi += p->get_mass()*dphi;
	}
    }

    return phi;
}

#define DIAG 0

real integrate_multiple(hdyn *b, int mode)	// default mode = 1
{
    // Integrate the internal motion of the multiple system lying
    // under node b, stopping when the interaction is "over," in some
    // sense that depends on the specific system under study.  Here
    // the assumption is that the system is a hard hierarchical one in
    // which the innermost binary has just become perturbed (mode =
    // 1), or where unperturbed motion in a very eccentric system is
    // prevented from becoming unperturbed by the presence of a
    // neighbor (mode = 2).  We end when unperturbed hierarchical
    // motion is restored, or when the problematic close encounter is
    // over.

    //-----------------------------------------------------------------
    // TO DO (SLWM, 05/07):
    //
    // 1. Modify the termination condition in mode 2 to allow for the
    //    (unlikely) possibility that the clump contains another close
    //    pair that might be left in an undesirable state once the
    //    original close encounter is over.
    //
    // 2. Absorb the tidal error into the relative velocity of the
    //    top-level pair.
    //
    // 3. Allow for the possibiliy of collisions, and transmit that
    //    outcome back to the calling function.
    //
    //-----------------------------------------------------------------

    // Determine the potential due to perturbers, for use in computing
    // the tidal error at end.

    real phi0 = perturber_pot(b);

    // Determine some properties of the multiple, for use in setting the
    // output energy scale and termination conditions in smallN_evolve.

    hdyn *inner;				// inner CM
    kepler k_in = get_multiple_params(b, inner);

    real mass = k_in.get_total_mass();
    real period = k_in.get_period();
    real sma = k_in.get_semi_major_axis();

    cerr << "in integrate_multiple, "; PRL(mode);
    if (DIAG) {
	PRC(mass); PRC(sma); PRL(inner->format_label());
    }

    // Inner is the center of mass of the inner binary.

    hdyn *sis = inner->get_binary_sister();

    // Flag and ignore incorrect call (Steve, 2.08).

    if (sis->is_top_level_node()) {		// shouldn't happen
        cerr << "integrate_multiple: binary system!" << endl;
	return 0;
    }
    kepler k_out;
    initialize_kepler_from_dyn_pair(k_out, inner, sis, false);
    real rnn = k_out.get_separation();

    // Determine convenient energy and length scales for the integration.

    real energy = k_in.get_energy();
    bool bound = (energy < 0);
    real pot = -k_in.get_total_mass()/k_in.get_separation();

    smallN_energy_scale = -min(energy, pot/2)/100;
    smallN_time_scale = period;

    // Correct all nodes potentially affected by the multiple
    // integration -- those having an nn, coll, or perturber which is
    // a node (not a leaf, since those are preserved) under b.  Also
    // create a list (nlist) of nodes elsewhere in the system with a
    // perturber within b.

    int ntot = b->get_root()->n_daughters();
    hdynptr *nlist = new hdynptr[ntot];

    int ilist = correct_external_nodes(b, nlist);

    //-----------------------------------------------------------------
    //
    // Remove all hierarchical structure within b.  Note that all
    // center of mass nodes will vanish, so the integration list must
    // check for invalid nodes.  Also, from here until we restore the
    // hierarchy below, the data in b are INCOMPATIBLE with the rest
    // of kira.

    flatten_hierarchy(b);

    real diameter = get_diameter(b);
//    PRC(smallN_energy_scale); PRC(period); PRL(bound);
//    PRC(rnn); PRL(diameter);

    // Save, set, and restore at the end the root node time and system time.
    // Appears that the root node really has to have a NULL parent...
    // *Might* be necessary to suppress GRAPE use too (check).

    xreal save_system_time = b->get_system_time();
    hdyn *save_root = b->get_root();
    hdyn *save_parent = b->get_parent();

    b->set_system_time(0);
    b->set_root(b);
    b->set_parent(NULL);

    b->set_time(0);
    set_all_times(b);

    char save_name[1024];
    strcpy(save_name, b->get_name());
    b->set_name("top");

    // Move to the center of mass frame.

    vec cmpos, cmvel;
    b->set_system_time(-1);
    compute_com(b, cmpos, cmvel);
    b->inc_pos(-cmpos);
    b->inc_vel(-cmvel);
    // cerr << "moved to CM frame" << endl;

    b->set_system_time(0);

    // Set integration parameters before calling smallN_evolve.  For a
    // sufficiently bound innermost system, the preferred termination
    // condition is unperturbed approaching motion with a receding
    // nearest neighbor.

    // Original use (mode = 1) of this function was to get a close
    // multiple around the pericenter of the outer orbit, during which
    // the inner binary is perturbed and the limited dynamic range
    // makes the time step insignificant.  In that case, the
    // integration should proceed for up to the outer orbital period,
    // but will preferably stop when the innermost binary is again
    // unperturbed.

    // New use (mode = 2, 5/07) is more modest, and simply integrates
    // the inner binary past its pericenter.  A typical situation here
    // is that the inner binary is not necessarily tight, but is very
    // eccentric, and the time step becomes insignificant before the
    // unperturbed criterion kicks in.

    // Probably a more unified termination condition should be
    // developed...

    real r_end = 2*rnn, t_end = VERY_LARGE_NUMBER, t_log = VERY_LARGE_NUMBER;
    bool end_on_unpert = false;

    if (mode == 1) {

	if (bound) end_on_unpert = true;

	// Set a fallback termination time.

	real t = k_out.get_time();
	if (k_out.get_rel_pos()*k_out.get_rel_vel() > 0) {

	    // Outer system is separating.  Integrate to apastron or for
	    // a time equal to the time since periastron..

	    real t_peri = t - k_out.get_time_of_periastron_passage();
	    real t_apo = VERY_LARGE_NUMBER;
	    if (k_out.get_energy() < 0)
	        t_apo = k_out.pred_advance_to_apastron() - t;

	    t_end = min(t_peri, t_apo);

	} else {

	    // Outer system is approaching.  Integrate to the reflection
	    // point.

	    t_end = 2*(k_out.pred_advance_to_periastron() - t);

	}

	t_log = 10*period;

    } else {

        // Just want to integrate the inner binary to its reflection
        // point.  We could also (maybe only?) require that the
        // minimum interparticle separation at termination exceed the
        // initial separation of k_in (would satisfy TO DO #1 above).
        // This might allow us to devote less effort to computing the
        // kepler structures.

	real t = k_in.get_time();
	if (k_in.get_rel_pos()*k_in.get_rel_vel() > 0) {

	    // Inner system is separating.  Not clear why this
	    // function was called...  Integrate to apastron or for a
	    // time equal to the time since periastron..

	    real t_peri = t - k_in.get_time_of_periastron_passage();
	    real t_apo = VERY_LARGE_NUMBER;
	    if (k_in.get_energy() < 0)
	        t_apo = k_in.pred_advance_to_apastron() - t;

	    t_end = min(t_peri, t_apo);

	} else {

	    // Inner system is approaching.  Integrate to the reflection
	    // point.

	    t_end = 2*(k_in.pred_advance_to_periastron() - t);

	}

	t_log = t_end;

    }

    // PRC(t_end); PRC(r_end); PRL(end_on_unpert);

    // Must also set dt_crit and r_crit appropriately.  For r_crit,
    // use twice the semi-major axis of the inner binary.  For
    // dt_crit, use the kepler step at r_crit.
    //
    // May need fine-tuning...

    real rcrit = 2*min(k_in.get_separation(), sma);
    real dtcrit = sqrt(eta*kepler_step_sq(rcrit, mass,
					  2*(energy + mass/rcrit)));

    // PRC(eta); PRC(rcrit); PRC(mass); PRL(energy);
    // PRL(energy + mass/rcrit);
    // PRL(kepler_step_sq(rcrit, mass, 2*(energy + mass/rcrit)));
    
    set_smallN_dtcrit(dtcrit);
    set_smallN_rcrit(rcrit);

    // PRC(get_smallN_dtcrit()); PRL(get_smallN_rcrit());
    local_sys_stats = true;

    //=================================================================
    //
    // Integrate the clump.

    cerr << "entering smallN_evolve()" << endl << flush;
    int status = smallN_evolve(b, t_end, r_end, end_on_unpert, DIAG*t_log);
    cerr << "exiting smallN_evolve(), ";
    real t = b->get_time();
    PRL(t); cerr << flush;

    //=================================================================

    // Restore saved global data.

    b->set_root(save_root);
    b->set_parent(save_parent);
    b->set_system_time(save_system_time);
    b->set_time(save_system_time);
    set_all_times(b);

    // Set an appropriate time step that fits the block structure.

    real dt = b->get_timestep();
    real dt2 = 1;
    while (dt2 > dt) dt2 /= 2;
    while (fmod(save_system_time, dt2) != 0) dt2 /= 2;
    set_all_timesteps(b, dt2);

    // Restore binary tree structure.

    restore_binary_tree(b);
    // b->set_name(save_name);
    label_binary_node(b);	// reflect new structure, in case it changed
    PRL(b->format_label());

    // Restore the proper center of mass position and velocity.

    b->inc_pos(cmpos);
    b->inc_vel(cmvel);

    // Check the tidal error.  Currently we only monitor the error,
    // applying no correction -- should do so as in hdyn_unpert.C.
    // See TO DO #2 above.

    real phi1 = perturber_pot(b);

    int np = b->get_n_perturbers();
    PRC(np); PRC(phi0); PRL(phi1); PRL(phi1-phi0);

    if (status != 3) {

	// See if the closest pair is unperturbed (wasn't the
	// termination condition.

	if (fully_unperturbed(b)) status = 3;
    }

    PRL(status);

    if (status == 3) {

	// Officially restart unperturbed motion, and replace bi_min
	// and bj_min by CM on perturber lists listed in nlist.

	cerr << "startup_unperturbed_motion for " << bi_min->format_label()
	     << endl << flush;
	bi_min->startup_unperturbed_motion();

	for (int i = 0; i < ilist; i++) {
	    hdyn *bb = nlist[i];
	    if (bb) {
		int np = bb->get_n_perturbers();
		if (np > MAX_PERTURBERS) np = MAX_PERTURBERS;
		int count = 0;
		for (int ip = 0; ip < np; ip++) {
		  hdyn *pert = bb->get_perturber_list()[ip];
		  if (pert == bi_min || pert == bj_min) {
		    pert = bi_min->get_parent();
		    count++;
		  }
		  if (count > 1 && ip+count-1 < np) {
		    bb->get_perturber_list()[ip]
		      = bb->get_perturber_list()[ip+count-1];
		    ip--;
		  }
		}
		np -= count-1;
		bb->set_n_perturbers(np);
	    }
	}
    }

    cerr << "leaving integrate_multiple" << endl << flush;

    delete [] nlist;
    return t;
}



//----------------------------------------------------------------------

#else

main(int argc, char *argv[])
{
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "a:d:D:E:g:n:r:t:";

    real dt_log = 0;
    bool d_set = false;
    real dt_snap = 0;
    bool D_set = false;
    real dt_energy = 0;
    bool E_set = false;
    real break_r2 = VERY_LARGE_NUMBER;
    real t_end = 200;

    while ((c = pgetopt(argc, argv, param_string,
			"$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'a': set_smallN_eta(atof(poptarg));
		      break;
	    case 'd': dt_log = atof(poptarg);
	    	      d_set = true;
		      break;
	    case 'D': dt_snap = atof(poptarg);
	    	      D_set = true;
		      break;
	    case 'E': dt_energy = atof(poptarg);
	    	      E_set = true;
		      break;
	    case 'g': set_smallN_gamma(atof(poptarg));
		      break;
	    case 'n': set_smallN_niter(atoi(poptarg));
		      break;
	    case 'r': break_r2 = pow(atof(poptarg),2);
		      break;
	    case 't': t_end = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      exit(1);
        }            

    PRL(get_smallN_eta());
    PRL(get_smallN_gamma());
    PRL(get_smallN_niter());
    PRL(break_r2);
    PRL(t_end);

    hdyn *b = get_hdyn();
    b->log_history(argc, argv);

    b->set_label("root");
    b->set_root(b);			// bad!!  necessary??

    kira_options ko;
    ko.perturber_criterion = 2;
    b->set_kira_options(&ko);

    for_all_nodes(hdyn, b, bi) bi->set_time(0);

    // Use smallN_evolve() as the integrator.  Later, may want to add an
    // option to use the autoscaling features of integrate_multiple.

    cerr.precision(8);
    int status = smallN_evolve(b, t_end, break_r2,
			       false, dt_log, dt_energy, dt_snap);

    if (dt_log == 0 && !d_set) {
	real t_sys = b->get_system_time();
	b->set_system_time(b->get_time());
	sys_stats(b);
	b->set_system_time(t_sys);
    }

    if (dt_energy == 0 && !E_set)
	print_recalculated_energies((dyn*)b, 1);

    if (dt_snap == 0 && !D_set)	put_node(b);
}

#endif
