

       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                       
//=======================================================//              /|\ ~


//// Self-contained few-body integrator, using a fourth-order Hermite scheme
//// incorporating a modified unperturbed treatment of close approaches and
//// time symmetrization. This version is intended ultimately to be incorporated
//// into kira to handle multiple systems containing perturbed hard binaries.
////
//// Usage: kira_smallN [OPTIONS]
////
//// Options:
////	     -a    set accuracy parameter [0.03]
////	     -g    set unperturbed limit [1.e-5]
////	     -n    set number of symmetrization iterations [1]
////	     -r    set termination radius [infinite]
////	     -t    set termination time [200]
////
//// Written by Steve McMillan and Fan-Chi Lin.
////
//// Report bugs to starlab@sns.ias.edu.

// As currently envisaged, on entry we will have an isolated (that is,
// relatively unperturbed) few-body (N <~ 5) system in which a hard
// binary has just become perturbed by a neighbor.  A typical scanario
// might be that a stable triple has been perturbed into an unstable
// orbit, and now the outer component has become a perturber of the
// inner pair.  Fow now, we
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
//	  structure will cause the multiple to be integrated to the
//	  specified time.  Correction will have no additional effect.
//	  Need to check that the normally flat but noe binary tree used
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
typedef	  hdyn DYN;
#define   GET get_hdyn

// Make these arguments global later...
// Note that eta and gamma2_unpert can be set from the command line.

real eta = 0.03;		// empirical
real gamma2_unpert = 1.e-10;	// conservative, given first-order correction

const bool print = false;
real t_snap = 0;
const real dt_snap = 0.05;

const real dt_crit = 1.e-4;
const real r2_crit = 1.e-6;

DYN *bi_min = NULL, *bj_min = NULL, *bk_nn = NULL;		// temporary

bool dbg_allow_binary = true;	// for debugging, set false to disable binaries
real rbin = 0;			// use to keep track of close pair



local inline real get_pairwise_acc_and_jerk(DYN *bi, DYN *bj,
					    vec &force, vec &jerk)
{
    // Compute the force and jerk on node bi due to node bj by summing
    // all pairwise component forces.  Return the unscaled time step
    // appropriate to the minimum distance between any component of bi
    // and any component of bj.

    force = jerk = 0;
    real min_distance2 = VERY_LARGE_NUMBER, min_distance3;

    vec delx = 0, delv = 0;
    if (bi->is_parent()) {		// NOTE: assuming for now that our trees
	delx -= bi->get_nopred_pos();	//	 are just a single level deep
	delv -= bi->get_nopred_vel();
    }
    if (bj->is_parent()) {
	delx += bj->get_nopred_pos();
	delv += bj->get_nopred_vel();
    }

    for_all_leaves(DYN, bi, bbi) {

	// Compute the force and jerk on leaf bbi due to all leaves under bj.

	vec iforce = 0, ijerk = 0;
	for_all_leaves(DYN, bj, bbj) {

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

    // Alternate time step criterion:

    // real timestep2 = min_distance3 / (bj->get_mass() + bi->get_mass());

    real dtff2 = 0.5 * min_distance3 / (bj->get_mass() + bi->get_mass());
    real dtv2  = min_distance2
			/ square(bj->get_nopred_vel() - bi->get_nopred_vel());

    return Starlab::min(dtff2, dtv2);
}



local inline real calculate_top_level_acc_and_jerk(DYN *b)
{
    // Compute the acc and jerk on all top-level nodes.  We could use the
    // hdyn functions in hdyn_ev.C to do this more simply, but they will
    // double-count all pairwise interactions.  All nodes must be resolved
    // into components for purposes of computing the acc and jerk.  Return
    // the minimum time step associated with any top-level pair.

    for_all_daughters(DYN, b, bi) {
	bi->set_acc(0);
	bi->set_jerk(0);
    }

    real min_timestep2 = VERY_LARGE_NUMBER;
    bi_min = bj_min = NULL;

    for_all_daughters(DYN, b, bi) {
	for (DYN *bj = bi->get_younger_sister();
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

local inline void advance_components_to_time(DYN *bi, real t)	  // unperturbed
{								  // "predictor"
    DYN *od = bi->get_oldest_daughter();

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

	    DYN *yd = od->get_younger_sister();
	    yd->set_pred_pos((1-fac)*k->get_rel_pos());
	    yd->set_pred_vel((1-fac)*k->get_rel_vel());
	}
    }
}

local inline void update_components_from_pred(DYN *bi)		  // unperturbed
{								  // "corrector"
    // Called after a step is completed.

    DYN *od = bi->get_oldest_daughter();

    if (od) {	// redundant

	od->set_time(bi->get_time());
	od->set_pos(od->get_nopred_pos());
	od->set_vel(od->get_nopred_vel());
	DYN *yd = od->get_younger_sister();
	yd->set_time(bi->get_time());
	yd->set_pos(yd->get_nopred_pos());
	yd->set_vel(yd->get_nopred_vel());
    }
}

local inline void correct_timestep_for_components(DYN *bi, real t, real &dt)
{
    // Reduce the time step if necessary to end at the termination
    // time for unperturbed motion.

    DYN *od = bi->get_oldest_daughter();

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

local inline bool unperturbed_and_approaching(DYN *b,DYN *bi, DYN *bj)
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
    for_all_leaves(DYN, b, bk)
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



local inline bool create_binary(DYN *bi, DYN *bj)
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

    DYN *cm = new DYN;
    cm->set_mass(mi+mj);
    cm->set_pos((mi*bi->get_pos()+mj*bj->get_pos())/cm->get_mass());
    cm->set_vel((mi*bi->get_vel()+mj*bj->get_vel())/cm->get_mass());
    cm->set_time(time);

    // Compute the tidal potential (components - CM) due to the rest
    // of the system.  Store it in bi->pot.

    real dphi_tidal = 0;
    for_all_daughters(DYN, bi->get_parent(), bk)
	if (bk != bi && bk != bj) {
	    real dphi_k = -bi->get_mass()/abs(bi->get_pos()-bk->get_pos())
			  -bj->get_mass()/abs(bj->get_pos()-bk->get_pos())
			  +cm->get_mass()/abs(cm->get_pos()-bk->get_pos());
	    dphi_tidal += bk->get_mass()*dphi_k;
	}
    bi->set_pot(dphi_tidal);

    if (!dbg_allow_binary) {

	// Don't actually create the binary.  Print info and return.

	bj->set_pot(dphi_tidal);
	PRL(bi->format_label());
	PRL(bj->format_label());
	PRL(dphi_tidal);
	print_binary_from_dyn_pair(bi, bj, 0, 0, true, true);
	return true;
    }

    // Offset components to the center of mass frame.

    bi->inc_pos(-cm->get_pos());
    bi->inc_vel(-cm->get_vel());
    bj->inc_pos(-cm->get_pos());
    bj->inc_vel(-cm->get_vel());

    cerr << endl << "too close: " << bi->format_label();
    cerr << " and " << bj->format_label() << endl;

    DYN *p = bi->get_parent();
    pp(p), cerr << endl;

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

    cerr << "created new binary " << cm->format_label() << endl;
    pp(p), cerr << endl;
    k->print_all();
    real dt_unpert = bi->get_t_pred() - time;
    PRC(dt_unpert); PRL(bi->get_t_pred());

    // Make sure pos and pred_pos agree.

    cm->set_t_pred(time);
    for_all_nodes(DYN, cm, bb) {
	bb->set_pred_pos(bb->get_pos());
	bb->set_pred_vel(bb->get_vel());
    }

    return true;
}



local inline bool terminate_binary(DYN*& bi)
{
    if (!bi) return false;
    DYN *od = bi->get_oldest_daughter();
    if (!od) return false;
    DYN *yd = od->get_younger_sister();
    if (!yd) return false;

    kepler *k = od->get_kepler();
    if (!k) err_exit("smallN: binary with no kepler.");

    cerr << endl << "terminating binary " << bi->format_label() << endl;
    k->print_all();
    PRL(od->get_t_pred());

    DYN *p = bi->get_parent();
    pp(p), cerr << endl;

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
    for_all_daughters(DYN, bi->get_parent(), bk)
	if (bk != bi) {
	    real dphi_k = -od->get_mass()/abs(od->get_pos()-bk->get_pos())
			  -yd->get_mass()/abs(yd->get_pos()-bk->get_pos())
			  +bi->get_mass()/abs(bi->get_pos()-bk->get_pos());
	    dphi_tidal += bk->get_mass()*dphi_k;
	}
    PRC(dphi_tidal); PRC(od->get_pot());
    dphi_tidal -= od->get_pot();
    PRL(dphi_tidal);

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

    cerr << "corrected component velocities: "; PRL(vfac);

    // Update the tree.

    add_node_before(od, bi);
    add_node_before(yd, bi);
    detach_node_from_general_tree(bi);

    // Delete the CM node (careful to detach the daughters to
    // avoid recursive deletion!).

    bi->set_oldest_daughter(NULL);
    delete(bi);

    pp(p), cerr << endl;

    // Make sure pos and pred_pos agree.  Note that init_pred sets
    // t_pred = time as well as pred_pos = pos.

    od->init_pred();
    yd->init_pred();

    // Change bi to yd, so next node is correct in loop below.

    bi = yd;
    return true;
}



local inline real acc_and_jerk_and_get_dt(DYN *b)
{
    // Compute acc and jerk, and return the top-level time step.

    real dt = calculate_top_level_acc_and_jerk(b);

    // Adjust dt to take unperturbed systems into account.

    real t = b->get_time();
    for_all_daughters(DYN, b, bi)
	if (bi->is_parent()) correct_timestep_for_components(bi, t, dt);

    // Don't store the time steps yet, as the old time step may
    // be needed for the corrector.

    return dt;
}

local inline void set_all_timesteps(DYN *b, real dt)
{
    // Update all time steps (root and top-level nodes).

    for_all_daughters(DYN, b, bi) bi->set_timestep(dt);
    b->set_timestep(dt);
}

local void set_all_times(DYN *b)
{
    // Update times of all nodes in the system from system_time.

    real sys_t = b->get_time();
    for_all_nodes(DYN, b, bi) bi->set_time(sys_t);
}



// Routines transplanted from sdyn/evolve/sdyn_ev.C.  Naming scheme is
// that pred_pos here is new_pos there (similarly vel), old_acc here is
// acc there, acc here is new_acc there (similarly jerk).

void correct_acc_and_jerk(DYN *bi,
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

void correct_pos_and_vel(DYN *bi, const real new_dt)
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



local void print_pred_energies(DYN *b)
{
    // Compute and print energies based on "pred" quantities.

    int n = 0;
    for_all_nodes(DYN, b, bb) n++;
    vec *p = new vec[n];
    vec *v = new vec[n];

    n = 0;
    for_all_nodes(DYN, b, bb) {
	p[n] = bb->get_pos();
	bb->set_pos(bb->get_nopred_pos());
	v[n] = bb->get_vel();
	bb->set_vel(bb->get_nopred_vel());
	n++;
    }

    int pr = cerr.precision(HIGH_PRECISION);
    cerr << b->get_time() << " (pred): ";
    cerr.precision(pr);
    print_recalculated_energies((dyn*)b, 1, 0);

    n = 0;
    for_all_nodes(DYN, b, bb) {
	bb->set_pos(p[n]);
	bb->set_vel(v[n]);
	n++;
    }
}

// Take a step to time t.  Return the actual system time at the end
// of the step, after symmetrization if specified.

local real take_a_step(DYN *b,		// root node
		       real &dt,	// natural step at start/end
		       int n_iter = 1)	// number of iterations to symmetrize
					// 0 == > explicit, > 0 ==> implicit
					// default of 1 should be adequate
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

    for_all_daughters(DYN, b, bi) {
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
	// apply the corrector for this iteraation.

	for_all_daughters(DYN, b, bi) {
	    if (new_dt != prev_new_dt)
		correct_acc_and_jerk(bi, new_dt, prev_new_dt);
	    correct_pos_and_vel(bi, new_dt);
	    if (bi->is_parent()) advance_components_to_time(bi, t+new_dt);
	}
     }

    // Complete the step.

    for_all_daughters(DYN, b, bi) {
	bi->set_pos(bi->get_pred_pos());
	bi->set_vel(bi->get_pred_vel());
	if (bi->is_parent())
	    update_components_from_pred(bi);
    }

    dt = end_point_dt;
    return t + new_dt;
}



bool is_fully_unperturbed(DYN *bi, DYN *bk)
{
    // Check if the binary whose elder component is bi is perturbed by bk.
    // Use the same criteria as in the rest of kira.  To come...

    return false;
}

// Evolve the system to time t_end.  Stop if any particle gets too far
// from the origin.

real smallN_evolve(DYN *b,
		   int n_iter = 1,
		   real t_end = 200,
		   real break_r2 = VERY_LARGE_NUMBER)
{
    b->init_pred();
    for_all_daughters(DYN, b, bi) bi->init_pred();

    int n_steps = 0;

    real dt = acc_and_jerk_and_get_dt(b);	// initialization
    set_all_timesteps(b, dt);
    set_all_times(b);

    // Note that we correctly maintain all times and time steps, but
    // we do not modify system_time.  It is not strictly necessary to
    // update all times, as times and time steps are shared, but it is
    // convenient to keep the system coherent for output purposes.

    // Optional output:

    if (print) {
	t_snap += dt_snap;
	real t_sys = b->get_system_time();
	b->set_system_time(b->get_time());
	put_node(b);
	sys_stats(b, 1, true, true, false, 2, true, true, true);
	b->set_system_time(t_sys);
    }

    print_recalculated_energies((dyn*)b, 0, 0);

    while (b->get_time() < t_end) {

	// Check termination criterion:

	for_all_daughters(DYN, b, bi)
	    if (square(bi->get_pos()) > break_r2) return t_end + 1;

	// Take a step.  Don't set time in advance, as the symmetrization
	// process will determine the actual time step.  Thus, during the
	// entire step, all times refer to the time at the *start* of the
	// step.  The time step dt on entry is the natural time step for
	// the system at the current time.  On return it is the new natural
	// time step for the system at the end of the step.  The return
	// value of the function is the new system time.

	b->set_time(take_a_step(b, dt, n_iter));

	// Update all times and time steps.

	set_all_timesteps(b, dt);
	set_all_times(b);
	n_steps++;

#if 0
	cerr << b->get_time() << " (" << n_steps << "): ";
	print_recalculated_energies((dyn*)b, 1, 0);
#endif

	// The step is over and all times have been advanced.
	// Check for the start/end of unperturbed motion.

	bool tree_changed = false;

	// Time step dt was set by bi_min and bj_min during the last acc
	// and jerk calculation.

	if (rbin > 0 && abs(bi_min->get_pos()-bj_min->get_pos()) > rbin) {

	    // Print debugging info on the binary we didn't actually create.

	    DYN * od = bi_min;
	    DYN * yd = bj_min;
	    PRL(od->format_label());
	    PRL(yd->format_label());
	    real mbin = od->get_mass() + yd->get_mass();
    	    vec cmpos = (od->get_mass()*od->get_pos()
			 + yd->get_mass()*yd->get_pos())/mbin;
	    real dphi_tidal = 0;
	    for_all_daughters(DYN, od->get_parent(), bk)
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

	if (dt < dt_crit && unperturbed_and_approaching(b, bi_min, bj_min)) {

	    tree_changed = create_binary(bi_min, bj_min);
	    // PRL(bk_nn->format_label());

	    if (!dbg_allow_binary && tree_changed) {

		// Didn't actually create a binary.  Save the separation
		// for later use.

		rbin = abs(bi_min->get_pos()-bj_min->get_pos());
		tree_changed = false;
	    }

	    if (tree_changed) {

		// See if the new binary triggers the termination criterion.
		// Nearest neighbor is bk_nn.

		DYN *cm = bi_min->get_parent();
		kepler *k = bi_min->get_kepler();
		PRC(cm); PRL(k);

		if (k->get_energy() < 0) {

		    real gamma = estimated_perturbation(bi_min->get_parent(),
							bk_nn);
		    PRL(gamma);

		    // Conservatively advance the perturbation to apastron.

		    gamma *= pow(k->get_apastron()/k->get_separation(), 3);
		    PRL(gamma);
		    vec dr = hdyn_something_relative_to_root(cm, &hdyn::get_pos)
		       - hdyn_something_relative_to_root(bk_nn, &hdyn::get_pos);
		    PRL(dr);
		    vec dv = hdyn_something_relative_to_root(cm, &hdyn::get_vel)
		       - hdyn_something_relative_to_root(bk_nn, &hdyn::get_vel);
		    PRC(dv); PRL(dr*dv);

		    if (gamma < 1.e-4 && dr*dv > 0) {
			cerr << "new binary is fully unperturbed -- exiting"
			    << endl;
			exit(1);
		    }
		}
	    }
	}

	// Check for termination of unperturbed motion.

	for_all_daughters(DYN, b, bi) {
	    DYN *od = NULL;
	    if (od = bi->get_oldest_daughter())
		if (od->get_t_pred() <= b->get_time()) {
		    int p = cerr.precision(HIGH_PRECISION);
		    cerr << b->get_time() << " (xxxxx): ";
		    cerr.precision(p);
		    print_recalculated_energies((dyn*)b, 1, 0);
		    tree_changed |= terminate_binary(bi);
		}
	}

	if (tree_changed) {

	    // Recompute accs, jerks, and the time step.

	    dt = acc_and_jerk_and_get_dt(b);

	    cerr << endl << "new structure: "; pp(b); cerr << ",  ";
	    PRL(dt);

	    set_all_timesteps(b, dt);
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << b->get_time() << " (yyyyy): ";
	    cerr.precision(p);
	    print_recalculated_energies((dyn*)b, 1, 0);
	}

	// Optional output:

	if (print && b->get_time() > t_snap) {
	    real t_sys = b->get_system_time();
	    b->set_system_time(b->get_time());
	    put_node(b);
	    t_snap += dt_snap;
	    // sys_stats(b, 1, true, true, false, 2, true, true, true);
	    b->set_system_time(t_sys);
	}

	// Debugging: print out the energy.

	int p = cerr.precision(HIGH_PRECISION);
	cerr << b->get_time() << " (" << n_steps << "): ";
	cerr.precision(p);
	print_recalculated_energies((dyn*)b, 1, 0);
    }

    return t_end;
}

local void restore_binary_tree(DYN *b)
{
    // Reinstate the binary tree structure below bode b.
    // Use a simple distance/potential criterion.

    // To come...  Needed for the kira interface.
}

real integrate_multiple(DYN *b)
{
    // Only termination condition for now is unperturbed motion
    // with a receding nearest neighbor.

    // b is the top-level node of the clunmp in question.  Begin
    // by flattening everything to a single-level tree.

    int n_mult = b->flatten_node();
    cerr << "integrate_multiple: "; PRC(n_mult); PRL(b->n_daughters());

    real t = smallN_evolve(b);
    restore_binary_tree(b);

    return t;
}


//----------------------------------------------------------------------

main(int argc, char *argv[])
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:g:n:r:t:";

    int n_iter = 1;
    real t_end = 100; // VERY_LARGE_NUMBER;	// 100;
    real break_r2 = VERY_LARGE_NUMBER;	// 400;

    while ((c = pgetopt(argc, argv, param_string,
			"$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'a': eta = atof(poptarg);
		      break;
	    case 'g': gamma2_unpert = pow(atof(poptarg), 2);
		      break;
	    case 'n': n_iter = atoi(poptarg);
		      break;
	    case 'r': break_r2 = pow(atof(poptarg),2);
		      break;
	    case 't': t_end = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      exit(1);
        }            

    PRC(eta); PRC(gamma2_unpert); PRL(n_iter);
    PRC(break_r2); PRL(t_end);

    DYN *b = GET();
    b->set_label("root");
    b->set_root(b);			// bad!!

    cerr.precision(10);

    kira_options ko;
    ko.perturber_criterion = 2;
    b->set_kira_options(&ko);

    b->set_time(0);
    set_all_times(b);

    // real t = smallN_evolve(b, n_iter, t_end, break_r2);
    real t = integrate_multiple(b);

    sys_stats(b, 1, true, true, false, 2, true, true, true);

    if (t != t_end)
	cerr << "Interaction over!" << endl;
    else
	cerr << "Interaction not over!" << endl;

    put_node(b);
}
