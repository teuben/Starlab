
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//
//  _dyn_ev.C: functions related to orbit integration within the _dyn_ class.
//    		     - member functions for FLAT trees only
//.............................................................................
//    version 1:  Jul 1996   Steve McMillan
//    version 2:
//.............................................................................

#include "_dyn_.h"

//-----------------------------------------------------------------------------
//  flat_accumulate_acc_and_jerk -- calculates the contribution to the
//                                  acceleration & jerk on this node from
//                                  the node leaf.
//                                  NOTE: for flat trees only
//-----------------------------------------------------------------------------

void _dyn_::flat_accumulate_acc_and_jerk(_dyn_ * leaf,	// attracting leaf node
					real eps2)	// softening length^2
{
    vector d_pos = leaf->get_pred_pos() - get_pred_pos();
    vector d_vel = leaf->get_pred_vel() - get_pred_vel();
    real r2inv = 1.0 / (d_pos * d_pos + eps2);
    real a3 = -3.0 * (d_pos * d_vel) * r2inv;
    real mrinv = leaf->get_mass() * sqrt(r2inv);
    pot -= mrinv;
    real mr3inv = mrinv * r2inv;
    acc += mr3inv * d_pos;
    jerk += mr3inv * (d_vel + a3 * d_pos);
}

//-----------------------------------------------------------------------------
//  flat_calculate_acc_and_jerk -- calculates acceleration & jerk on this node:
//                                 from all nodes under p (if p is root ptr);
//                                 or from the node p (if p is a leaf ptr).
//
//                                 NOTE: for flat trees only
//
//-----------------------------------------------------------------------------

#ifdef RECURSIVE_LOOP

  void _dyn_::flat_calculate_acc_and_jerk(_dyn_ * p,   // root or leaf
					  real eps2)   // softening length^2
  {
    if (p->is_root()) {			// p is root, so
					// invoke for all leaves
	for_all_daughters(_dyn_, p, d)	
	    flat_calculate_acc_and_jerk(d, eps2);

    } else {				// p is leaf, but
					// not this leaf
	if (p != this)
	    flat_accumulate_acc_and_jerk(p, eps2);

    }
  }

#else

  // Looping by recursion (as above) is relatively expensive on many systems.
  // Remove the recursion for the sake of efficiency (and readability?) by
  // using an explicit loop.

  void _dyn_::flat_calculate_acc_and_jerk(_dyn_ * p,   // root or leaf
					  real eps2)   // softening length^2
  {
    if (p->is_root()) {			// p is root, so
					// invoke for all leaves
	for_all_daughters(_dyn_, p, d)	
	    if (d != this)
		flat_accumulate_acc_and_jerk(d, eps2);

    } else {				// NOTE: else clause currently not used

	if (p != this)
	    flat_accumulate_acc_and_jerk(p, eps2);

    }
  }

#endif

#define SAFETY_FACTOR 16

void _dyn_::flat_set_first_timestep(real eta_for_firststep,
				    real max_step_size)
{
    // Don't have high-order derivatives, so do the best we can
    // and include an extra safety factor.

    real a1scaled2 = jerk * jerk;
    real a2 = acc * acc;
    real newstep;

    if (a1scaled2 > 0.0)
	newstep = eta_for_firststep * sqrt(a2 / a1scaled2);
    else
	newstep = eta_for_firststep / SAFETY_FACTOR;

    newstep = adjust_number_to_power(newstep, max_step_size);
    timestep = newstep;
}

local real flat_new_timestep(vector & at3,	// third order term
			     vector & bt2,	// 2nd order term
			     vector & jerk,	// 1st order term
			     vector & acc,	// 0th order term
			     real timestep,	// old timestep
			     real time,		// present time
			     real eta,		// accuracy parameter
			     real dtmax)	// maximum stepsize
{
    real a3scaled2 = 36 * at3 * at3;
    real a2scaled2 = 4 * bt2 * bt2;
    real a1scaled2 = timestep * jerk * jerk;
    real a2 = acc * acc;

    // Simple criterion:

    real newstep = eta * sqrt(sqrt(a2 / a2scaled2)) * timestep;

    // Probably should use the Aarseth criterion:

    // (see ../hdyn/evolve/hdyn_ev.C...)

    // Force step into the block scheme.

    if (newstep < timestep)
	return timestep * 0.5;
    else if (newstep < 2 * timestep)
	return timestep;
    else if (fmod(time, 2 * timestep) != 0.0 || 2 * timestep > dtmax)
	return timestep;
    else
	return timestep * 2;
}

void _dyn_::flat_update(const real eta, const real dtmax)
{
    vector at3 = 2 * (old_acc - acc) + timestep * (old_jerk + jerk);
    vector bt2 = -3 * (old_acc - acc) - timestep * (2 * old_jerk + jerk);
    time = time + timestep;
    timestep = flat_new_timestep(at3, bt2, jerk, acc,
				 timestep, time, eta, dtmax);
}

//-----------------------------------------------------------------------------
// flat_correct:  Apply the corrector for 4th-order hermite scheme.
//                See Makino & Aarseth (PASJ 1992) for derivation.
//-----------------------------------------------------------------------------

#define ONE12 0.0833333333333333333333
#define ONE3  0.3333333333333333333333

bool _dyn_::flat_correct()
{
    vector at3 = 2 * (old_acc - acc) + timestep * (old_jerk + jerk);
    vector bt2 = -3 * (old_acc - acc) - timestep * (2 * old_jerk + jerk);
    real dt2 = timestep * timestep;

    // Reminder:	at3  =  a''' dt^3 / 6		(a' = j)
    //			bt2  =  a''  dt^2 / 2

    vector new_pos = pred_pos + (0.05 * at3 + ONE12 * bt2) * dt2;
    vector new_vel = pred_vel + (0.25 * at3 + ONE3 * bt2) * timestep;

    pos = new_pos;
    vel = new_vel;

    return true;		// normal return value
}

//-----------------------------------------------------------------------------
// predict_loworder_all:  Predict the positions of all nodes below b at time t.
//-----------------------------------------------------------------------------

//#define FIFTH_ORDER_PRED		// doesn't work yet...

void predict_loworder_all(_dyn_ * b,
			  xreal t)
{
    if (!b) return;			// shouldn't happen -- flag as error?

    // Basic recursion:

    if (b->is_parent())
	for_all_daughters(_dyn_, b, d)
	    predict_loworder_all(d, t);

    if (b->is_top_level_node()) {

#ifdef FIFTH_ORDER_PRED
	real dt = b->get_timestep();
	vector kdt = 18*dt*b->get_k_over_18();
	PRL(b->get_jerk());
	PRC(dt); PRL(kdt);
	b->predict_loworder5(t);	// new -- for GRAPE-6 kira...
					// still being tested...
#else
	b->predict_loworder(t);
#endif

    } else {

	// For low-level nodes, the secondary pred_pos and pred_vel are set
	// once the primary is known.  May as well take advantage of that...

	_dyn_ *s = b->get_elder_sister();
	if (s && s->get_t_pred() == t)

	    // This is a younger binary sister.  Use the elder sister
	    // data to compute pred_pos and pred_vel.

	    b->predict_from_elder_sister();

	else

	    // This is an elder binary sister, or the elder sister is
	    // somehow not up to date.

	    b->predict_loworder(t);

    }
}





