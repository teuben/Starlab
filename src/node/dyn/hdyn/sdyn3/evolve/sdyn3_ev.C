
// sdyn3_ev.C

// Member functions for the sdyn3 class that are related to orbit integration.
// Also include some useful general functions (move elsewhere?...)

#include "sdyn3.h"

real potential_energy(sdyn3 * b)
{
    real u = 0;

    for_all_daughters(sdyn3, b, bi) {
	for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	  u -= bi->get_mass() * bj->get_mass()
	       / abs(bi->get_pos() - bj->get_pos());
    }

    return u;
}

real energy(sdyn3 * b)
{
    real k = 0, u = 0, dissipation = 0;

    for_all_daughters(sdyn3, b, bi) {

	// PRC(bi), PRL(bi->get_pos());

	k += bi->get_mass() * bi->get_vel() * bi->get_vel();
	dissipation += bi->get_energy_dissipation();

	for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister()) {

	    // PRC(bj), PRL(bj->get_pos());

	    u -= bi->get_mass() * bj->get_mass()
		  / abs(bi->get_pos() - bj->get_pos());
	}
    }

    return 0.5*k + u + dissipation;
}

real total_angular_momentum(sdyn3* b)
{
    vector am = vector(0);
    for_all_daughters(sdyn3, b, bb)
	am += bb->get_mass() * (bb->get_pos() ^ bb->get_vel());
    return abs(am);
}

void sdyn3::accumulate_new_acc_and_jerk_from_new(
			sdyn3* bj,                  // n-body system pointer
			real eps2,                  // softening length squared
			int no_ssd_flag,
			int& collision_flag)
{
    if (bj->get_oldest_daughter() != NULL)
        for_all_daughters(sdyn3, bj, bb)
	    accumulate_new_acc_and_jerk_from_new(bb, eps2,
						 no_ssd_flag, collision_flag);
    else
	if (this != bj) {

	    vector d_pos = new_pos - bj->get_new_pos();
	    vector d_vel = new_vel - bj->get_new_vel();
	    real r2 = d_pos*d_pos;

	    if (r2 < (radius + bj->get_radius()) * (radius + bj->get_radius()))
		collision_flag = 1;

	    real r2inv = 1.0/(r2 + eps2);
	    real a3 = -3.0*(d_pos*d_vel)*r2inv;
	    real rinv = sqrt(r2inv);
	    real mrinv = bj->get_mass() * rinv;
	    new_pot -= mrinv;
	    real mr3inv = mrinv * r2inv;
	    new_acc -= mr3inv*d_pos;
            new_jerk -= mr3inv*(d_vel + a3*d_pos);

	    if (!no_ssd_flag) {

		// Update nn and oscillation data.

		if (r2 < nn_dr2) {
		    nn_dr2 = r2;
		    nn_label = bj->get_index();
		    nn = bj;
		}
		if (r2 < min_nn_dr2) {
		    min_nn_dr2 = r2;
		    min_nn_label = bj->get_index();
		}

		((sdyn3 *) parent)->inc_ssd(0.5*r2);
	    }                                // correct for double counting

	    // Determine characteristic time scales:

	    real encounter_time_squared = 1 / (r2inv * d_vel * d_vel);
	    if (min_encounter_time_sq > encounter_time_squared)
		min_encounter_time_sq = encounter_time_squared;

	    real inverse_free_fall_time_squared = 
		r2inv * rinv * (mass + bj->get_mass());
	    if (min_free_fall_time_sq * inverse_free_fall_time_squared > 1)
		min_free_fall_time_sq = 1/inverse_free_fall_time_squared;
	}
}

void sdyn3::calculate_new_acc_and_jerk_from_new(sdyn3* b,
						real eps_squared,
						int  no_ssd_flag,
						int& collision_flag)
{
    if (oldest_daughter != NULL) {

	collision_flag = 0;

	if (!no_ssd_flag){
	    ssd = 0;
	    nn_label = 0;
	    for_all_daughters(sdyn3, this, c)
		c->set_nn_dr2(VERY_LARGE_NUMBER);
	}

        for_all_daughters(sdyn3, this, bb)
	    bb->calculate_new_acc_and_jerk_from_new(b, eps_squared,
						    no_ssd_flag,
						    collision_flag);
	if (!no_ssd_flag)
	    for_all_daughters(sdyn3, this, d) {
		if (d->get_init_nn_label() <= 0)
		    d->set_init_nn_label(d->get_nn_label());
		if (d->get_init_nn_label() != d->get_nn_label())
			d->set_nn_change_flag(1);
	    }

	if (!no_ssd_flag) {

	    if (min_min_ssd > ssd) min_min_ssd = ssd;
	    
	    if (ssd_ingoing_flag) {
		if (min_ssd > ssd) min_ssd = ssd;
		if (ssd > min_ssd * SSD_HYSTERESIS_FACTOR)
		    ssd_ingoing_flag = 0, max_ssd = ssd, n_ssd_osc++;
	    } else {
		if (max_ssd < ssd) max_ssd = ssd;
		if (ssd < max_ssd / SSD_HYSTERESIS_FACTOR)
		    ssd_ingoing_flag = 1, min_ssd = ssd;
	    }
	    
	    real e_tot = 0;
	    for_all_daughters(sdyn3, this, bbb)
		e_tot += 0.5 * bbb->get_mass() * (bbb->get_new_pot()
				    + bbb->get_new_vel() * bbb->get_new_vel());
	    if (de_tot_abs_max < abs(e_tot - e_tot_init))
	        de_tot_abs_max = abs(e_tot - e_tot_init);
	}

    } else {

        new_pot = 0;
	new_acc = new_jerk = 0;
	min_encounter_time_sq = min_free_fall_time_sq = VERY_LARGE_NUMBER;
	accumulate_new_acc_and_jerk_from_new(b, eps_squared, no_ssd_flag,
					     collision_flag);
    }
}

void  sdyn3::taylor_pred_new_pos_and_vel(const real dt)
{
    taylor_pred_new_pos(dt);
    taylor_pred_new_vel(dt);
}

void  sdyn3::correct_new_acc_and_jerk(const real new_dt, const real prev_new_dt)
{
    if (new_dt == prev_new_dt) return;

    real dt_off = new_dt - 0.5 * prev_new_dt;
                                    // offset from midpoint of prev_new_dt step
    real theta = 0.5 * prev_new_dt;
    real tau = dt_off / theta;          // equals 1 if new_dt = prev_new_dt

    real inv_theta = 1 / theta;
    real tau2 = tau * tau;
    real tau3 = tau2 * tau;
    vector prev_new_acc = new_acc;
    vector prev_new_jerk = new_jerk;

    new_acc = 0.25 * (
		      acc * (2 - 3 * tau + tau3)
		      + prev_new_acc * (2 + 3 * tau - tau3)
		      + jerk * theta * (1 - tau - tau2 + tau3)
		      + prev_new_jerk * theta * (-1 - tau + tau2 + tau3)
		      );

    new_jerk = 0.25 * (
		       acc * inv_theta * 3 * (-1 + tau2)
		       + prev_new_acc * inv_theta * 3 * (1 - tau2)
		       + jerk * (-1 - 2*tau + 3*tau2)
		       + prev_new_jerk * (-1 + 2*tau + 3*tau2)
		       );
}

void  sdyn3::correct_new_pos_and_vel(const real new_dt)  // timestep
{
    real new_dt2 = new_dt * new_dt;
    real new_dt3 = new_dt2 * new_dt;

    new_vel = vel + new_dt * (new_acc + acc)/2
	          - new_dt2 * (new_jerk - jerk)/12;

    new_pos = pos + new_dt * (new_vel + vel)/2 
	          - new_dt2 * (new_acc - acc)/12;          // alpha = -5/3
}

//----------------------------------------------------------------------

// Useful functions (used by low_n_evolve and scatter3):

void set_kepler_from_sdyn3(kepler& k, sdyn3* b1, sdyn3* b2)
{
    if (b1->get_time() != b2->get_time()) 
	err_exit("set_kepler_from_sdyn3: inconsistent times");

    k.set_time(b1->get_time());
    k.set_total_mass( b1->get_mass() + b2->get_mass() );
    k.set_rel_pos( b2->get_pos() - b1->get_pos() );
    k.set_rel_vel( b2->get_vel() - b1->get_vel() );

    k.initialize_from_pos_and_vel();
}

void kepler_pair_to_triple(kepler & k1,	// Inner binary (b1 + b2)
			   kepler & k2,	// Outer binary
			   sdyn3 * b1,
			   sdyn3 * b2,
			   sdyn3 * b3)
{
    // Construct positions and velocities from the specified kepler structures.

    // Assume that labels, masses, times, and the center of mass are handled
    // elsewhere.

    // Set up the outer orbit first.  Note that the center of mass of
    // the three-body system is taken to be at rest at the origin.

    // Make no assumptions about masses!

    real m1 = b1->get_mass();
    real m2 = b2->get_mass();
    real m3 = b3->get_mass();

    real m12 = m1 + m2;
    real m123 = m12 + m3;

    // Relative position and velocity from (1,2) to 3

    b3->set_pos(m12 * k2.get_rel_pos() / m123);
    b3->set_vel(m12 * k2.get_rel_vel() / m123);
    
    b1->set_pos(-m3 * k2.get_rel_pos() / m123);
    b1->set_vel(-m3 * k2.get_rel_vel() / m123);

    b2->set_pos(b1->get_pos());
    b2->set_vel(b1->get_vel());

    // Then set up the inner binary.

    b1->inc_pos(-m2 * k1.get_rel_pos() / m12);	    // rel_pos is from 1 to 2
    b1->inc_vel(-m2 * k1.get_rel_vel() / m12);

    b2->inc_pos( m1 * k1.get_rel_pos() / m12);
    b2->inc_vel( m1 * k1.get_rel_vel() / m12);

}
