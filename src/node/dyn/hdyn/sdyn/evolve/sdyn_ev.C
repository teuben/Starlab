
// sdyn_ev.C
// member functions for the sdyn class that are related to orbit integration.
//

#include "sdyn.h"

void sdyn::accumulate_new_acc_and_jerk_from_new(
			sdyn * bj,                  // n-body system pointer
			real eps2,                  // softening length squared
			int no_diag_flag,    // if true, intermediate iteration
			int & collision_flag)
    {
    if(bj->get_oldest_daughter() != NULL)
        for_all_daughters(sdyn, bj, bb)
	    accumulate_new_acc_and_jerk_from_new(bb, eps2,
						 no_diag_flag, collision_flag);
    else
	if(this != bj)
	    {
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

	    if (! no_diag_flag)
		{
		if (nn_dr2 > r2)
		    {
		    nn_dr2 = r2;
		    nn_label = bj->get_index();
		    }
		if (min_nn_dr2 > r2)
		    {
		    min_nn_dr2 = r2;
		    min_nn_label = bj->get_index();
		    }
		}

	    real encounter_time_squared = 1 / (r2inv * d_vel * d_vel);
	    if (min_encounter_time_sq > encounter_time_squared)
		min_encounter_time_sq = encounter_time_squared;

	    real inverse_free_fall_time_squared = 
		r2inv * rinv * (mass + bj->get_mass());
	    if (min_free_fall_time_sq * inverse_free_fall_time_squared > 1)
		min_free_fall_time_sq = 1/inverse_free_fall_time_squared;
	    }
    }

void sdyn::calculate_new_acc_and_jerk_from_new(sdyn * b, real eps_squared,
					       int  no_diag_flag,
					       int  & collision_flag)
    {
    if(oldest_daughter != NULL)
	{
	collision_flag = 0;

	if (! no_diag_flag)
	    {
	    nn_label = 0;
	    for_all_daughters(sdyn, this, c)
		c->set_nn_dr2(VERY_LARGE_NUMBER);
	    }

        for_all_daughters(sdyn, this, bb)
	    bb->calculate_new_acc_and_jerk_from_new(b, eps_squared,
						    no_diag_flag,
						    collision_flag);
	if (! no_diag_flag)
	    {
	    for_all_daughters(sdyn, this, d)
		{
		if (d->get_init_nn_label() <= 0)
		    d->set_init_nn_label(d->get_nn_label());
		if (d->get_init_nn_label() != d->get_nn_label())
			d->set_nn_change_flag(1);
		}

	    real e_tot = 0;
	    for_all_daughters(sdyn, this, bbb)
		e_tot += 0.5 * bbb->get_mass() * (bbb->get_new_pot()
				    + bbb->get_new_vel() * bbb->get_new_vel());
	    if (de_tot_abs_max < abs(e_tot - e_tot_init))
	        de_tot_abs_max = abs(e_tot - e_tot_init);
	    }
	}
    else
	{
        new_pot = 0;
	new_acc = new_jerk = 0;
	min_encounter_time_sq = min_free_fall_time_sq = VERY_LARGE_NUMBER;
	accumulate_new_acc_and_jerk_from_new(b, eps_squared, no_diag_flag,
					     collision_flag);
	}
    }

void  sdyn::taylor_pred_new_pos_and_vel(const real dt)
    {
    taylor_pred_new_pos(dt);
    taylor_pred_new_vel(dt);
    }

void  sdyn::correct_new_acc_and_jerk(const real new_dt, const real prev_new_dt)
    {
    if (new_dt == prev_new_dt)
	return;

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

void  sdyn::correct_new_pos_and_vel(const real new_dt)  // timestep
    {
    real new_dt2 = new_dt * new_dt;
    real new_dt3 = new_dt2 * new_dt;

    new_vel = vel + new_dt * (new_acc + acc)/2
	         - new_dt2 * (new_jerk - jerk)/12;

    new_pos = pos + new_dt * (new_vel + vel)/2 
	          - new_dt2 * (new_acc - acc)/12;          // alpha = -5/3
    }


