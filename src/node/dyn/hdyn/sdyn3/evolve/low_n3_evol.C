
//// low_n3_evol:  three-body time-symmetrized Hermite integrator.
////               (Called low_n3_evolve internally; name shortened
////               to satisfy some library managers.)
////
//// Options:      -A    specify accuracy parameter [0.02]
////               -c    specify CPU time check interval (s) [3600]
////               -C    specify cube size for snapshot output [10]
////               -d    specify log output interval [none]
////               -D    specify snapshot output interval [none]
////               -e    specify softening parameter [0.05]
////               -f    specify constant or dynamic timesteps [constant]
////               -n    specify number of iterations in symmetrization [0]
////               -q    quiet output [verbose]
////               -s    symmetric timestep [false]
////               -t    specify time span of integration [1]
////               -x    terminate at exact end time [no]
////               -z    specify maximum number of steps to take [unspecified]

//		  9/7/95:  Added "double unperturbed" motion near
//			   close binary periastron.  SMcM

// Starlab library function.

#include "sdyn3.h"

#ifndef TOOLBOX

#define UNPERTURBED_DIAG false
#define PRECISION 12

local void restore_pos_and_vel(kepler& inner, kepler& outer, real time,
			       int n_transform,
			       sdyn3* b1, sdyn3* b2, sdyn3* b3,
			       vector& cmpos, vector& cmvel)
{
    inner.transform_to_time(time);
    if (n_transform == 2) outer.transform_to_time(time);

    // Reconstruct the new sdyn3s.  Note that kepler_pair_to_triple
    // will set pos and vel, but we really want to change new_pos
    // and new_vel.

    kepler_pair_to_triple(inner, outer, b1, b2, b3);

    for_all_daughters(sdyn3, b1->get_parent(), bb) {
	bb->inc_pos(cmpos);
	bb->inc_vel(cmvel);
	bb->store_old_into_new();
    }
}

local real total_energy(sdyn3* b)
{
    real k = 0, u = 0, dissipation = 0;

    for_all_daughters(sdyn3, b, bi) {
	k += bi->get_mass() * bi->get_vel() * bi->get_vel();
	dissipation += bi->get_energy_dissipation();
	for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	    u -= bi->get_mass() * bj->get_mass()
		  / abs(bi->get_pos() - bj->get_pos());
    }

    return 0.5*k + u + dissipation;
}

static int save_flag = 1;
static real last_error = 0, last_test3 = 0;

local void print_unperturbed_diag(sdyn3* b, sdyn3* b1, sdyn3* b2, sdyn3* b3,
				  kepler& inner, kepler& outer,
				  char* header)
{
    if (UNPERTURBED_DIAG) {

	int p = cerr.precision(PRECISION);	// nonstandard precision

	if (header) cerr << header << endl;

	cerr << "  inner binary:"
	     << "    time = "
	     << inner.get_time() + (real)b->get_time_offset()
	     << "  t_peri = " << inner.get_time_of_periastron_passage()
		 			+ (real)b->get_time_offset() << endl
	     << "    r = " << inner.get_separation()
	     << "  a = " << inner.get_semi_major_axis()
	     << "  e = " << inner.get_eccentricity() << endl
	     << "    mass = " << inner.get_total_mass() << endl
	     << "    pos = " << inner.get_rel_pos() << endl
	     << "    vel = " << inner.get_rel_vel() << endl
	     << "    th  = " << inner.get_true_anomaly() << endl
	     << "    l = " << inner.get_longitudinal_unit_vector() << endl
	     << "    n = " << inner.get_normal_unit_vector() << endl;

	cerr << "  energy test 1:  " << inner.get_energy()
	     << " ?= " << 0.5 * square(b1->get_vel() - b2->get_vel())
		 	 - (b1->get_mass() + b2->get_mass())
			     / abs(b1->get_pos() - b2->get_pos()) << endl;

	cerr << "  outer binary:"
	     << "    R = " << outer.get_separation()
	     << "  a = " << outer.get_semi_major_axis()
	     << "  e = " << outer.get_eccentricity() << endl
	     << "    mass = " << outer.get_total_mass() << endl
	     << "    pos = " << outer.get_rel_pos() << endl
	     << "    vel = " << outer.get_rel_vel() << endl
	     << "    th  = " << outer.get_true_anomaly() << endl
	     << "    l = " << outer.get_longitudinal_unit_vector() << endl
	     << "    n = " << outer.get_normal_unit_vector() << endl;

	vector bcm_pos = (b1->get_mass() * b1->get_pos()
			   + b2->get_mass() * b2->get_pos())
	    		     / (b1->get_mass() + b2->get_mass());
	vector bcm_vel = (b1->get_mass() * b1->get_vel()
			   + b2->get_mass() * b2->get_vel())
	    		     / (b1->get_mass() + b2->get_mass());

	cerr << "  energy test 2:  " << outer.get_energy()
	     << " ?= " << 0.5 * square(b3->get_vel() - bcm_vel)
		 	 - (b1->get_mass() + b2->get_mass() + b3->get_mass())
			     / abs(b3->get_pos() - bcm_pos) << endl;

	real r = inner.get_separation();
	real R = outer.get_separation();
	cerr << "  check r: " << r << " ?= "
	     << abs(b2->get_pos() - b1->get_pos()) << endl;
	cerr << "  check R: " << R << " ?= "
	     << abs(b3->get_pos() - bcm_pos) << endl;

	// This term is basically the error incurred by the assumption
	// of double unperturbed motion (apart from higher-order terms):

	real m1 = b1->get_mass();
	real m2 = b2->get_mass();
	real m3 = b3->get_mass();
	real test3 = - m1*m3 / abs(b3->get_pos() - b1->get_pos())
		     - m2*m3 / abs(b3->get_pos() - b2->get_pos())
		     + (m1+m2)*m3 / outer.get_separation();

	real coeff = (m1*m2/(m1+m2)) * (0.5*m3/pow(R,3));
	real term1 = coeff * r*r;
	real term2 = coeff * (-3) * pow(inner.get_rel_pos()
				     * outer.get_rel_pos() / R, 2);

	cerr << "  energy test 3:  " << test3
	     << " =? " << term1 + term2 << endl;

	real error = total_energy(b) - b->get_e_tot_init();
	cerr << "  energy error = " << error << endl;

	if (save_flag) {
	    last_test3 = test3;
	    last_error = error;
	} else {
	    cerr << "  delta(test3) = "
		 << test3 - last_test3 << endl;
	    cerr << "  delta(error) = "
		 << error - last_error << endl;
	}

	save_flag = 1 - save_flag;
	cerr.precision(p);

	cerr << flush;
    }
}

local void print_perturbed_error(kepler& inner, kepler& outer,
				 sdyn3* b1, real dt)
{
    // On entry, inner is in its initial state; outer has been
    // advanced to the final time.

    cerr << "  inner, outer time = " << inner.get_time()
	 << "  " << outer.get_time() << endl;

    // Evolve outer binary to the time of inner periastron and get
    // relevant parameters.

    outer.transform_to_time(outer.get_time() - 0.5*dt);

    real M = outer.get_total_mass() - inner.get_total_mass();
    real R = abs(outer.get_separation());
    real X = outer.get_rel_pos() * inner.get_longitudinal_unit_vector();
    real Y = outer.get_rel_pos() * inner.get_transverse_unit_vector();
    real Z = outer.get_rel_pos() * inner.get_normal_unit_vector();

    // Evaluate the inner binary integral.

    real E = inner.get_energy();
    real e = inner.get_eccentricity();

    if (E >= 0 || e <= 0 || e >= 1) return;	// Safety check!

    real a = inner.get_semi_major_axis();
    real h = inner.get_angular_momentum();
    real r = inner.get_separation();
    real th = inner.get_true_anomaly();

    real rp = a * (1 - e);
    real ra = a * (1 + e);

    // Attempt to predict the energy change in the quadrupole approximation.

    real red_mass = b1->get_mass()
			* (1 - b1->get_mass()/inner.get_total_mass());
    real coeff = 6*red_mass*M*X*Y/pow(R,5);

    // Energy terms come from algebra and Maple, and are highly suspect!

    real mu = cos(th);
    real I1 = a*a*pow(1-e*e,2) *
		(-sin(th) * (2 - e*e + mu*e*(3-2*e*e))
		          / (e*(1-e*e)*pow(1+mu*e,2))
		 + (asin((mu+e)/(1+mu*e))-M_PI/2)
		 	  * (2-3*e*e) / (e*e*pow(1-e*e,1.5))
		 + 2*th / (e*e));

    real f1 = pow(a*e,2) - pow(a-r,2);
    real f2 = atan((r-a*(1-e*e))/sqrt((1-e*e)*f1)) - M_PI/2;
    real f3 = asin((a-r)/(a*e)) - M_PI/2;
    real I2 = (h/(e*e*sqrt(2*abs(E)))) *
		((e*e-2)*sqrt(f1) + 2*a*pow(1-e*e,1.5)*f2 - a*f3*(2-3*e*e));

    cerr << "  predicted potential change = "
	 << coeff * r*r * cos(th) * sin(th)
	 << endl;
    cerr << "  predicted error = " << coeff * (I1 + I2) << "\n";
    cerr << "  terms:  I1, I2 = " << I1 << "  " << I2 << endl;

    // The terms for dv are more reliable...(?)

    real j0 = asin((a-r)/(a*e)) - M_PI/2;
    real term1 = 1.5*a*a*e*e * j0;
    real term2 = (a*(0.5+e*e)+0.5*r) * sqrt(a*a*e*e-(a-r)*(a-r));
    real I = (term1 + term2) / (e*sqrt(2*abs(E)));
    real coeffdv = -2*M/pow(R,5) * I;

    vector dv = coeffdv * (-(R*R - 3*X*X) * inner.get_longitudinal_unit_vector()
			         - 3*X*Y  * inner.get_transverse_unit_vector()
			         - 3*X*Z  * inner.get_normal_unit_vector());

    // We have added a "-" to the longitudinal term here because inner
    // is not yet advanced to the final time (v_long --> -v_long).

    cerr << "  predicted dv = " << dv << endl;
    cerr << "  associated dE = " << red_mass*dv*inner.get_rel_vel() << endl;
    cerr << "  terms:  " << term1 << "  " << term2 << endl;
}

static real tol_23 = -1;

local bool unperturbed_step(sdyn3* b,		// n-body system pointer
			    real tidal_tolerance,
			    real& true_dt,	// new timestep
			    int&  collision_flag)
{
    // Check for and implement unperturbed motion.

    // Only take an unperturbed step if the unperturbed criterion is met
    // by both the old and the final systems.  (This is necessary for time
    // symmetry, and is probably a good idea anyway.)

    // Return TRUE iff an unperturbed step was taken.
    // If so,
    //    1. replace new_pos and new_vel by updated quantities,
    //    2. set true_dt and collision_flag,

    // Need to make the initial cut reasonably efficient.  Use the
    // nearest-neighbor information returned by the integrator: each
    // body has nn, nn_dr2, nn_label set.  We may even want to copy
    // SJA and have some sort of time-step criterion!

    if (tidal_tolerance <= 0) return FALSE;

    sdyn3* bi = b->get_oldest_daughter();
    sdyn3* bj = bi->get_younger_sister();
    sdyn3* bk = bj->get_younger_sister();

    real min_sep_sq = min(min(bi->get_nn_dr2(), bj->get_nn_dr2()),
			  bk->get_nn_dr2());
    real max_sep_sq = max(max(bi->get_nn_dr2(), bj->get_nn_dr2()),
			  bk->get_nn_dr2());

    // The following mass factor represents a "worst-case" scenario
    // (see the comments about tidal factors in scatter3.C).  For
    // binary components 1 and 2, the tidal term should really be
    //
    //		(m1 + m2) * tolerance / (m1 + m2 + m3).
    //
    // For now, at least, we take the minimum possible binary mass.

    if (tol_23 < 0) {
	real total_mass = bi->get_mass() + bj->get_mass() + bk->get_mass();
	real max_mass = max(max(bi->get_mass(), bj->get_mass()),
			    bk->get_mass());
	tol_23 = pow((1 - max_mass/total_mass) * tidal_tolerance, 0.6666667);
    }

    if (min_sep_sq > max_sep_sq * tol_23) return FALSE;

    // -------------------------------------------------------------------
    // Passed the first cut!

    // Now determine the binary components, and calculate the separation
    // between the binary center of mass and the third body more carefully.

    sdyn3 *b1, *b2, *b3;

    // Find a component of the binary:

    if (min_sep_sq == bi->get_nn_dr2())
	b1 = bi;
    else if (min_sep_sq == bj->get_nn_dr2())
	b1 = bj;
    else
	b1 = bk;

    // Other component:

    b2 = (sdyn3*)b1->get_nn();
    if (!b2) return FALSE;	// No neighbor defined yet.

    // Third star:

    if (bi != b1 && bi != b2)
	b3 = bi;
    else if (bj != b1 && bj != b2)
	b3 = bj;
    else
	b3 = bk;

    // Binary components must be approaching.

    vector r = b2->get_pos() - b1->get_pos();

    if (r * (b2->get_vel() - b1->get_vel()) >= 0) return FALSE;

    // Binary components must be close compared to the distance
    // to the third body

    real m12 = b1->get_mass() + b2->get_mass();
    real m123 = m12 + b3->get_mass();

    // Center of mass position of the inner binary:
    
    vector binary_cmpos = (b1->get_mass() * b1->get_pos()
			    + b2->get_mass() * b2->get_pos()) / m12;

    vector R = b3->get_pos() - binary_cmpos;

    real true_tol_23 = 	pow((m12/m123) * tidal_tolerance, 0.6666667);
    if (square(r) > square(R) * true_tol_23) return FALSE;

    // -------------------------------------------------------------------
    // We have unperturbed motion.  Keep track of the system center of
    // mass (which will be lost by the kepler structures below).

    vector system_cmpos = (m12 * binary_cmpos
			   + b3->get_mass() * b3->get_pos()) / m123;

    // Center of mass velocity of the inner binary:

    vector binary_cmvel = (b1->get_mass() * b1->get_vel()
			    + b2->get_mass() * b2->get_vel()) / m12;

    vector system_cmvel = (m12 * binary_cmvel
			    + b3->get_mass() * b3->get_vel()) / m123;

    // Make kepler structures out of the inner and outer orbits.

    kepler inner, outer;

    // Inner orbit:

    set_kepler_from_sdyn3(inner, b1, b2);

    // Outer orbit:

    outer.set_time(b3->get_time());
    outer.set_total_mass(m123);    
    outer.set_rel_pos(R);
    outer.set_rel_vel(b3->get_vel() - binary_cmvel);
 
    outer.initialize_from_pos_and_vel();

    real t_init = inner.get_time();
    real t_peri = inner.get_time_of_periastron_passage();

    // Update the closest-approach variables.  Assume that the pointers
    // and labels don't need to be changed.

    real peri_sq = inner.get_periastron() * inner.get_periastron();
    b1->set_nn_dr2(min(b1->get_nn_dr2(), peri_sq));
    b2->set_nn_dr2(min(b2->get_nn_dr2(), peri_sq));
    b1->set_min_nn_dr2(min(b1->get_min_nn_dr2(), peri_sq));
    b2->set_min_nn_dr2(min(b2->get_min_nn_dr2(), peri_sq));

    print_unperturbed_diag(b, b1, b2, b3,
			   inner, outer,
			   "Beginning unperturbed motion");

    // Check for a collision at binary periastron.

    if (inner.get_periastron() <= b1->get_radius() + b2->get_radius()) {

	// Move to periastron and return.

	restore_pos_and_vel(inner, outer, t_peri, 2,
			    b1, b2, b3,
			    system_cmpos, system_cmvel);
	collision_flag = 1;
	true_dt = t_peri - t_init;

	return TRUE;
    }

    real dt = 2*(t_peri - t_init);

    // Require that both the old and the final configurations pass
    // the "unperturbed" test.

    outer.transform_to_time(t_init + dt);

    real sep = outer.get_separation();
    if (square(r) > sep * sep * true_tol_23) {

	// Restore the original system.

	if (UNPERTURBED_DIAG)
	    cerr <<
    "unperturbed_step: failed unperturbed test after unperturbed step!\n";

	restore_pos_and_vel(inner, outer, t_init, 2, b1, b2, b3,
			    system_cmpos, system_cmvel);

	save_flag = 1;
	return FALSE;
    }

    if (UNPERTURBED_DIAG) print_perturbed_error(inner, outer, b1, dt);

    // Advance the inner orbit to an outgoing inner binary at the
    // same separation (outer orbit has already been transformed).

    restore_pos_and_vel(inner, outer, t_init + dt, 1, b1, b2, b3,
			system_cmpos, system_cmvel);
    true_dt = dt;

    print_unperturbed_diag(b, b1, b2, b3,
			   inner, outer,
			   "Ending unperturbed motion");

    // Note that old and new quantities are the same on exit.

    return TRUE;
}

local void predict(sdyn3* b,		// n-body system pointer
		   real dt)		// timestep
{
    // Predict the entire system from the current time to time + dt.

    if(b->get_oldest_daughter() !=NULL)
        for_all_daughters(sdyn3, b, bb)
	    predict(bb, dt);
    else
        b->taylor_pred_new_pos_and_vel(dt);
}

local void correct(sdyn3* b,		// n-body system pointer
		   real new_dt,		// new timestep
		   real prev_new_dt)	// previous new timestep
{
    // Apply Hermite corrector step to the entire system.

    if(b->get_oldest_daughter() !=NULL)
        for_all_daughters(sdyn3, b, bb)
	    correct(bb, new_dt, prev_new_dt);
    else {
	b->correct_new_acc_and_jerk(new_dt, prev_new_dt);
        b->correct_new_pos_and_vel(new_dt);
    }
}

local void calculate_energy(sdyn3* b, real& ekin, real& epot)
{
    ekin = epot = 0;
    for_all_daughters(sdyn3, b, bb) {
	epot += bb->get_mass() * bb->get_pot();
	ekin += 0.5 * bb->get_mass() * (bb->get_vel() * bb->get_vel());
    }
    epot *= 0.5;
}

typedef real (*tfp)(sdyn3*, real);

local bool step(sdyn3* b,	// sdyn3 array
		real& t,	// time
		real eps,	// softening length
		real eta,	// time step parameter
		real dt,	// time step of the integration 
		real max_dt,	// maximum time step (to end of the integration)
		int& end_flag,  // to flag that integration end is reached
		tfp  the_tfp,	// timestep function pointer
		int  n_iter,	// number of iterations
		int  x_flag,	// exact-time termination flag
		int  s_flag)	// symmetric timestep ?
{
    // Take a timestep, symmetrizing (n_iter iterations) if requested.
    // Note that the actual time step taken will in general not equal dt.

    bool unpert_flag = FALSE;

    real true_dt = dt; 
    int collision_flag = 0;

    if (eps > 0 ||
	!(unpert_flag = unperturbed_step(b, DEFAULT_TIDAL_TOL_FACTOR,
					 true_dt, collision_flag))) {

	// Predict all particles to time + dt.

	predict(b, dt);

	real new_dt = dt; 
	for (int i = 0; i <= n_iter; i++) {

	    // Get acc and jerk from current "predicted" quantities.

	    b->calculate_new_acc_and_jerk_from_new(b, eps*eps,
						   n_iter - i,	   // hack
						   collision_flag);

	    real prev_new_dt = new_dt;
	    if (s_flag && i < n_iter) {

		// Iterate to enforce time symmetry.

		real end_point_dt = the_tfp(b, eta);

		// Piet's mysterious accelerator:

		new_dt = 0.5 * (dt + end_point_dt);
		new_dt = dt + 0.5 * (end_point_dt - dt) * (new_dt/prev_new_dt);

		// See if we must force termination at a particular time.

		if (x_flag) {
		    if (new_dt >= max_dt) {
			end_flag = 1;
			new_dt = max_dt;
		    } else
			end_flag = 0;
		}
	    }

	    // Correct all particles for this (sub-)timestep.

	    correct(b, new_dt, prev_new_dt);
	}

	true_dt = new_dt;

    } else	// Need acc and jerk if we aren't about to quit.

	if (!collision_flag)
	    b->calculate_new_acc_and_jerk_from_new(b, eps*eps,
						   0, collision_flag);
	
    // Update all dynamical quantities:

    for_all_daughters(sdyn3, b, bb) bb->store_new_into_old();
    t += true_dt;

    if (collision_flag)	end_flag = 1;

    return unpert_flag;
}

local void initialize(sdyn3* b,	// sdyn3 array                   
		      real eps)	// softening length             
{
    predict(b, 0);

    int collision_flag;
    b->calculate_new_acc_and_jerk_from_new(b, eps*eps, 1, collision_flag);

    for_all_daughters(sdyn3, b, bb)
	bb->store_new_into_old();
}

// NOTE: If we want to use the_tfp as a variable function pointer,
// we MUST NOT make these timestep functions local!

real constant_timestep(sdyn3* b, real eta)
{
    if (b == NULL) eta = 0; // To keep an HP compiler happy...
    return  eta;
}

real dynamic_timestep(sdyn3* b, real eta)
{
    real global_min_encounter_time_sq = VERY_LARGE_NUMBER;
    real global_min_free_fall_time_sq = VERY_LARGE_NUMBER;

    for_all_daughters(sdyn3, b, bb) {
	global_min_encounter_time_sq =
	    min(global_min_encounter_time_sq, bb->get_min_encounter_time_sq());
	global_min_free_fall_time_sq =
	    min(global_min_free_fall_time_sq, bb->get_min_free_fall_time_sq());
    }

    return eta *
	sqrt(min(global_min_encounter_time_sq, global_min_free_fall_time_sq));
}

local tfp get_timestep_function_ptr(char* timestep_name)
{
    if (streq(timestep_name, "constant_timestep"))
	return constant_timestep;
    else if (streq(timestep_name, "dynamic_timestep"))
	return dynamic_timestep;
    else {
	cerr << "get_timestep_function_ptr: no timestep function implemented"
	     << " with name `" << timestep_name << "'" << endl;
	exit(1);
    }
    return (tfp)NULL; // To keep an HP g++ happy!
}

local void start_up(sdyn3* b, real& n_steps)
{
    if (b->get_oldest_daughter() !=NULL) {

	if ((n_steps = b->get_n_steps()) == 0)
	    b->prepare_root();

	for_all_daughters(sdyn3, b, bb)
	    start_up(bb, n_steps);

    } else
	if (n_steps == 0) b->prepare_branch();
}

local void clean_up(sdyn3 * b, real n)
{
    b->set_n_steps(n);
}

int system_in_cube(sdyn3* b, real cube_size)
{
    for_all_daughters(sdyn3, b, bb)
	for (int k = 0; k < 3; k++)
	    if (abs(bb->get_pos()[k]) > cube_size) return 0;
    return 1;
}

#define N_STEP_CHECK 1000	   // Interval for checking CPU time

#define PERT_OUTPUT 2
real last_unpert = -VERY_LARGE_NUMBER;

void low_n3_evolve(sdyn3* b,	   // sdyn3 array
		   real delta_t,   // time span of the integration
		   real dt_out,	   // output time interval
		   real dt_snap,   // snapshot output interval
		   real snap_cube_size,
		   real eps,	   // softening length 
		   real eta,	   // time step parameter
		   int  x_flag,	   // exact-time termination flag
		   char* timestep_name,
		   int  s_flag,	   // symmetric timestep ?
		   int  n_iter,	   // number of iterations
		   real n_max,	   // if > 0: max. number of integration steps
		   real cpu_time_check,
		   real dt_print,  // external print interval
		   sdyn3_print_fp
		        print)	   // pointer to external print function
{
    real t = b->get_time();

    real t_end = t + delta_t;	   // final time, at the end of the integration
    real t_out = t + dt_out;	   // time of next diagnostic output
    real t_snap = t + dt_snap;	   // time of next snapshot;
    real t_print = t + dt_print;   // time of next printout;

    bool unpert;

    real n_steps;
    int  count_steps = 0;
    real cpu_init = cpu_time();
    real cpu_save = cpu_init;

    tfp the_tfp = get_timestep_function_ptr(timestep_name);

    start_up(b, n_steps);
    initialize(b, eps);

    int p = cerr.precision(PRECISION);

    real ekin, epot;
    calculate_energy(b, ekin, epot);

    if (t_out <= t_end && n_steps == 0) {
	cerr << "Time = " << t + (real)b->get_time_offset()
	     << "  n_steps = " << n_steps
	     << "  Etot = " << ekin + epot << endl;
    }

    if (b->get_n_steps() == 0) {	   // should be better interfaced
					   // with start_up: to be done
	b->set_e_tot_init(ekin + epot);
	b->clear_de_tot_abs_max();
    }

    int end_flag = 0;
    while (t < t_end && !end_flag) {

        real max_dt = t_end - t;
        real dt = the_tfp(b, eta);

        end_flag = 0;
        if (dt > max_dt && x_flag) {
	    end_flag = 1;
	    dt = max_dt;
	}

        unpert = step(b, t, eps, eta, dt, max_dt, end_flag, the_tfp, n_iter,
		      x_flag, s_flag);

	// Time may be offst.  System time will be correct.

	b->set_system_time(b->get_system_time()+dt);

	b->set_time(t);			   // should be prettified some time
	for_all_daughters(sdyn3, b, bi)
	    bi->set_time(t);

	n_steps += 1;			   // (Note that n_steps is real)
	count_steps++;

	if (unpert) last_unpert = n_steps;

//      Check for (trivial) output to cerr...

	if (t >= t_out
	    || (UNPERTURBED_DIAG && n_steps - last_unpert < PERT_OUTPUT)) {
	    calculate_energy(b, ekin, epot);
	    int pp = cerr.precision(PRECISION);
	    cerr << "Time = " << t + (real)b->get_time_offset()
		 << "  n_steps = " << n_steps
		 << "  dE = " << ekin + epot - b->get_e_tot_init() << endl;
	    while (t >= t_out) t_out += dt_out;
	    cerr.precision(pp);
	}

//      ...and (not-so-trivial) output handled elsewhere.

	if (print && t >= t_print) {
	    (*print)(b);
	    t_print += dt_print;
	}

//      Output a snapshot to cout at the scheduled time, or at end of run.

        if (n_max > 0 && n_steps >= n_max) end_flag = 1;

	if (end_flag || t >= t_snap) {
	    if ((t >= t_snap) && system_in_cube(b, snap_cube_size)) {
		put_node(cout, *b);
		cout << flush; 
		t_snap += dt_snap;         // too early to get clean_up info?
	    }
	}

//      Check the number of steps and the CPU time every N_STEP_CHECK steps.
//      Note that the printed CPU time is the time since this routine was
//      entered.

	if (count_steps >= N_STEP_CHECK) {

	    count_steps = 0;

	    if (b->get_n_steps() > MAX_N_STEPS) return;

	    if (cpu_time() - cpu_save > cpu_time_check) {
		cpu_save = cpu_time();
		calculate_energy(b, ekin, epot);
		int p = cerr.precision(STD_PRECISION);
		cerr << "low_n3_evolve:  CPU time = " << cpu_save - cpu_init;
		cerr.precision(PRECISION);
		cerr << "  time = " << t + (real)b->get_time_offset();
		cerr.precision(STD_PRECISION);
		cerr << "  offset = " << b->get_time_offset() << endl;
		cerr << "                n_steps = " << n_steps
		     << "  Etot = " << ekin + epot
		     << "  dt = " << dt
		     << endl << flush;
		cerr.precision(p);
	    }
	}
    }

    cerr.precision(p);
    clean_up(b, n_steps);       // too late for snapshot?
}

#else

main(int argc, char **argv)
{
    sdyn3* b;		   // pointer to the nbody system
    int  n_iter = 0;	   // number of iterations (0: explicit; >=1: implicit)
    real n_max = -1;	   // if > 0: maximum number of integration steps
    
    real delta_t = 1;	   // time span of the integration
    real eta = 0.02;	   // time step parameter (for fixed time step,
			   //   equal to the time step size; for variable
			   //   time step, a multiplication factor)

    real dt_out = VERY_LARGE_NUMBER;
			   // output time interval
    real dt_snap = VERY_LARGE_NUMBER;
			   // snap output interval
    real snap_cube_size = 10;

    real cpu_time_check = 3600;

    real eps = 0.05;	   // softening length 	       	   
    char* timestep_name = "constant_timestep";

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:e:f:n:qr:st:xz:";

    bool a_flag = FALSE;
    bool d_flag = FALSE;
    bool D_flag = FALSE;
    bool e_flag = FALSE;
    bool f_flag = FALSE;
    bool n_flag = FALSE;
    bool q_flag = FALSE;
    bool r_flag = FALSE;
    bool s_flag = FALSE;   // symmetric timestep ?
    bool t_flag = FALSE;
    bool x_flag = FALSE;   // if true: termination at the exact time of
                           //          of the final output, by
                           //          adjustment of the last time step;
                           // if false: no adjustment of the last time step,
                           //           as a consequence the time of final
                           //           output might be slightly later than
                           //           the time specified.
    bool  z_flag = FALSE;  // to specify termination after n_max steps

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'A': a_flag = TRUE;
		      eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = atof(poptarg);
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': d_flag = TRUE;
		      dt_out = atof(poptarg);
		      break;
	    case 'D': D_flag = TRUE;
		      dt_snap = atof(poptarg);
		      break;
	    case 'e': e_flag = TRUE;
		      eps = atof(poptarg);
		      break;
	    case 'f': f_flag = TRUE;
		      timestep_name = poptarg;
		      break;
	    case 'n': n_flag = TRUE;
		      n_iter = atoi(poptarg);
		      break;
	    case 'q': q_flag = TRUE;
		      break;
	    case 's': s_flag = TRUE;
		      break;
	    case 't': t_flag = TRUE;
		      delta_t = atof(poptarg);
		      break;
	    case 'x': x_flag = TRUE;
		      break;
	    case 'z': z_flag = TRUE;
		      n_max = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
	}            

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag) cerr << "default delta_t = " << delta_t << "\n";
	if (!a_flag) cerr << "default eta = " << eta << "\n";
	if (!d_flag) cerr << "default dt_out = " << dt_out << "\n";
	if (!e_flag) cerr << "default eps = " << eps << "\n";
	if (!f_flag) cerr << "default timestep_name = " << timestep_name
	                  << "\n";
	if (!n_flag) cerr << "default n_iter = " << n_iter << "\n";
	if (!s_flag) cerr << "s_flag (symmetric timestep ?) = FALSE" << "\n";
	if (!x_flag) cerr << "default termination: not at exact t_end" << "\n";
	if (!z_flag) cerr << "default n_max = " << n_max << "\n";
    }

    if (!D_flag) dt_snap = delta_t;

    b = get_sdyn3(cin);
    b->log_history(argc, argv);
    cpu_init();

    low_n3_evolve(b, delta_t, dt_out, dt_snap, eps, eta, snap_cube_size,
		  x_flag, timestep_name, s_flag, n_iter, n_max,
		  cpu_time_check);
}

#endif
