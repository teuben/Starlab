
// scatter.C: Perform three-body scattering experiments.
//
// Input is an initial state, output is intermediate and final states.
// No reference is made to scatter profiles in this file.

// Note the function scatter3 is completely deterministic.
// No randomization is performed at this level.

// The program scatter3 randomizes all phase angles prior to invoking
// the scatter3 function.  Other parameters are fully determined, either
// by default or on the command line.

#include "scatter3.h"

// Set up the orientation of the outer orbit with respect to the inner
// binary [which always lies in the (x-y) plane].

// See "sigma3.h" for details.

local void set_orientation(kepler &k, phase3 &p)
{
    real mu = p.cos_theta;
    real sin_theta = sqrt(1 - mu * mu);

    // Construct the normal vector:

    vec n = vec(sin_theta*cos(p.phi), sin_theta*sin(p.phi), mu);

    // Construct unit vectors a and b perpendicular to n:

    vec temp = vec(1, 0, 0);
    if (abs(n[0]) > .5) temp = vec(0, 1, 0);	// temp is not parallel to n
    if (n[2] < 0) temp = -temp;

    vec b = n ^ temp;
    b /= abs(b);
    vec a = b ^ n;
    if (n[2] < 0) a = -a;	// Force (a, b) to be (x, y) for n = +/-z
    
    // Construct *random* unit vectors l and t perpendicular to each
    // other and to n (psi = 0 ==> periastron along a):

    vec l = cos(p.psi)*a + sin(p.psi)*b;
    vec t = n ^ l;

    k.set_orientation(l, t, n);
    k.initialize_from_shape_and_phase();
}

local real potential_energy(sdyn3 * b)
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

local real energy(sdyn3 * b)
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

local void kepler_pair_to_triple(kepler & k1,	// Inner binary (b1 + b2)
			         kepler & k3,	// Outer binary
			         sdyn3 * b1, sdyn3 * b2, sdyn3 * b3)

// Construct positions and velocities from the specified kepler structures.
// Assume that labels. masses and times are handled elsewhere.

{
    // Set up the outer orbit first.  Note that the center of mass of
    // the three-body system is taken to be at rest at the origin.

    // Make no assumptions about masses!

    real m1 = b1->get_mass();
    real m2 = b2->get_mass();
    real m3 = b3->get_mass();

    real m12 = m1 + m2;
    real m123 = m12 + m3;

    // relative position and velocity from (1,2) to 3

    b3->set_pos(m12 * k3.get_rel_pos() / m123);
    b3->set_vel(m12 * k3.get_rel_vel() / m123);
    
    b1->set_pos(-m3 * k3.get_rel_pos() / m123);
    b1->set_vel(-m3 * k3.get_rel_vel() / m123);

    b2->set_pos(b1->get_pos());
    b2->set_vel(b1->get_vel());

    // Then set up the inner binary.

    b1->inc_pos(-m2 * k1.get_rel_pos() / m12);	    // rel_pos is from 1 to 2
    b1->inc_vel(-m2 * k1.get_rel_vel() / m12);

    b2->inc_pos( m1 * k1.get_rel_pos() / m12);
    b2->inc_vel( m1 * k1.get_rel_vel() / m12);

}

local sdyn3 * set_up_dynamics(real m2,      // mass of secondary in binary
			      real m3,      // mass of projectile (m1 + m2 = 1)
			      kepler & k1,  // inner binary
			      kepler & k3)  // outer binary
    {
    // Create bodies for integration.

    sdyn3 * b = mksdyn3(3);            	   // pointer to the 3-body system

    // Initialize the dynamics.

    b->set_time(k3.get_time());
    for_all_daughters(sdyn3, b, bb)
	bb->set_time(k3.get_time());

    sdyn3 * b1 = b->get_oldest_daughter();
    sdyn3 * b2 = b1->get_younger_sister();
    sdyn3 * b3 = b2->get_younger_sister();

    b1->set_label(1);
    b2->set_label(2);
    b3->set_label(3);

    b->set_mass(1 + m3);
    b1->set_mass(1 - m2);
    b2->set_mass(m2);
    b3->set_mass(m3);

    kepler_pair_to_triple(k1, k3, b1, b2, b3);

    return b;
    }

// init_to_sdyn3: convert an initial state to an sdyn3 for integration.

local sdyn3 * init_to_sdyn3(initial_state3 & init, final_state3 & final)
{
    set_kepler_tolerance(2);

    kepler k1;	// Inner binary.

    real peri = 1; // Default value (unimportant unless init.ecc = 1).
    if (init.ecc == 1) peri = 0;

    make_standard_kepler(k1, 0, 1, -0.5 / init.sma, init.ecc,
			 peri, init.phase.mean_anomaly);

//    cout << endl << "Inner binary:" << endl;
//    k1.print_all(cout);

    // Multiply the incoming velocity by vcrit.

    real m2 = init.m2;
    real m3 = init.m3;
    real mtotal = 1 + m3;
    real v_inf = init.v_inf * sqrt( (1 - m2) * m2 * mtotal / m3 );
    
    real energy3 = .5 * v_inf * v_inf;
    real ang_mom3 = init.rho * v_inf;

    real ecc3 = (energy3 == 0 ? 1
		              : sqrt( 1 + 2 * energy3 
				            * pow(ang_mom3/mtotal, 2)));

    // Special case: if v_inf = 0, assume rho is periastron.

    real virial_ratio = init.rho;

    kepler k3;	// Outer binary.

    make_standard_kepler(k3, 0, mtotal, energy3, ecc3, virial_ratio, 0);

    // Radius for "unperturbed" inner binary:

    real r_unp = (init.sma + init.r3)
	             * pow(TIDAL_TOL_FACTOR / mtotal, -1/3.0);

    if (r_unp <= k3.get_periastron()) {
	final.descriptor = preservation;
	final.sma = k1.get_semi_major_axis();
	final.ecc = k1.get_eccentricity();
	final.outer_separation = k3.get_periastron();
	final.escaper = 3;
	final.error = 0;
	final.time = 0;
	final.n_steps = 0;
	final.virial_ratio = k3.get_energy() * k3.get_periastron() 
	                                        / k3.get_total_mass() - 1;

	// Sufficiently large stellar radii (sum larger than the binary
	// periastron) will cause an immediate collision to occur, so we
	// must check for it explicitly here.

	if (k1.get_periastron() < init.r1 + init.r2) {
	    final.descriptor = merger_escape_3;
	    final.sma = final.ecc = -1;
	}

	return NULL;
    }

    init.r_init = max(init.r_init_min, min(init.r_init_max, r_unp));

    k1.transform_to_time(k3.return_to_radius(init.r_init));

    set_orientation(k3, init.phase);

//    cout << endl << "Outer binary:" << endl;
//    k3.print_all(cout);

    sdyn3 * b;
    b = set_up_dynamics(m2, m3, k1, k3);

//  Set up radii:

    sdyn3 * bb = b->get_oldest_daughter();
    bb->set_radius(init.r1);
    bb = bb->get_younger_sister();
    bb->set_radius(init.r2);
    bb = bb->get_younger_sister();
    bb->set_radius(init.r3);

    return  b;
}

local void extend_orbits(sdyn3 * b1, sdyn3 * b2, sdyn3 * b3)

// Analytically extend the orbits of the [[1,2],3] hierarchical system.

    {
    if ( (b1->get_time() != b2->get_time())
	 || (b1->get_time() != b3->get_time()) )
	err_exit("extend_orbits: inconsistent times");
    
    real time = b3->get_time();
    
    kepler inner, outer;
    
    set_kepler_from_sdyn3(inner, b1, b2);
    
    real m1 = b1->get_mass();
    real m2 = b2->get_mass();
    real m3 = b3->get_mass();
    real m12 = m1 + m2;
    real m123 = m12 + m3;
    
    outer.set_time(time);
    outer.set_total_mass(m123);
    
    // Center of mass position and velocity of the inner binary:
    
    vec cmr = (m1 * b1->get_pos() + m2 * b2->get_pos()) / m12;
    vec cmv = (m1 * b1->get_vel() + m2 * b2->get_vel()) / m12;
    
    outer.set_rel_pos(b3->get_pos() - cmr);
    outer.set_rel_vel(b3->get_vel() - cmv);
    
    outer.initialize_from_pos_and_vel();
    
    // Map forward in time by an integral number of inner binary periods
    // to an incoming configuration as close as possible to the outgoing one.
    
    int inner_orbits = 2 * (int) ((outer.pred_advance_to_apastron() - time)
				  / inner.get_period() + 0.5);
    time += inner_orbits * inner.get_period();
    outer.transform_to_time(time);
    
    // Note that no modification is needed to the inner orbit
    // except to increase its time.
    
    inner.set_time(time);
    
    // Reconstruct the new sdyn3s:
    
    b1->get_parent()->set_time(time);
    b1->set_time(time);
    b2->set_time(time);
    b3->set_time(time);
    
    kepler_pair_to_triple(inner, outer, b1, b2, b3);
    
//    cerr << "extend_orbits:  " << inner_orbits << " inner orbits" << endl;
//    cerr << "outer orbit: " << endl;
//    outer.print_elements(cerr);
//    cerr << "inner orbit: " << endl;
//    inner.print_elements(cerr);
    
    }

local int escape(sdyn3 * bi, sdyn3 * bj, sdyn3 * bk,
		 real ejk, 				//(for convenience)
		 real r_stop, 
		 intermediate_state3 & inter,
		 final_state3 & final)

// Return true iff bi is escaping from the center of mass of bj and bk.
// Perform analytic continuation if appropriate.
// On entry, we already know that i is "widely" separated from j and k.

{
    final.descriptor = unknown_final;
    if (ejk >= 0) return 0;		// Inner pair not bound.

    real mi = bi->get_mass();
    real mj = bj->get_mass();
    real mk = bk->get_mass();
    real mjk = mj + mk;

    // Center of mass position and velocity of j-k "binary":

    vec cmr = (mj * bj->get_pos() + mk * bk->get_pos()) / mjk;
    vec cmv = (mj * bj->get_vel() + mk * bk->get_vel()) / mjk;

    // Return immediately if third particle is approaching binary.

    vec dv = bi->get_vel() - cmv;
    vec dr = bi->get_pos() - cmr;

    if (dr*dv <= 0) return 0;

    real sep = abs(dr);
    real virial_ratio = .5 * dv * dv * sep / (mi + mj + mk);

    if (sep > r_stop) {    // Emergency stop: should really be put elsewhere!

        real ajk = 0.5 * mjk / abs(ejk);

	final.descriptor = stopped;

	final.sma = ajk;
	real ang_mom = abs((bk->get_pos() - bj->get_pos())
			   ^ (bk->get_vel() - bj->get_vel()));

	real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
	if (ecc2 < 1.e-10) {

	    // Take care with small eccentricities! 

	    kepler k;
	    set_kepler_from_sdyn3(k, bj, bk);
	    final.ecc = k.get_eccentricity();

	} else
	    final.ecc = sqrt(ecc2);

	final.escaper = bi->get_index();
	final.virial_ratio = virial_ratio;
	final.outer_separation = sep;

	return 1;
    }

    // Test for sufficiently unperturbed inner binary.
    
    real ajk = 0.5 * mjk / abs(ejk);

    // Determine limiting radius for "negligible" tidal perturbation
    // Note that the term in the denominator is mjk + mi, not just
    // mi, in order to prevent problems with low-mass third stars.

    real rlimit = (ajk + bi->get_radius())
	             * pow(TIDAL_TOL_FACTOR * mjk / (mi + mjk), -1/3.0);

    if (sep < rlimit) return 0;

    // Test for sufficiently non-perturbed outer binary. The safety
    // factor is reduced at large radii.

    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / sep;

//    int p = cerr.precision(STD_PRECISION);
//    cerr << "inner unperturbed, sep, rlim, factor, |v - 1| = " << endl;
//    cerr << "   " << sep <<" "<< rlimit <<" "<< scaled_safety_factor
//	 <<" "<< abs(virial_ratio - 1) << endl;
//    cerr.precision(p);

    if (abs(virial_ratio - 1) < scaled_safety_factor) return 0;

    // Now outer binary is either clearly bound or clearly unbound.

    if (virial_ratio > 1) {

	// Body bi is escaping.

	if (bi->get_index() == 1)	// Initial binary was (1,2).
	    final.descriptor = exchange_1;
	else if (bi->get_index() == 2)
	    final.descriptor = exchange_2;
	else if (bi->get_index() == 3)
	    final.descriptor = preservation;
	else
	    final.descriptor = unknown_final;

	final.sma = ajk;
	real ang_mom = abs((bk->get_pos() - bj->get_pos())
			   ^ (bk->get_vel() - bj->get_vel()));

	real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
	if (ecc2 < 1.e-10) {

	    // Take care with small eccentricities! 

	    kepler k;
	    set_kepler_from_sdyn3(k, bj, bk);
	    final.ecc = k.get_eccentricity();

	} else
	    final.ecc = sqrt(ecc2);

	final.escaper = bi->get_index();
	final.virial_ratio = virial_ratio;
	final.outer_separation = sep;

	return 1;

    } else {

	// Analytically extend the orbits by an integral number
	// of binary periods.

//	real e = energy(bi->get_parent());

	extend_orbits(bj, bk, bi);
	inter.n_kepler++;

//	cerr << "energy error = " << energy(bi->get_parent()) - e << endl;

    }
    
    return 0;
}

// set_merger_mass_and_radius: determine mass and radius of merger product

void set_merger_mass_and_radius(sdyn3 * bn, sdyn3 * bi, sdyn3 * bj)
    {
    bn->set_mass(bi->get_mass() + bj->get_mass());	 // No mass loss
    bn->set_radius(bi->get_radius() + bj->get_radius()); // R \propto M
    }

// set_merger_dyn: determine pos, vel, acc and jerk of merger product,
//                 and determine energy dissipation

void set_merger_dyn(sdyn3 * bn, sdyn3 * bi, sdyn3 * bj)
    {
    real  mi = bi->get_mass();
    real  mj = bj->get_mass();
    real  m_inv = 1/(mi + mj);
    
    bn->set_pos((mi*bi->get_pos() + mj*bj->get_pos())*m_inv);
    bn->set_vel((mi*bi->get_vel() + mj*bj->get_vel())*m_inv);
    bn->set_acc((mi*bi->get_acc() + mj*bj->get_acc())*m_inv);
    bn->set_jerk((mi*bi->get_jerk() + mj*bj->get_jerk())*m_inv);

    vec d_pos = bi->get_pos() - bj->get_pos();
    real rij = sqrt(d_pos*d_pos);

    vec d_vel = bi->get_vel() - bj->get_vel();
    real vij2 = d_vel*d_vel;

    real eij_pot = -mi * mj / rij;
    real eij_kin = 0.5 * mi * mj * m_inv * vij2;

    // Include tidal term from spectator star, if any:

    sdyn3 * b = bi->get_parent();
    sdyn3 * bk = NULL;

    for_all_daughters(sdyn3, b, bb)
	if (bb != bi && bb != bj) bk = bb;

    real tidal_pot = 0;
    if (bk)
	{
	real mn = mi + mj;
	real mk = bk->get_mass();
	tidal_pot = -mn * mk / abs(bn->get_pos() - bk->get_pos())
		    + mi * mk / abs(bi->get_pos() - bk->get_pos())
		    + mj * mk / abs(bj->get_pos() - bk->get_pos());
	}

    bn->set_energy_dissipation(bi->get_energy_dissipation()
			       + bj->get_energy_dissipation()
			       + eij_pot + eij_kin
			       - tidal_pot);

    }

// merge: replace two particles by their center of mass.

void merge(sdyn3 * bi, sdyn3 * bj)
{
    sdyn3 * b = bi->get_parent();
    if (b != bj->get_parent()) err_exit("merge: parent conflict...");

    // Note: any stories attached to the particles are lost.

    sdyn3 * bn = new sdyn3();

//    int p = cerr.precision(STD_PRECISION);
//    cerr << "entering merge(" << bi->get_index() << ", " << bj->get_index()
//         << ") with r1 = " << bi->get_pos() << "\n                     "
//         << " and r2 = " << bj->get_pos()   << endl
//         << "         at t = " << b->get_time() << endl;
//    cerr.precision(p);

    set_merger_mass_and_radius(bn, bi, bj);
    set_merger_dyn(bn, bi, bj);
    bn->set_time(bi->get_time());

    int max_index = 0;
    for_all_daughters(sdyn3, b, bb)
	if (max_index < bb->get_index())
	    max_index = bb->get_index();

    bn->set_label(max_index+1);

    detach_node_from_general_tree(bi);
    detach_node_from_general_tree(bj);
  
    add_node(bn, b);
}

// merge_collisions: recursively merge any stars in contact.

local void merge_collisions(sdyn3 * b)
    {
    int coll = 1;
    while (coll)
	{
	coll = 0;
	for_all_daughters(sdyn3, b, bi)
	    {
	    if (coll) break;
	    for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL; 
		 bj = bj->get_younger_sister())
		{
		if ( abs(bi->get_pos() - bj->get_pos())
		        < (bi->get_radius() + bj->get_radius()) )
		    {
		    merge(bi, bj);
		    coll = 1;
		    break;
		    }
		}
	    }
	}
    }

// triple_escape: return 1 if star 3 is unbound with respect to the center
//                of mass of stars 1 and 2 (distance r3) and is a negligible
//                perturber of the (1,2) hyperbolic orbit.

local int triple_escape(real e12, real r3, real m12, real m3)
{
    if (e12 < 0) return 0; // e12 = energy per unit reduced mass of pair (1,2)

    real a12 = 0.5 * m12 / e12;

    // The m12 in the denominator here includes the tidal effect of (1,2) on 3.

    if (r3 > a12 * pow(TIDAL_TOL_FACTOR * m12 / (m3 + m12), -1/3.0))
        return 1;
    else
        return 0;
}

local int extend_or_end_scatter3(sdyn3 * b, real r_stop,
				 intermediate_state3 & inter,
				 final_state3 & final)
    {

    b->to_com();          // map the whole system to center-of-mass frame

    sdyn3 * b1 = b->get_oldest_daughter();
    sdyn3 * b2 = b1->get_younger_sister();
    sdyn3 * b3 = b2->get_younger_sister();

    real r12 = abs(b1->get_pos() - b2->get_pos());
    real r23 = abs(b2->get_pos() - b3->get_pos());
    real r31 = abs(b3->get_pos() - b1->get_pos());

    real m1 = b1->get_mass();
    real m2 = b2->get_mass();
    real m3 = b3->get_mass();

    real k1 = 0.5 * m1 * b1->get_vel() * b1->get_vel();
    real k2 = 0.5 * m2 * b2->get_vel() * b2->get_vel();
    real k3 = 0.5 * m3 * b3->get_vel() * b3->get_vel();

    real phi12 = m1 * m2 / r12;
    real phi23 = m2 * m3 / r23;
    real phi31 = m3 * m1 / r31;

    real mu12 = m1 * m2 / (m1 + m2);
    real mu23 = m2 * m3 / (m2 + m3);
    real mu31 = m3 * m1 / (m3 + m1);

    vec dv = b2->get_vel() - b1->get_vel();
    real k12 = 0.5 * mu12 * dv * dv;
    real vr12 = dv * (b2->get_pos() - b1->get_pos());

    dv = b3->get_vel() - b2->get_vel();
    real k23 = 0.5 * mu23 * dv * dv;
    real vr23 = dv * (b3->get_pos() - b2->get_pos());

    dv = b1->get_vel() - b3->get_vel();
    real k31 = 0.5 * mu31 * dv * dv;
    real vr31 = dv * (b1->get_pos() - b3->get_pos());

    // First test for ionization.

    if (k1 + k2 + k3 >= phi12 + phi23 + phi31
	&& r12 > LARGE_SEPARATION 
	&& r23 > LARGE_SEPARATION
	&& r31 > LARGE_SEPARATION
	&& k12 > phi12 && k23 > phi23 && k31 > phi31
	&& vr12 > 0 && vr23 > 0 && vr31 > 0
	&& triple_escape(k12 - phi12, min(r23, r31), m1 + m2, m3)
	&& triple_escape(k23 - phi23, min(r31, r12), m2 + m3, m1)
	&& triple_escape(k31 - phi31, min(r12, r23), m3 + m1, m2)) {

	      final.descriptor = ionization;
	      final.virial_ratio = min(min(k12/phi12, k23/phi23), k31/phi31);
	      final.outer_separation = -1;
	      final.escaper = -1;                // -1 means all escaping
	      return 1;
    }

    // Now test the closest pair and the third star for escape.
    // Function "escape" also does analytic extension, if appropriate.
    // Note that the emergency stopping criterion R > r_stop returns
    // looking like an escape has occurred. This should be decoupled from
    // the escape test (someday).

    if (r12 * LARGE_SEPARATION_FACTOR < (r23 + r31))
	return escape(b3, b1, b2, (k12 - phi12)/mu12, r_stop, inter, final);

    if (r23 * LARGE_SEPARATION_FACTOR < (r31 + r12))
	return escape(b1, b2, b3, (k23 - phi23)/mu23, r_stop, inter, final);

    if (r31 * LARGE_SEPARATION_FACTOR < (r12 + r23))
	return escape(b2, b3, b1, (k31 - phi31)/mu31, r_stop, inter, final);

    return 0;
    }

local int extend_or_end_scatter2(sdyn3 * b, final_state3 & final)
{
    sdyn3 * d1 = b->get_oldest_daughter();
    sdyn3 * d2 = d1->get_younger_sister();
    kepler k;

    sdyn3 * merger = NULL;
    sdyn3 * single = NULL;

    if (d1->get_index() > 3)
	{
	merger = d1;
	single = d2;
	}
    else
	{
	merger = d2;
	single = d1;
	}

    set_kepler_from_sdyn3(k, d1, d2);

    real closest_approach = k.get_periastron();
    if (k.get_energy() >= 0
	&& (d1->get_pos() - d2->get_pos())
	    * (d1->get_vel() - d2->get_vel()) > 0)
	closest_approach = k.get_separation();

    if (closest_approach < d1->get_radius() + d2->get_radius()) {

        merge(d1, d2);
	d1 = b->get_oldest_daughter();

	final.descriptor = triple_merger;
	final.escaper = 0;                // 0 means no escaper
	final.outer_separation = 0;
	final.virial_ratio = -1;

    } else {

        real closest_approach_squared = closest_approach * closest_approach;

        merger->set_min_nn_dr2(closest_approach_squared);
        merger->set_min_nn_label(single->get_index());

	if (closest_approach_squared < single->get_min_nn_dr2())
	    {
	    single->set_min_nn_dr2(closest_approach_squared);
	    single->set_min_nn_label(merger->get_index());
	    }


	if (k.get_energy() < 0) {

	    if (single->get_index() == 1)
		final.descriptor = merger_binary_1;
	    else if (single->get_index() == 2)
		final.descriptor = merger_binary_2;
	    else if (single->get_index() == 3)
		final.descriptor = merger_binary_3;
	    else
		final.descriptor = unknown_final;

	    final.sma = k.get_semi_major_axis();
	    final.ecc = k.get_eccentricity();
	    final.outer_separation = -1;
	    final.escaper = 0;
	    final.virial_ratio = -1;

	} else {

	    if (single->get_index() == 1)
		final.descriptor = merger_escape_1;
	    else if (single->get_index() == 2)
		final.descriptor = merger_escape_2;
	    else if (single->get_index() == 3)
		final.descriptor = merger_escape_3;
	    else
		final.descriptor = unknown_final;

	    final.outer_separation = k.get_separation();
	    final.escaper = single->get_index();
	    final.virial_ratio = k.get_energy() * k.get_separation() 
	                                        / k.get_total_mass() - 1;
	}
    }

    return 1;
}

local int extend_or_end_scatter1(final_state3 & final)
{
    final.descriptor = triple_merger;
    final.escaper = 0;                      // 0 denotes no escaper
    final.outer_separation = 0;
    final.virial_ratio = -1;

    return 1;
}

// extend_or_end_scatter:
// 	set the final state and return 1 iff the scattering is over
//      analytically extend the scattering if appropriate

local int extend_or_end_scatter(sdyn3 * b, real r_stop,
				intermediate_state3 & inter,
				final_state3 & final)
{
    // Place an absolute limit on the total number of steps:

    if (b->get_n_steps() > MAX_N_STEPS) {
	final.descriptor = stopped;
	return 1;
    }

    int n = 0;
    for_all_daughters(sdyn3, b, bb) n++;

    // Set some default values for these quantities:

    final.sma = -1;
    final.ecc = -1;

    if (n == 3) 
	return extend_or_end_scatter3(b, r_stop, inter, final);
    else if (n == 2)
	return extend_or_end_scatter2(b, final);
    else if (n == 1)
	return extend_or_end_scatter1(final);
    else {
	cerr << "extend_or_end_scatter: n = " << n << endl;
	exit(1);
    }
    return 0; // To keep HP g++ happy!
}

// check_init: make sure an initial state is sensible.

local void check_init(initial_state3 & init)
    {
    if (init.m2 > 1) err_exit("check_init: m1 < 0");
    if (init.m2 < 0) err_exit("check_init: m2 < 0");
    if (init.m3 < 0) err_exit("check_init: m3 < 0");
    if (init.sma <= 0) err_exit("check_init: semi-major axis <= 0");
    if (init.ecc < 0) err_exit("check_init: eccentricity < 0");
    if (init.ecc > 1) err_exit("check_init: eccentricity > 1");
    if (init.r1 < 0) err_exit("check_init: r1 < 0");
    if (init.r2 < 0) err_exit("check_init: r2 < 0");
    if (init.r3 < 0) err_exit("check_init: r3 < 0");
    if (init.v_inf < 0) err_exit("check_init: v_inf < 0");
    if (init.rho < 0) err_exit("check_init: rho < 0");
    if (init.eta <= 0) err_exit("check_init: eta <= 0");
    }

// sdyn3_to_system: save sdyn3 dynamical information in the specified array

local void sdyn3_to_system(sdyn3 * root, body * system)
{
    int n = 0;
    for_all_daughters(sdyn3, root, bb) {
	system[n].index = bb->get_index();
	system[n].mass  = bb->get_mass();
	for (int kp = 0; kp < 3; kp++) system[n].pos[kp] = bb->get_pos()[kp];
	for (int kv = 0; kv < 3; kv++) system[n].vel[kv] = bb->get_vel()[kv];
	n++;
    }
}

// scatter3: perform a three-body scattering, initializing from a specified
//           state and returning intermediate- and final-state structures.

void scatter3(initial_state3 & init,
	      intermediate_state3 & inter,
	      final_state3 & final,
	      real cpu_time_check,
	      real dt_out,         // diagnostic output interval
	      real dt_snap,        // snapshot output interval
	      real snap_cube_size, // limit output to particles within cube
	      real dt_print,       // print output interval
	      sdyn3_print_fp p)    // function to print output
{
    static int count = 0;

    check_init(init);
    sdyn3 * b = init_to_sdyn3(init, final);


    //--------------------------------------------------


    // Just print out info and return...

    sdyn3* b1 = b->get_oldest_daughter();
    sdyn3* b2 = b1->get_younger_sister();
    sdyn3* b3 = b2->get_younger_sister();
    
    // First look at 1-2 vector:

    vec dx = b1->get_pos() - b2->get_pos();
    cout << ++count
	 << "  " << 180*atan2(dx[1], dx[0])/PI
	 << "  " << dx[2]/abs(dx);

    // Then the outer orbit:

    vec cm = (b1->get_mass()*b1->get_pos()
	           + b2->get_mass()*b2->get_pos())
	             / (b1->get_mass() + b2->get_mass());
    dx = b3->get_pos() - cm;
	
    cout << "  " << 180*atan2(dx[1], dx[0])/PI
	 << "  " << dx[2]/abs(dx);

    cout << endl;


    //--------------------------------------------------


    // Delete the 3-body sytem.

    sdyn3 * bi = b->get_oldest_daughter();
    while (bi) {
	sdyn3 * tmp = bi->get_younger_sister();
	delete bi;
	bi = tmp;
    }
    delete b;
}

main(int argc, char **argv)
    {

    initial_state3 init;
    make_standard_init(init);

    int  seed 	    = 0;    	// seed for random number generator
    int n_rand      = 0;        // number of times to invoke the generator
                                // before starting for real
    int  n_experiments = 1;     // default: only one run
    real dt_out     =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	// output time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int planar_flag = 0;

    bool  b_flag = FALSE;
    bool  q_flag = FALSE;
    bool  Q_flag = FALSE;

    extern char *poptarg;
    int  c;

    while ((c = pgetopt(argc, argv,
		"A:bc:C:d:D:e:L:m:M:n:N:pPqQr:R:s:S:U:v:x:y:z:",
				"$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'b': b_flag = 1 - b_flag;
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      break;
	    case 'D': dt_snap = atof(poptarg);
		      break;
	    case 'e': init.ecc = atof(poptarg);
		      break;
	    case 'L': init.r_init_min = atof(poptarg);
		      break;
	    case 'm': init.m2 = atof(poptarg);
		      break;
	    case 'M': init.m3 = atof(poptarg);
		      break;
	    case 'n': n_experiments = atoi(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'p': planar_flag = 1;
		      break;
	    case 'P': planar_flag = -1;
		      break;
	    case 'q': q_flag = 1 - q_flag;
		      break;
	    case 'Q': Q_flag = 1 - Q_flag;
		      break;
	    case 'r': init.rho = atof(poptarg);
		      break;
	    case 'R': init.r_stop = init.r_init_min
				  = init.r_init_max
				  = atof(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 'S': init.r_stop = atof(poptarg);
		      break;
	    case 'U': init.r_init_max = atof(poptarg);
		      break;
	    case 'v': init.v_inf = atof(poptarg);
		      break;
	    case 'x': init.r1 = atof(poptarg);
		      break;
	    case 'y': init.r2 = atof(poptarg);
		      break;
	    case 'z': init.r3 = atof(poptarg);
		      break;
            case '?': cerr <<
		     "usage: scatter3 [-A #] [-b] [-c #] [-C #] [-d #] [-D #] "
		      << "[-e #] [-L #] [-m #] [-M #] [-n #] [-N #] [-p] "
		      << "[-P] [-q] [-Q] [-r #] [-R #] [-s #] "
		      << "[-S #] [-U #] [-v #] [-x #] [-y #] [-z #]" << endl;
		      exit(1);
	    }            

    if (Q_flag) q_flag = TRUE;

    if (init.m2 > 1)
	{
	cerr << "sigma3: init.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
	}

    cpu_init();
    int random_seed = srandinter(seed, n_rand);

    while (n_experiments--) {

//	cerr << "Random seed = " << get_initial_seed()
//	     << "  n_rand = " << get_n_rand() << flush;

	randomize_angles(init.phase);

	if (planar_flag == 1)
	    init.phase.cos_theta = 1;	// Planar prograde
	else if (planar_flag == -1)
	    init.phase.cos_theta = -1;	// Planar retrograde

	intermediate_state3 inter;
	final_state3 final;

	real cpu = cpu_time();
	scatter3(init, inter, final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size);
    }
}
