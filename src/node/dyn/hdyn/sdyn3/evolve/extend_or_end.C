
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Code shared by scatter3 and hier3 relating to possible termination
// or analytic extension of motion, or merging two or all stars.
//
// Globally visible functions:
//
//	merge_collisions
//	extend_or_end_scatter

#include "scatter3.h"

#define UNPERTURBED_DIAG false

local void extend_orbits(sdyn3 * b1, sdyn3 * b2, sdyn3 * b3, real& apo)

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
    
    vector cmr = (m1 * b1->get_pos() + m2 * b2->get_pos()) / m12;
    vector cmv = (m1 * b1->get_vel() + m2 * b2->get_vel()) / m12;
    
    outer.set_rel_pos(b3->get_pos() - cmr);
    outer.set_rel_vel(b3->get_vel() - cmv);

    outer.initialize_from_pos_and_vel();

    // Determine the outer apocenter, to check the stopping criterion
    // in the calling routine.

    apo = outer.get_semi_major_axis() * (1 + outer.get_eccentricity());
    
    if (UNPERTURBED_DIAG) {

	cerr << "  inner binary:  r = " << inner.get_separation()
	     << "  a = " << inner.get_semi_major_axis()
	     << "  e = " << inner.get_eccentricity()
	     << endl
	     << "  mass = " << inner.get_total_mass() << endl
	     << "  pos = " << inner.get_rel_pos() << endl
	     << "  vel = " << inner.get_rel_vel() << endl;
	cerr << "  outer binary:  R = " << outer.get_separation()
	     << "  a = " << outer.get_semi_major_axis() << endl
	     << "  e = " << outer.get_eccentricity()
	     << endl
	     << "  mass = " << outer.get_total_mass() << endl
	     << "  pos = " << outer.get_rel_pos() << endl
	     << "  vel = " << outer.get_rel_vel() << endl;

    }

    // Map forward in time by an integral number of inner binary periods
    // to an incoming configuration as close as possible to the outgoing one.

    real orbits = 2 * (outer.pred_advance_to_apastron() - time)
				/ inner.get_period();

    // Note: orbits may be too big to fit in an int...

    real time_0 = time;

    real big = 1.e9;
    while (orbits > big) {
	time += big*inner.get_period();
	orbits -= big;
    }

    int inner_orbits = (int) (orbits + 0.5);
    time += inner_orbits * inner.get_period();

    outer.transform_to_time(time);
    
    // Note that no modification is needed to the inner orbit
    // except to increase its time.
    
    inner.set_time(time);
    
    if (UNPERTURBED_DIAG) {

	cerr << "  outer binary:  R = " << outer.get_separation()
	     << "  a = " << outer.get_semi_major_axis() << endl
	     << "  e = " << outer.get_eccentricity()
	     << endl
	     << "  mass = " << outer.get_total_mass() << endl
	     << "  pos = " << outer.get_rel_pos() << endl
	     << "  vel = " << outer.get_rel_vel() << endl;

    }

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

local void stop_integration(sdyn3 * bi, sdyn3 * bj, sdyn3 * bk,
			    real ejk, real mjk,
			    real sep, real virial_ratio,
			    final_state3 * final = NULL)
{
    if (!final) return;

    real ajk = 0.5 * mjk / abs(ejk);

    final->descriptor = stopped;
    final->sma = ajk;

    real ang_mom = abs((bk->get_pos() - bj->get_pos())
		       ^ (bk->get_vel() - bj->get_vel()));

    real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
    
    if (ecc2 < 1.e-10) {

	// Take care with small eccentricities! 

	kepler k;
	set_kepler_from_sdyn3(k, bj, bk);
	final->ecc = k.get_eccentricity();

    } else
	final->ecc = sqrt(ecc2);

    final->escaper = bi->get_index();
    final->virial_ratio = virial_ratio;
    final->outer_separation = sep;
}

local int escape(sdyn3 * bi, sdyn3 * bj, sdyn3 * bk,
		 real ejk, 				//(for convenience)
		 initial_state3& init,
		 intermediate_state3 * inter = NULL,
		 final_state3 * final = NULL)

// Return true iff bi is escaping from the center of mass of bj and bk.
// Perform analytic continuation if appropriate.
// On entry, we already know that i is "widely" separated from j and k.

{
    if (final) final->descriptor = unknown_final;

    if (ejk >= 0) return 0;		// Inner pair not bound.

    real mi = bi->get_mass();
    real mj = bj->get_mass();
    real mk = bk->get_mass();
    real mjk = mj + mk;

    // Center of mass position and velocity of j-k "binary":

    vector cmr = (mj * bj->get_pos() + mk * bk->get_pos()) / mjk;
    vector cmv = (mj * bj->get_vel() + mk * bk->get_vel()) / mjk;

    // Return immediately if the third particle is approaching the binary,
    // unless r_stop < 0 (meaning that the integration is to be terminated
    // if we pass apastron).  Note that NO analytic extension is performed
    // on the initial incoming orbit -- it is assumed that the initial
    // conditions take care of this.

    vector dv = bi->get_vel() - cmv;
    vector dr = bi->get_pos() - cmr;

    // PRL(abs(dr));
    // PRL(abs(dv));
    // PRL(dr*dv);
    // PRL(init.r_stop);

    if (dr*dv <= 0 && init.r_stop > 0) return 0;

    real sep = abs(dr);
    real virial_ratio = .5 * dv * dv * sep / (mi + mj + mk);

    // Test stopping criterion.

    // PRL(sep);

    if (sep > abs(init.r_stop)
	|| (inter && init.r_stop < 0 && dr*dv <= 0 && inter->n_osc > 0)) {

//	cerr << "stopping... \n";
//	cerr << "  sep = " << sep << "  r_stop = " << init.r_stop
//	     << "  dr*dv = " << dr*dv;
//      if (inter) cerr << "  n_osc = " << inter->n_osc;
//      cerr << endl;

	stop_integration(bi, bj, bk, ejk, mjk,
			 sep, virial_ratio, final);
	return 1;
    }

    // Test for sufficiently unperturbed inner binary.
    
    real ajk = 0.5 * mjk / abs(ejk);

    // PRL(ajk);

    // Determine limiting radius for "negligible" tidal perturbation
    // Note that the term in the denominator is mjk + mi, not just
    // mi, in order to prevent problems with low-mass third stars.

    // The use of r_stop is a kludge to allow "flyby" experiments where
    // we catch the outcoming particle at some specific radius.  For the
    // strictest flyby experiments (see flyby3.C), r_stop < 0, meaning
    // we stop at |r_stop| or apocenter.  In this case ONLY, let r_stop
    // override the tidal limit.  If r_stop > 0, assume that the built-in
    // tidal parameters are OK.

    real rlimit = (ajk + bi->get_radius())
		      * pow(init.tidal_tol_factor * mjk / (mi + mjk), -1/3.0);
    if (init.r_stop < 0) rlimit = max(rlimit, abs(init.r_stop));

    // PRL(rlimit);

    if (sep < rlimit) return 0;

    // Test for sufficiently non-perturbed outer binary. The safety
    // factor is reduced at large radii.

    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / sep;

//    int p = cerr.precision(STD_PRECISION);
//    cerr << "inner unperturbed, sep, rlim, factor, |v - 1| = " << endl;
//    cerr << "   " << sep <<" "<< rlimit <<" "<< scaled_safety_factor
//	 <<" "<< abs(virial_ratio - 1) << endl;
//    cerr.precision(p);

    // PRL(virial_ratio);

    if (abs(virial_ratio - 1) < scaled_safety_factor) return 0;

    // Now outer binary is either clearly bound or clearly unbound.

    if (virial_ratio > 1) {

	// Body bi is escaping.

	if (final) {

	    if (bi->get_index() == 1)	// Initial binary was (1,2).
		final->descriptor = exchange_1;
	    else if (bi->get_index() == 2)
		final->descriptor = exchange_2;
	    else if (bi->get_index() == 3)
		final->descriptor = preservation;
	    else
		final->descriptor = unknown_final;

	    final->sma = ajk;
	    real ang_mom = abs((bk->get_pos() - bj->get_pos())
			       ^ (bk->get_vel() - bj->get_vel()));

	    real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
	    if (ecc2 < 1.e-10) {

		// Take care with small eccentricities! 

		kepler k;
		set_kepler_from_sdyn3(k, bj, bk);
		final->ecc = k.get_eccentricity();

	    } else
		final->ecc = sqrt(ecc2);

	    final->escaper = bi->get_index();
	    final->virial_ratio = virial_ratio;
	    final->outer_separation = sep;
	}

	return 1;

    } else {

	// Analytically extend the orbits by an integral number
	// of binary periods.

	real e = energy(bi->get_parent());

	if (UNPERTURBED_DIAG) 
	    cerr << "Beginning analytic extension at time "
		 << bi->get_parent()->get_time()
		           + bi->get_parent()->get_time_offset()
		 << endl;

	real apo;
	extend_orbits(bj, bk, bi, apo);

	if (inter) inter->n_kepler++;

	if (UNPERTURBED_DIAG)
	    cerr << "Done.  Energy error = "
		 << energy(bi->get_parent()) - e
		 << endl;

	// See if we exceeded the stopping condition during the extension.

	if (apo >= init.r_stop) {
	    stop_integration(bi, bj, bk, ejk, mjk,
			     min(apo, abs(init.r_stop)), virial_ratio, final);
	    return 1;
	}
    }
    
    return 0;
}

// set_merger_mass_and_radius: determine mass and radius of merger product

local void set_merger_mass_and_radius(sdyn3 * bn, sdyn3 * bi, sdyn3 * bj)
{
    bn->set_mass(bi->get_mass() + bj->get_mass());	 // No mass loss
    bn->set_radius(bi->get_radius() + bj->get_radius()); // R \propto M
}

// set_merger_dyn: determine pos, vel, acc and jerk of merger product,
//                 and determine energy dissipation

local void set_merger_dyn(sdyn3 * bn, sdyn3 * bi, sdyn3 * bj)
{
    real  mi = bi->get_mass();
    real  mj = bj->get_mass();
    real  m_inv = 1/(mi + mj);
    
    bn->set_pos((mi*bi->get_pos() + mj*bj->get_pos())*m_inv);
    bn->set_vel((mi*bi->get_vel() + mj*bj->get_vel())*m_inv);
    bn->set_acc((mi*bi->get_acc() + mj*bj->get_acc())*m_inv);
    bn->set_jerk((mi*bi->get_jerk() + mj*bj->get_jerk())*m_inv);

    vector d_pos = bi->get_pos() - bj->get_pos();
    real rij = sqrt(d_pos*d_pos);

    vector d_vel = bi->get_vel() - bj->get_vel();
    real vij2 = d_vel*d_vel;

    real eij_pot = -mi * mj / rij;
    real eij_kin = 0.5 * mi * mj * m_inv * vij2;

    // Include tidal term from spectator star, if any:

    sdyn3 * b = bi->get_parent();
    sdyn3 * bk = NULL;

    for_all_daughters(sdyn3, b, bb)
	if (bb != bi && bb != bj) bk = bb;

    real tidal_pot = 0;
    if (bk) {
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

local void merge(sdyn3 * bi, sdyn3 * bj)
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

    // Construct the new index as one more than the maximum existing index.

//    int max_index = 0;
//    for_all_daughters(sdyn3, b, bb)
//	if (max_index < bb->get_index())
//	    max_index = bb->get_index();
//
//    bn->set_label(max_index+1);

    // Construct the new index by adding the component indices,
    // then adding 1 if neither component is a merger product
    // (so 1 + 2 --> 4, 1 + 3 --> 5, 2 + 3 --> 6, 1 + 2 + 3 --> 7.

    int new_index = bi->get_index() + bj->get_index();
    if (bi->get_index() < 4 && bj->get_index() < 4) new_index++;
    bn->set_label(new_index);

    detach_node_from_general_tree(*bi);
    detach_node_from_general_tree(*bj);
  
    add_node(*bn, *b);
}

// merge_collisions: recursively merge any stars in contact.

void merge_collisions(sdyn3 * b)
{
    int coll = 1;
    while (coll) {

	coll = 0;
	for_all_daughters(sdyn3, b, bi)
	    {
	    if (coll) break;
	    for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL; 
		 bj = bj->get_younger_sister()) {

		if ( abs(bi->get_pos() - bj->get_pos())
		      < (bi->get_radius() + bj->get_radius()) ) {

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

local int triple_escape(real e12, real r3, real m12, real m3,
			real tidal_tol_factor)
{
    if (e12 < 0) return 0; // e12 = energy per unit reduced mass of pair (1,2)

    real a12 = 0.5 * m12 / e12;

    //----------------------------------------------------------------------
    // Hmmm... Still a problem with pow on RH 5.2?

    // if (r3 > a12 * pow(tidal_tol_factor * m12 / (m3 + m12), -1/3.0))

    // The m12 in the denominator here includes the tidal effect of (1,2) on 3.

    real tid_fac = tidal_tol_factor * m12 / (m3 + m12);
    real ratio = a12/r3;

    if (ratio*ratio*ratio > tid_fac)
    //----------------------------------------------------------------------
        return 1;
    else
        return 0;
}

local int extend_or_end_scatter3(sdyn3 * b,
				 initial_state3& init,
				 intermediate_state3 * inter = NULL,
				 final_state3 * final = NULL)
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

    vector dv = b2->get_vel() - b1->get_vel();
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
	&& triple_escape(k12 - phi12, min(r23, r31), m1 + m2, m3,
			 init.tidal_tol_factor)
	&& triple_escape(k23 - phi23, min(r31, r12), m2 + m3, m1,
			 init.tidal_tol_factor)
	&& triple_escape(k31 - phi31, min(r12, r23), m3 + m1, m2,
			 init.tidal_tol_factor)) {

	if (final) {
	    final->descriptor = ionization;
	    final->virial_ratio = min(min(k12/phi12, k23/phi23), k31/phi31);
	    final->outer_separation = -1;
	    final->escaper = -1;                // -1 means all escaping
	}

	return 1;
    }

    // Now test the closest pair and the third star for escape.
    // Function "escape" also does analytic extension, if appropriate.
    // Note that the emergency stopping criterion R > |r_stop| returns
    // looking like an escape has occurred. This should be decoupled from
    // the escape test (someday).

    if (r12 * LARGE_SEPARATION_FACTOR < (r23 + r31))
	return escape(b3, b1, b2, (k12 - phi12)/mu12, init, inter, final);

    if (r23 * LARGE_SEPARATION_FACTOR < (r31 + r12))
	return escape(b1, b2, b3, (k23 - phi23)/mu23, init, inter, final);

    if (r31 * LARGE_SEPARATION_FACTOR < (r12 + r23))
	return escape(b2, b3, b1, (k31 - phi31)/mu31, init, inter, final);

    return 0;
}

local int extend_or_end_scatter2(sdyn3 * b, final_state3 * final = NULL)
{
    sdyn3 * d1 = b->get_oldest_daughter();
    sdyn3 * d2 = d1->get_younger_sister();

    // Consistency (just in case):

    d2->set_younger_sister(NULL);

    kepler k;

    sdyn3 * merger = NULL;
    sdyn3 * single = NULL;

    // One of the two stars must be a merger product.

    if (d1->get_index() > 3) {
	merger = d1;
	single = d2;
    } else {
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
	d1->set_younger_sister(NULL);	// Again, just in case.

	if (final) {
	    final->descriptor = triple_merger;
	    final->escaper = 0;                // 0 means no escaper
	    final->outer_separation = 0;
	    final->virial_ratio = -1;
	}

    } else {

        real closest_approach_squared = closest_approach * closest_approach;

        merger->set_min_nn_dr2(closest_approach_squared);
        merger->set_min_nn_label(single->get_index());

	if (closest_approach_squared < single->get_min_nn_dr2()) {
	    single->set_min_nn_dr2(closest_approach_squared);
	    single->set_min_nn_label(merger->get_index());
	}


	if (final) {

	    if (k.get_energy() < 0) {

		if (single->get_index() == 1)
		    final->descriptor = merger_binary_1;
		else if (single->get_index() == 2)
		    final->descriptor = merger_binary_2;
		else if (single->get_index() == 3)
		    final->descriptor = merger_binary_3;
		else
		    final->descriptor = unknown_final;

		final->sma = k.get_semi_major_axis();
		final->ecc = k.get_eccentricity();
		final->outer_separation = -1;
		final->escaper = 0;
		final->virial_ratio = -1;

	    } else {

		if (single->get_index() == 1)
		    final->descriptor = merger_escape_1;
		else if (single->get_index() == 2)
		    final->descriptor = merger_escape_2;
		else if (single->get_index() == 3)
		    final->descriptor = merger_escape_3;
		else
		    final->descriptor = unknown_final;

		final->outer_separation = k.get_separation();
		final->escaper = single->get_index();
		final->virial_ratio = k.get_energy() * k.get_separation() 
		    				     / k.get_total_mass() - 1;
	    }
	}
    }

    return 1;
}

local int extend_or_end_scatter1(final_state3 * final = NULL)
{
    if (final) {
	final->descriptor = triple_merger;
	final->escaper = 0;                      // 0 denotes no escaper
	final->outer_separation = 0;
	final->virial_ratio = -1;
    }

    return 1;
}

// extend_or_end_scatter:
// 	set the final state and return 1 iff the scattering is over
//      analytically extend the motion if appropriate

int extend_or_end_scatter(sdyn3 * b,
			  initial_state3& init,
			  intermediate_state3 * inter,	// default = NULL
			  final_state3 * final)		// default = NULL
{
    // Place an absolute limit on the total number of steps:

    if (init.r_stop == VERY_LARGE_NUMBER && b->get_n_steps() > MAX_N_STEPS) {
	if (final) final->descriptor = stopped;
	return 1;
    }

    int n = 0;
    for_all_daughters(sdyn3, b, bb) n++;

    // Set some default values for these quantities:

    if (final) {
	final->sma = -1;
	final->ecc = -1;
    }

    if (n == 3) 
	return extend_or_end_scatter3(b, init, inter, final);
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
