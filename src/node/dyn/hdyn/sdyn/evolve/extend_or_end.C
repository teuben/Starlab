
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

#include "scatter.h"
#include "kepler.h"

#define UNPERTURBED_DIAG false

void kepler_pair_to_triple(kepler & k1,	// Inner binary (b1 + b2)
			   kepler & k2,	// Outer binary
			   sdyn * b1,
			   sdyn * b2,
			   sdyn * b3)
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

// tree_is_unbound: return TRUE iff all the top-level nodes of the specified
//                  tree are unbound, widely separated, and receding.
//
// Specifically, require that all top-level pairs be unbound and receding,
// and that no node be a significant perturber of any other pair.
//
// Since deeper levels of structure are already taken into account in the
// tree structure, we only need to check that all nodes at the top level
// are unbound and separating.
//
// Note: If we want to apply analytic extension of nearly unbound top-level
//	 nodes, it should be done here.

bool tree_is_unbound(sdyn* root, real ttf, int debug) {

    if (debug) {
	cerr << "Top level nodes: ";
	for_all_daughters(sdyn, root, bb) cerr << " " << id(bb);
	cerr << endl;
    }

    //    root->flatten_node();		// make the node flat.
    root->to_com();      		// Move to the center-of-mass frame.

    real kin = 0;
    real pot = 0;

    for_all_daughters(sdyn, root, bi) {
	for (sdyn* bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {

	    if (debug) cerr << "checking i = " << id(bi)
			    << "  j = " << id(bj) << endl;
	    
	    // Test radial velocity, separation, and relative energy of (i,j):

	    if ((bi->get_pos() - bj->get_pos())
		 * (bi->get_vel() - bj->get_vel()) < 0) return FALSE;

	    real rij = abs(bi->get_pos() - bj->get_pos());

	    if (debug) cerr << "    rij = " << rij << endl;

	    real mij = bi->get_mass() + bj->get_mass();
	    real mu_scale = ttf * bi->get_mass() / mij;

	    // (mij here for the same reason as in scatter3/triple_escape.)

	    real rlimit = Starlab::max(LARGE_SEPARATION,
	         bi->get_radius() * pow(mu_scale, -1/3.0));

	    if (debug) cerr << "    rlimit = " << rlimit << endl;

	    if (rij < rlimit) return FALSE;

	    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / rij;

	    real kij = 0.5 * square(bi->get_vel() - bj->get_vel());
	    real pij = mij / rij;
	    real eij = kij - pij;

	    if (debug) cerr << "    kij = " << kij
		 << "  pij = " << pij
		 << "  eij = " << kij - pij
		 << endl;

	    if (eij < 0) return FALSE;
	    if (abs(kij/pij - 1) < scaled_safety_factor) return FALSE;

	    real aij = 0.5 * mij / eij;
	    vec cmij = (bi->get_mass() * bi->get_pos()
			     + bj->get_mass() * bj->get_pos()) / mij;

	    // Check the perturbations of all other particles on (i,j).

	    for_all_daughters(sdyn, root, bk)
		if (bk != bi && bk != bj) {

		    real rk = abs(bk->get_pos() - cmij);

		    if (debug) cerr << "    checking perturber " << id(bk)
			 << " at distance " << rk << "...";

		    //PRC(bi->format_label());
		    //PRC(bj->format_label());
		    //PRC(bk->format_label());
		    //PRL(rk/aij);

		    if (rk < aij * pow(ttf * mij
				       / (bk->get_mass() + mij), -1/3.0)) {
			if (debug) cerr << "too close" << endl;
			return FALSE;
		    }
		    if (debug) cerr << endl;
		}

	    if (debug) cerr << "    done" << endl;

	    pot -= bi->get_mass() * bj->get_mass() / rij;
	}
	kin += 0.5 * bi->get_mass() * square(bi->get_vel());
    }

    // Finally, check total energy.

    // cerr << "total system energy = " << kin + pot << endl;

    if (kin + pot <= 0) return FALSE;
    return TRUE;
}

#if 0

bool tree_is_unbound(sdyn* root, int debug) {

    if (debug) {
	cerr << "Top level nodes: ";
	for_all_daughters(sdyn, root, bb) cerr << " " << id(bb);
	cerr << endl;
    }

    root->to_com();      // Move to the center-of-mass frame.

    real kin = 0;
    real pot = 0;

    for_all_daughters(sdyn, root, bi) {
	for (sdyn* bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {

	    if (debug) cerr << "checking i = " << id(bi)
			    << "  j = " << id(bj) << endl;
	    
	    // Test radial velocity, separation, and relative energy of (i,j):

	    if ((bi->get_pos() - bj->get_pos())
		 * (bi->get_vel() - bj->get_vel()) < 0) return FALSE;

	    real rij = abs(bi->get_pos() - bj->get_pos());

	    if (debug) cerr << "    rij = " << rij << endl;

	    real mij = bi->get_mass() + bj->get_mass();
	    real mu_scale = TIDAL_TOL_FACTOR * bi->get_mass() / mij;

	    // (mij here for the same reason as in scatter3/triple_escape.)

	    real rlimit = Starlab::max(LARGE_SEPARATION,
	         bi->get_radius() * pow(mu_scale, -1/3.0));

	    if (debug) cerr << "    rlimit = " << rlimit << endl;

	    if (rij < rlimit) return FALSE;

	    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / rij;

	    real kij = 0.5 * square(bi->get_vel() - bj->get_vel());
	    real pij = mij / rij;
	    real eij = kij - pij;

	    if (debug) cerr << "    kij = " << kij
		 << "  pij = " << pij
		 << "  eij = " << kij - pij
		 << endl;

	    if (eij < 0) return FALSE;
	    if (abs(kij/pij - 1) < scaled_safety_factor) return FALSE;

	    real aij = 0.5 * mij / eij;
	    vec cmij = (bi->get_mass() * bi->get_pos()
			     + bj->get_mass() * bj->get_pos()) / mij;

	    // Check the perturbations of all other particles on (i,j).

	    for_all_daughters(sdyn, root, bk)
		if (bk != bi && bk != bj) {

		    real rk = abs(bk->get_pos() - cmij);

		    if (debug) cerr << "    checking perturber " << id(bk)
			 << " at distance " << rk << "...";

		    if (rk < aij * pow(TIDAL_TOL_FACTOR * mij
				       / (bk->get_mass() + mij), -1/3.0)) {
			if (debug) cerr << "too close" << endl;
			return FALSE;
		    }
		    if (debug) cerr << endl;
		}

	    if (debug) cerr << "    done" << endl;

	    pot -= bi->get_mass() * bj->get_mass() / rij;
	}
	kin += 0.5 * bi->get_mass() * square(bi->get_vel());
    }

    // Finally, check total energy.

    // cerr << "total system energy = " << kin + pot << endl;

    if (kin + pot <= 0) return FALSE;
    return TRUE;
}
#endif

local real potential_energy(sdyn * b) {

    real u = 0;

    for_all_daughters(sdyn, b, bi) {
	for (sdyn * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	  u -= bi->get_mass() * bj->get_mass()
	       / abs(bi->get_pos() - bj->get_pos());
    }

    return u;
}

local real energy(sdyn * b)
{
    real k = 0, u = 0, dissipation = 0;

    for_all_daughters(sdyn, b, bi) {

	// PRC(bi), PRL(bi->get_pos());

	k += bi->get_mass() * bi->get_vel() * bi->get_vel();
	dissipation += bi->get_energy_dissipation();

	for (sdyn * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister()) {

	    // PRC(bj), PRL(bj->get_pos());

	    u -= bi->get_mass() * bj->get_mass()
		  / abs(bi->get_pos() - bj->get_pos());
	}
    }

    return 0.5*k + u + dissipation;
}

void set_kepler_from_sdyn(kepler& k, sdyn* b1, sdyn* b2)
{
    if (b1->get_time() != b2->get_time()) 
	err_exit("set_kepler_from_sdyn: inconsistent times");

    k.set_time(b1->get_time());
    k.set_total_mass( b1->get_mass() + b2->get_mass() );
    k.set_rel_pos( b2->get_pos() - b1->get_pos() );
    k.set_rel_vel( b2->get_vel() - b1->get_vel() );

    k.initialize_from_pos_and_vel();
}


local void stop_integration(sdyn * bi, sdyn * bj, sdyn * bk,
			    real ejk, real mjk,
			    real sep, real virial_ratio) {

    real ajk = 0.5 * mjk / abs(ejk);

    //    final->descriptor = stopped;
    //    final->sma = ajk;

    real ang_mom = abs((bk->get_pos() - bj->get_pos())
		       ^ (bk->get_vel() - bj->get_vel()));

    real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
    
    if (ecc2 < 1.e-10) {

	// Take care with small eccentricities! 

	kepler k;
	set_kepler_from_sdyn(k, bj, bk);
	//	final->ecc = k.get_eccentricity();

    } else
      //	final->ecc = sqrt(ecc2);
      ;

    //    final->escaper = bi->get_index();
    //    final->virial_ratio = virial_ratio;
    //    final->outer_separation = sep;

}


local void extend_orbits(sdyn * b1, sdyn * b2, sdyn * b3, real& apo)

// Analytically extend the orbits of the [[1,2],3] hierarchical system.

{
    if ( (b1->get_time() != b2->get_time())
	 || (b1->get_time() != b3->get_time()) )
	err_exit("extend_orbits: inconsistent times");
    
    real time = b3->get_time();
    
    kepler inner, outer;

    set_kepler_from_sdyn(inner, b1, b2);
    
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

    // Reconstruct the new sdyns:
    
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

local int escape(sdyn * bi, sdyn * bj, sdyn * bk, real ejk,
		 real ttf) {

// Return true iff bi is escaping from the center of mass of bj and bk.
// Perform analytic continuation if appropriate.
// On entry, we already know that i is "widely" separated from j and k.

    if (ejk >= 0) return 0;		// Inner pair not bound.

    real mi = bi->get_mass();
    real mj = bj->get_mass();
    real mk = bk->get_mass();
    real mjk = mj + mk;

    // Center of mass position and velocity of j-k "binary":

    vec cmr = (mj * bj->get_pos() + mk * bk->get_pos()) / mjk;
    vec cmv = (mj * bj->get_vel() + mk * bk->get_vel()) / mjk;

    // Return immediately if the third particle is approaching the binary,
    // unless r_stop < 0 (meaning that the integration is to be terminated
    // if we pass apastron).  Note that NO analytic extension is performed
    // on the initial incoming orbit -- it is assumed that the initial
    // conditions take care of this.

    vec dv = bi->get_vel() - cmv;
    vec dr = bi->get_pos() - cmr;

    // PRL(abs(dr));
    // PRL(abs(dv));
    // PRL(dr*dv);
    // PRL(init.r_stop);

    if (dr*dv <= 0) return 0;

    real sep = abs(dr);
    real virial_ratio = .5 * dv * dv * sep / (mi + mj + mk);

    // Test stopping criterion.

    // PRL(sep);

    if (dr*dv <= 0) {

	cerr << "stopping... \n";
	cerr << "  sep = " << sep << "  dr*dv = " << dr*dv;
      cerr << endl;

	stop_integration(bi, bj, bk, ejk, mjk, sep, virial_ratio);
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
		      * pow(ttf * mjk / (mi + mjk), -1/3.0);
    //    if (init.r_stop < 0) rlimit = Starlab::max(rlimit, abs(init.r_stop));

    // PRL(rlimit);

    if (sep < rlimit) return 0;

    // Test for sufficiently non-perturbed outer binary. The safety
    // factor is reduced at large radii.

    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / sep;

//    cerr.precision(6);
//    cerr << "inner unperturbed, sep, rlim, factor, |v - 1| = " << endl;
//    cerr << "   " << sep <<" "<< rlimit <<" "<< scaled_safety_factor
//	 <<" "<< abs(virial_ratio - 1) << endl;

    // PRL(virial_ratio);

    if (abs(virial_ratio - 1) < scaled_safety_factor) return 0;

    // Now outer binary is either clearly bound or clearly unbound.

    if (virial_ratio > 1) {

	// Body bi is escaping.


	if (false) {
#if 0

	    if (bi->get_index() == 1)	// Initial binary was (1,2).
		final->descriptor = exchange_1;
	    else if (bi->get_index() == 2)
		final->descriptor = exchange_2;
	    else if (bi->get_index() == 3)
		final->descriptor = preservation;
	    else
		final->descriptor = unknown_final;
#endif

	    //	    final->sma = ajk;
	    real ang_mom = abs((bk->get_pos() - bj->get_pos())
			       ^ (bk->get_vel() - bj->get_vel()));

	    real ecc2 = 1 - ang_mom * ang_mom / (mjk * ajk);
	    if (ecc2 < 1.e-10) {

		// Take care with small eccentricities! 

		kepler k;
		set_kepler_from_sdyn(k, bj, bk);
		//		final->ecc = k.get_eccentricity();

	    } else
	      //		final->ecc = sqrt(ecc2);
	      ;

	    //	    final->escaper = bi->get_index();
	    //	    final->virial_ratio = virial_ratio;
	    //	    final->outer_separation = sep;
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

	//	if (inter) inter->n_kepler++;

	if (UNPERTURBED_DIAG)
	    cerr << "Done.  Energy error = "
		 << energy(bi->get_parent()) - e
		 << endl;

	// See if we exceeded the stopping condition during the extension.

	//	if (apo >= init.r_stop) {
	//	    stop_integration(bi, bj, bk, ejk, mjk,
	//			     min(apo, abs(init.r_stop)), virial_ratio);
	//	    return 1;
	//	}
    }
    
    return 0;
}

#define UNPERTURBED_DIAG false

local void determine_triple_parameters(sdyn * b, real& apo) {

  cerr << "determine_triple_parameters" << endl;

    sdyn *b1, *b2, *b3;
    int i=0;
    for_all_leaves(sdyn, b, bi) { 
      if(i==0) b1 = bi;
      else if(i==1) b2 = bi;
      else if(i==2) b3 = bi;
      i++;
    }

    if ( (b1->get_time() != b2->get_time())
	 || (b1->get_time() != b3->get_time()) )
	err_exit("extend_orbits: inconsistent times");
    
    real time = b3->get_time();
    
    kepler inner, outer;

    set_kepler_from_sdyn(inner, b1, b2);
    
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

    // Determine the outer apocenter, to check the stopping criterion
    // in the calling routine.

    apo = outer.get_semi_major_axis() * (1 + outer.get_eccentricity());
    
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


    // Map forward in time by an integral number of inner binary periods
    // to an incoming configuration as close as possible to the outgoing one.

//    real orbits = 2 * (outer.pred_advance_to_apastron() - time)
//				/ inner.get_period();

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

local int extend_or_end_scatter4(sdyn * b) {

    b->flatten_node();
    b->to_com();          // map the whole system to center-of-mass frame
    real mt= b->get_mass();

    real rij[5], phi[5], mu[5], kij[5], vrij[5];
    real mij[5], mkl[5];
    sdyn* bj = NULL;

    int i=0;
    real phit=0;
    vec dv;
    for_all_leaves(sdyn, b, bi) {
      for (bj = (sdyn*)bi->next_node(b); bj != NULL;
	   bj = (sdyn*)bj->next_node(b)) {

	mij[i] = bi->get_mass() + bj->get_mass();
	mkl[i] = mt - mij[i];
	rij[i] = abs(bi->get_pos() - bj->get_pos());
	PRC(rij[i]);PRC(bi->format_label());PRL(bj->format_label());
	phi[i] = bi->get_mass() * bj->get_mass() / rij[i];
	phit += phi[i];
	mu[i]  = bi->get_mass() * bj->get_mass() / (mij[i]);

	dv = bj->get_vel() - bi->get_vel();
	kij[i] = 0.5 * mu[i] * dv * dv;
	vrij[i] = dv * (bj->get_pos() - bi->get_pos());

	i++;
      }
    }

    i=0;
    real m[4], k[4];
    real kt=0;
    {
    for_all_leaves(sdyn, b, bj) {
      m[i]  = bj->get_mass();
      k[i]  = 0.5 * m[i] * bj->get_vel() * bj->get_vel();
      kt += k[i]; 
      i++;
    }}

    bool unbound = false;
    if (kt >= phit)
      unbound = true;

    bool large_distance = false;
    bool receding = false;
    if(!unbound)
      for(i=0; i<5; i++) {
	if(kij[i]>phi[i]) {
	  unbound = true;
	  break;
	}
	if(rij[i]>LARGE_SEPARATION) {
	  large_distance = true;
	  break;
	}
	if(vrij[i]>0) {
	  receding = true;
	  break;
	}
      }

    real rmin;
    bool tescape = false;
    real ttf = TIDAL_TOL_FACTOR; //  1e-6;
    for(i=0; i<5; i++) {
      rmin = VERY_LARGE_NUMBER;
      for(int j=0; j<5; j++) 
	if(j!=i) rmin = Starlab::min(rmin, rij[j]);
      tescape = triple_escape(kij[i] - phi[i], rmin, mij[i], mkl[i], ttf);
      if(tescape) break;
    }

    PRC(unbound);PRC(large_distance);PRC(receding);PRL(tescape);
    // terminate if any of these is true.
    if(unbound || large_distance || receding || tescape)
      return 0;


    // Now test the closest pair and the third star for escape.
    // Function "escape" also does analytic extension, if appropriate.
    // Note that the emergency stopping criterion R > |r_stop| returns
    // looking like an escape has occurred. This should be decoupled from
    // the escape test (someday).

#if 0 
    bool bescape = false;
    real rtot;
    for(i=0; i<5; i++) {
      rtot = 0;
      for(j=0; j<5; j++) 
	if(j!=i) rtot += rij[j];
      if (rij[i] * LARGE_SEPARATION_FACTOR < rtot)
     bescape = escape(b3, b1, b2, (k12 - phi12)/mu12, ttf);
    if(bescape) return 0;
}
#endif

    cerr << "Make tree" << flush << endl;
    make_tree(b, DYNAMICS, STABILITY, K_MAX, true);
    cerr << "tree made" << endl;
    
    return 1;

}

local int extend_or_end_scatter3(sdyn * b) {

    b->to_com();          // map the whole system to center-of-mass frame
    sdyn *b1, *b2, *b3;
    int i=0;
    for_all_leaves(sdyn, b, bi) { 
      if(i==0) b1 = bi;
      else if(i==1) b2 = bi;
      else if(i==2) b3 = bi;
      i++;
    }

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

    real ttf = TIDAL_TOL_FACTOR; //  1e-6;
    if (k1 + k2 + k3 >= phi12 + phi23 + phi31
	&& r12 > LARGE_SEPARATION 
	&& r23 > LARGE_SEPARATION
	&& r31 > LARGE_SEPARATION
	&& k12 > phi12 && k23 > phi23 && k31 > phi31
	&& vr12 > 0 && vr23 > 0 && vr31 > 0
	&& triple_escape(k12 - phi12, Starlab::min(r23, r31), m1 + m2, m3, ttf)
	&& triple_escape(k23 - phi23, Starlab::min(r31, r12), m2 + m3, m1, ttf)
	&& triple_escape(k31 - phi31, Starlab::min(r12, r23), m3 + m1, m2, ttf)) {

	return 1;
    }

    // Now test the closest pair and the third star for escape.
    // Function "escape" also does analytic extension, if appropriate.
    // Note that the emergency stopping criterion R > |r_stop| returns
    // looking like an escape has occurred. This should be decoupled from
    // the escape test (someday).

    if (r12 * LARGE_SEPARATION_FACTOR < (r23 + r31))
      return escape(b3, b1, b2, (k12 - phi12)/mu12, ttf);

    if (r23 * LARGE_SEPARATION_FACTOR < (r31 + r12))
      return escape(b1, b2, b3, (k23 - phi23)/mu23, ttf);

    if (r31 * LARGE_SEPARATION_FACTOR < (r12 + r23))
      return escape(b2, b3, b1, (k31 - phi31)/mu31, ttf);

    return 0;
}

local int extend_or_end_scatter2(sdyn * b) {

    sdyn * d1 = b->get_oldest_daughter();
    sdyn * d2 = d1->get_younger_sister();

    // Consistency (just in case):

    d2->set_younger_sister(NULL);

    kepler k;

    sdyn * merger = NULL;
    sdyn * single = NULL;

    // One of the two stars must be a merger product.

    if (d1->get_index() > 3) {
	merger = d1;
	single = d2;
    } else {
	merger = d2;
	single = d1;
    }

    set_kepler_from_sdyn(k, d1, d2);
    //    k.print_all(cerr);

    real closest_approach = k.get_periastron();
    if (k.get_energy() >= 0
	&& (d1->get_pos() - d2->get_pos())
	    * (d1->get_vel() - d2->get_vel()) > 0)
	closest_approach = k.get_separation();
    else {
      return 2;   // indicate that two-body system is still bound
    }

    if (closest_approach < d1->get_radius() + d2->get_radius()) {
      PRL(closest_approach);
      cerr <<"Merge nodes"<<endl;
      merge(d1, d2);
      d1 = b->get_oldest_daughter();
      d1->set_younger_sister(NULL);	// Again, just in case.

#if 0
	if (final) {
	    final->descriptor = triple_merger;
	    final->escaper = 0;                // 0 means no escaper
	    final->outer_separation = 0;
	    final->virial_ratio = -1;
	}
#endif

    } else {

        real closest_approach_squared = closest_approach * closest_approach;


	merger->set_min_nn_dr2(closest_approach_squared);
	merger->set_min_nn_label(single->get_index());

	if (closest_approach_squared < single->get_min_nn_dr2()) {
	  // cerr <<" Merge NODES?"<<endl;
	  // PRC(closest_approach_squared);PRL(single->get_min_nn_dr2());
	  // PRC(single->format_label());PRL(merger->format_label());
	    single->set_min_nn_dr2(closest_approach_squared);
	    single->set_min_nn_label(merger->get_index());
	}
    }

    return 1;
}

int extend_or_end_scatter(sdyn * b, real ttf, bool debug) {

  if (b->get_n_steps() > MAX_N_STEPS) {
    return 1;
  }

  int unbound = 0;
  if(b->n_leaves()==4) {
    unbound = tree_is_unbound(b, ttf, debug);
  }
  else if(b->n_leaves()==3) {
    unbound = extend_or_end_scatter3(b);
    //    determine_triple_parameters(b, apo);
  } else if(b->n_leaves()==2) {
    unbound = extend_or_end_scatter2(b);
  } else if(b->n_leaves()==1) {
    //          return extend_or_end_scatter1(b);
    unbound = 1;
  }

    return unbound; 
}

