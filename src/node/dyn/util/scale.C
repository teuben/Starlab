 
//// scale:  (re)scale an N-body system to specified M, Q (=T/U), and E.
////         Note that if the required Q implies E > 0, no additional
////         energy scaling is applied.
////
//// Options:    -e    specify softening parameter [0]
////             -E    specify total energy [don't scale]
////             -m/M  specify total mass [don't scale]
////             -q/Q  specify virial ratio [don't scale]
////             -s/S  scale to "standard units [not set]
////
////         NOTE that only top-level nodes are considered in scale_virial
////         and scale_energy.
////
////         If standard units are chosen with "-S", then the center of mass
////         of the system is also set to zero.
////
////         As of 6/99, systems with embedded tidal fields are also properly
////         scaled (assuming they are Roche-lobe filling).

//  Steve McMillan, April 1993

#include "dyn.h"

#ifndef TOOLBOX

#define ALL_i ni = n->get_oldest_daughter(); ni != NULL; \
					     ni = ni->get_younger_sister()
#define j_ABOVE_i nj = ni->get_younger_sister(); nj != NULL; \
					     nj = nj->get_younger_sister()

void scale_mass(dyn* n, real m)
{
    real mass = 0;
    for_all_daughters(dyn, n, ni)
	mass += ni->get_mass();

    real mscale = m/mass;
    for_all_nodes(dyn, n, nj)
	nj->set_mass(mscale*nj->get_mass());	// no "scale_mass"!

    n->set_mass(m);
}

real get_top_level_kinetic_energy(dyn*b)	// top-level nodes only
{
    real kinetic_energy = 0;

    for_all_daughters(dyn, b, bb)
	kinetic_energy += bb->get_mass() * square(bb->get_vel());

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

real get_kinetic_energy(dyn*b)			// all nodes
{
    real kinetic_energy = 0;

    for_all_leaves(dyn, b, bb)
	kinetic_energy += bb->get_mass()
	    * square(something_relative_to_root(bb, &dyn::get_vel));

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

// NOTE: get_top_level_energies does *not* resolve binaries.

void get_top_level_energies(dyn * n, real eps2,
			    real& potential_energy,
			    real& kinetic_energy)
{
//    dyn  * ni, *nj;
// 
//     kinetic_energy = get_top_level_kinetic_energy(n);
// 
//     potential_energy = 0;
//     for (ALL_i) {
//         real dphi = 0;
// 	vector xi = ni->get_pos();
// 	for (j_ABOVE_i) {
// 	    vector xij = nj->get_pos() - xi;
//  	    dphi += nj->get_mass()/sqrt(xij*xij + eps2);
// 	}
// 	potential_energy -= ni->get_mass() * dphi;
//     }

    real total_energy;
    calculate_energies(n, eps2,
		       potential_energy, kinetic_energy, total_energy,
		       true);		// ==> CM approximation
}

real get_tidal_energy(dyn* b, real alpha1, real alpha3)
{
    real tidal_energy = 0;

    // Compute the tidal energy for top-level nodes only (i.e. in the
    // center of mass approximation) in the specified tidal field.

    for_all_daughters(dyn, b, bb) {

	real x = bb->get_pos()[0];
	real z = bb->get_pos()[2];

	real dp = alpha1*x*x + alpha3*z*z;
	tidal_energy += bb->get_mass() * dp;
    }

    tidal_energy *= 0.5;
    return tidal_energy;
}

void scale_virial(dyn * n, real q, real& kinetic_energy, real potential_energy)
{
    // Set the virial ratio by scaling the velocities.

    dyn  * ni;

    if (q > 0) q = -q;  // Only need specify |Q|.

    real vscale = sqrt(q*potential_energy/kinetic_energy);

    for (ALL_i) ni->scale_vel(vscale);

    kinetic_energy = q*potential_energy;
}

void scale_energy(dyn * n, real e, real& energy)
{
    // Set the energy by scaling positions and velocities, keeping
    // the virial ratio fixed.  Note that eps = 0 is implicit.

    dyn  * ni;

    if (energy >= 0) return;
    if (e > 0) e = -e;  // Only need specify |E| and only E < 0 makes sense!

    real xscale = energy/e;
    real vscale = sqrt(1./xscale);

    for (ALL_i) {
	ni->scale_pos(xscale);
	ni->scale_vel(vscale);
    }

    energy = e;
}

#else

#define  FALSE  0
#define  TRUE   1

main(int argc, char ** argv) {
    real m = 0, q = -1, e = 0;
    real eps = 0;
    int m_flag = FALSE;
    int q_flag = FALSE;
    int e_flag = FALSE;
    int eps_flag = FALSE;
    int s_flag = false;
    dyn *root;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "m:M:q:Q:e:E:sS";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'M':
	    case 'm': m_flag = TRUE;
		      m = atof(poptarg);
		      break;
	    case 'Q':
	    case 'q': q_flag = TRUE;
		      q = atof(poptarg);
		      break;
	    case 'E': e_flag = TRUE;
		      e = atof(poptarg);
		      break;
	    case 'e': eps_flag = TRUE;
		      eps = atof(poptarg);
		      break;
	    case 'S':
	    case 's': m_flag = TRUE;
		      q_flag = TRUE;
		      e_flag = TRUE;
	    	      s_flag = true;
		      m = 1;
	    	      q = 0.5;
		      e = -0.25;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    if (m_flag && m <= 0) warning("Specified mass <= 0");
    if (e_flag && e >= 0) warning("Specified energy >= 0");
    if (q_flag && q <  0) warning("Specified virial ratio < 0");

    root = get_dyn(cin);

    // Trust the data in the input snapshot.

    real initial_r_virial = getrq(root->get_log_story(), "initial_rvirial");
    real initial_mass = getrq(root->get_log_story(), "initial_mass");

    // Check for embedded tidal information (only from mk_aniso_king, and
    // scaling is OK only because the cluster is Roche-lobe filling).

    bool has_tidal = false;
    real r_jacobi_over_r_virial
	= getrq(root->get_log_story(), "initial_rtidal_over_rvirial");
    real alpha3_over_alpha1
	= getrq(root->get_log_story(), "alpha3_over_alpha1");

    real alpha1 = 0, initial_r_jacobi = 0, alpha3 = 0;

    if (r_jacobi_over_r_virial > -VERY_LARGE_NUMBER
	&& alpha3_over_alpha1 > -VERY_LARGE_NUMBER) {

	// (Saved initial_mass should be 1 in this case, set by mk_aniso_king.)

	if (initial_mass < 0) {
	    warning("scale: initial_mass not set");
	    initial_mass = 0;
	    for_all_daughters(dyn, root, bb)
		initial_mass += bb->get_mass();
	    
	    putrq(root->get_log_story(), "initial_mass", initial_mass);
	}

	has_tidal = true;

	initial_r_jacobi = r_jacobi_over_r_virial * initial_r_virial;

	// Tidal parameters in consistent units:

	alpha1 = -initial_mass / pow(initial_r_jacobi, 3);
	alpha3 = alpha3_over_alpha1 * alpha1;

	PRI(4);PRC(alpha1);PRL(alpha3);
	// (Note that r_jacobi and alpha1,3 should be OK even if
	//  r_virial will become suspect after setting com below.)

    }

    if (s_flag)			// should this be done routinely?
	root->to_com();		// -- changes the kinetic energy and
    				//    the tidal potential, if any

    if (m_flag)
        scale_mass(root, m);

    // NOTE: Scaling the total mass changes both the total energy and
    //	     the virial ratio.  No attempt is made to preserve either
    //	     in cases where they are not specified on the command line.

    real kinetic_energy = 0, potential_energy = 0, tidal_energy = 0;

    if (q_flag || e_flag) {

	get_top_level_energies(root, eps*eps, potential_energy, kinetic_energy);

	if (has_tidal)
	    tidal_energy = get_tidal_energy(root, alpha1, alpha3);
	    
	if (q_flag)
	    scale_virial(root, q, kinetic_energy,
			 potential_energy + tidal_energy);

	// NOTE: Changing the virial ratio changes the total kinetic energy of
	//	 the system, but preserves the potential energy and total mass.

	real energy = kinetic_energy + potential_energy + tidal_energy;

	if (e_flag) {

	    real fac = e/energy;

	    scale_energy(root, e, energy);

	    // Update all relevant quantities.

	    kinetic_energy *= fac;
	    potential_energy *= fac;
	    tidal_energy *= fac;

	    initial_r_virial /= fac;
	    initial_r_jacobi /= fac;

	    // For Roche-lobe filling systems, this also rescales
	    // alpha1 and alpha3 by fac^3:
	    //
	    //		Utotal	  -->  Utotal * fac
	    //		Utotal	  =    Uc + Ut
	    //		Uc ~ M/R	==> R	       -->  R / fac
	    //		Ut ~ alpha*R^2	==> alpha*R^2  -->  alpha*R^2 * fac
	    //				==> alpha      -->  alpha * fac^3

	    alpha1 *= pow(fac, 3);
	    alpha3 *= pow(fac, 3);

	    if (alpha1 != 0) {
		PRI(4);cerr << "scaled  "; PRC(alpha1);PRL(alpha3);
	    }
	}

	// NOTE: The energy scaling is done on the assumption that eps = 0,
	//	 and preserves the value of the virial ratio in that case.
	//	 The total mass is always left unchanged.

	putrq(root->get_log_story(), "initial_total_energy", energy);
    }

    root->log_history(argc, argv);

    // Update any initial data found in the root log story
    // -- probably best to do this only if system_time = 0.

    if (root->get_system_time() == 0) {

        real mass = 0;
	if (m_flag)
	    mass = m;
	else
	    for_all_daughters(dyn, root, bb)
		mass += bb->get_mass();

	if (m_flag && initial_mass > 0)
	    putrq(root->get_log_story(), "initial_mass", mass);

	if (kinetic_energy <= 0) {
	    get_top_level_energies(root, eps*eps,
				   potential_energy, kinetic_energy);
	    if (has_tidal)
		tidal_energy = get_tidal_energy(root, alpha1, alpha3);
	}

	real r_virial = -0.5 * mass / (potential_energy + tidal_energy);

	if (initial_r_virial > 0 && r_virial != initial_r_virial)
	    putrq(root->get_log_story(), "initial_rvirial", r_virial);

	// Better rewrite the tidal ratio if we recomputed the center of mass
	// -- even if r_jacobi is unchanged, r_virial may be different.

	if (has_tidal && s_flag)
	    putrq(root->get_log_story(), "initial_rtidal_over_rvirial",
		  initial_r_jacobi / r_virial,
		  10);		// (precision as in mk_aniso_king.C)

    }

    put_node(cout, *root);
}

#endif

// end of: scale.C
