 
//// scale:  (re)scale an N-body system to specified M, Q (=T/U), and E.
////         Note that if the required Q implies E > 0, no additional
////         energy scaling is applied.
////
//// Options:    -e    specify softening parameter [0]
////             -E    specify total energy [don't scale]
////             -m/M  specify total mass [don't scale]
////             -q/Q  specify virial ratio [don't scale]
////             -r/R  specify virial radius [don't scale]
////             -s/S  scale to "standard" units [not set]
////
////         NOTE that only top-level nodes are considered in scale_virial
////         and scale_energy.
////
////         The center of mass of the system is set to zero.
////         Option "-s" is equivalent to "-M 1 -R 1 -Q 0.5".
////
////         As of 7/01, systems with embedded tidal or other external fields
////         are also properly scaled (assuming that the Jacobi radius scales
////         with the virial radius).  Note that the virial radius is defined
////         in terms of the *internal* potential energy only.  For now, only
////         allow the energy to be specified if there are no external
////         (non-tidal) fields.

// Steve McMillan, April 1993
//
// Significant changes to both command-line options and internal operation
// by SMcM, July 2001 -- added external fields and redefined virial radius.

// ***  MUST make sure that the definitions of virial radius and virial  ***
// ***  equilibrium are consistent between scale and sys_stats.          ***

#include "dyn.h"

#ifndef TOOLBOX

real get_mass(dyn *b)
{
    real mass = 0;
    for_all_daughters(dyn, b, bb)
	mass += bb->get_mass();
    return mass;
}

void scale_mass(dyn* b, real mscale)
{
    for_all_nodes(dyn, b, bb)
	bb->scale_mass(mscale);
}

void scale_pos(dyn* b, real rscale)
{
    for_all_nodes(dyn, b, bb)
	bb->scale_pos(rscale);
}

void scale_vel(dyn* b, real vscale)
{
    for_all_nodes(dyn, b, bb)
	bb->scale_vel(vscale);
}

real get_top_level_kinetic_energy(dyn *b)	// top-level nodes only
{
    real kinetic_energy = 0;

    for_all_daughters(dyn, b, bb)
	kinetic_energy += bb->get_mass() * square(bb->get_vel());

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

real get_kinetic_energy(dyn *b)			// all nodes
{
    real kinetic_energy = 0;

    for_all_leaves(dyn, b, bb)
	kinetic_energy += bb->get_mass()
	    * square(something_relative_to_root(bb, &dyn::get_vel));

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

// NOTE: get_top_level_energies does *not* resolve binaries.

void get_top_level_energies(dyn *b, real eps2,
			    real& potential_energy,
			    real& kinetic_energy)
{
    real total_energy;
    calculate_energies(b, eps2,
		       potential_energy, kinetic_energy, total_energy,
		       true);		// ==> CM approximation
}

void scale_virial(dyn *b, real q, real potential_energy, real& kinetic_energy)
{
    // Set the virial ratio by scaling the velocities.
    // Also rescale the kinetic energy.

    if (q > 0) q = -q;  	// only need specify |Q|

    real vscale = sqrt(q*potential_energy/kinetic_energy);
    scale_vel(b, vscale);
    kinetic_energy = q*potential_energy;
}

real scale_energy(dyn * b, real e, real& energy)
{
    // Set the energy by scaling positions and velocities, keeping
    // the virial ratio fixed.  Note that eps = 0 is implicit.

    if (energy >= 0) return 1;
    if (e > 0) e = -e;  // only need specify |E| and only E < 0 makes sense!

    real xscale = energy/e;
    real vscale = sqrt(1./xscale);

    for_all_daughters(dyn, b, bi) {
	bi->scale_pos(xscale);
	bi->scale_vel(vscale);
    }

    energy = e;
    return xscale;
}

void scale(dyn *b, real eps,
	   bool e_flag, real e,
	   bool m_flag, real m,
	   bool q_flag, real q,
	   bool r_flag, real r,
	   void (*top_level_energies)(dyn*, real,    // default =
				      real&, real&)) // get_top_level_energies()
{
    // Another consistency check:

    if (q_flag && find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {
	warning("scale: can't reset virial ratio for aniso_king model");
	q_flag = false;
    }

    // Always transform to the center of mass frame.

    b->to_com();

    // Define various relevant quantities.  Trust the data in the input
    // snapshot, if current.

    // Always need to know the mass.

    real mass = 0;

    if (b->get_system_time() == 0
	&& find_qmatch(b->get_log_story(), "initial_mass"))
	mass = getrq(b->get_log_story(), "initial_mass");
    else
	mass = get_mass(b);


    // Need to know some energies if the E, Q, or R flags are set.
    // Note that the definition of r_virial now includes *only* the
    // internal potential energy.

    real r_virial = 0, pot_int = 0, pot_ext = 0, kin = 0;

    if (e_flag || r_flag || q_flag) {
	if (b->get_system_time() == 0
	    && find_qmatch(b->get_log_story(), "initial_rvirial")) {

	    // Avoid N^2 calculation if possible.

	    r_virial = getrq(b->get_log_story(), "initial_rvirial");
	    pot_int = -0.5*mass*mass/r_virial;
	    kin = get_top_level_kinetic_energy(b);

#if 0
	    // Check:

	    real pe, ke;
	    top_level_energies(b, eps*eps, pe, ke);
	    PRC(pot_int); PRL(pe);
	    PRC(kin); PRL(ke);
#endif

	} else {

	    // Compute the potential energy.

	    top_level_energies(b, eps*eps, pot_int, kin);
	    r_virial = -0.5*mass*mass/pot_int;
	}

	// Get external potential parameters (including any tidal field)
	// from the input data and compute the external potential energy
	// (excluding any tidal field).

	check_set_external(b);
	pot_ext = get_external_pot(b);
	// PRL(pot_ext);

	// Variable kira_initial_jacobi_radius is set by check_set_external,
	// but it will be rescaled, so just remove it.

	rmq(b->get_log_story(), "kira_initial_jacobi_radius");

	// If an external (non-tidal) field exists, we probably won't
	// naturally want to specify the energy.  In addition, the
	// procedure for doing this is complicated (iterative, and may
	// not converge).  For now, at least, only allow e to be set
	// in the case of no external fields.

	if (e_flag && pot_ext != 0)
	    err_exit("Can't set energy in case of external field.");
    }

    // First do the mass scaling, if any.  NOTE: Scaling the total mass
    // changes both the total energy and the virial ratio.  No attempt is
    // made to preserve either in cases where they are not specified on
    // the command line.

    if (m_flag) {
        real mfac = m/mass;
	// PRL(mfac);

	scale_mass(b, mfac);

	mass = m;
	pot_int *= mfac*mfac;
	pot_ext *= mfac;
	kin *= mfac;

	cerr << "scale:  "; PRL(mass);
    }

    // Simplest choice now is r_flag; e_flag is more complicated,
    // particularly in the presence of an external field.

    if (r_flag) {

	// Rescale all positions to force the virial radius to the
	// specified value.

	real oldvir = pot_int + get_external_virial(b);	  // denominator of
							  // virial ratio

	real rfac =  r/r_virial;
	// PRL(rfac);

	scale_pos(b, rfac);

	r_virial = r;
	pot_int /= rfac;
	pot_ext = get_external_pot(b);

	// Rescale all velocities to preserve the virial ratio.
	// Scaling is already OK in case of a tidal field, so
	// get_external_virial() includes only non-tidal terms.

	real vfac = sqrt((pot_int + get_external_virial(b))/oldvir);
	// PRL(vfac);

	scale_vel(b, vfac);
	kin *= vfac*vfac;

	cerr << "scale:  "; PRL(r_virial);

#if 0
	// Check:

	real pe, ke;
	top_level_energies(b, eps*eps, pe, ke);
	real pext = get_external_pot(b);
	real vir = get_external_virial(b);

	PRC(pot_int); PRL(pe);
	PRC(kin); PRL(ke);
	PRL(ke/(pe+vir));
	PRL(-0.5*mass*mass/pe);
#endif
    }

    if (e_flag) {

	// Attempt to set the energy while holding the virial ratio
	// constant, or else set both e and q.  Assume no external
	// fields (see above).  Procedure is OK as is in the presence
	// of a tidal field.

	if (q_flag) {
	    scale_virial(b, q, pot_int, kin);	// scales kin

	    // NOTE: Changing the virial ratio changes the total kinetic
	    // energy of the system, but preserves the potential energy
	    // and total mass.

	    // kin = q*pot_int;
	    q_flag = false;

	    cerr << "scale:  "; PRL(q);
	}

	// NOTE: The energy scaling is done on the assumption that eps = 0,
	//	 and preserves the value of the virial ratio in that case.
	//	 The total mass is always left unchanged.

	real ee = kin+pot_int;
	real fac = scale_energy(b, e, ee);

	// Update all relevant quantities.

	kin /= fac;
	pot_int /= fac;
	r_virial *= fac;

	// For Roche-lobe filling systems, this also rescales
	// alpha1 and alpha3 by fac^{-3}.

	cerr << "scale:  "; PRL(e);
    }

    if (q_flag) {

	// Scale the velocities to set q and preserve r_virial.

	real vir = get_external_virial(b);
	real qvir = -kin/(pot_int + vir);
	real vfac = sqrt(q/qvir);
	// PRC(vir); PRC(qvir); PRL(vfac);

	scale_vel(b, vfac);
	kin *= vfac*vfac;

	cerr << "scale:  "; PRL(q);
    }

    // Update the root log story -- probably best to do this only if
    // system_time = 0.

    if (b->get_system_time() == 0) {
	putrq(b->get_log_story(), "initial_mass", mass, HIGH_PRECISION);
	putrq(b->get_log_story(), "initial_rvirial", r_virial);
	putrq(b->get_dyn_story(), "initial_total_energy", kin+pot_int+pot_ext);
    }
}

bool parse_scale_main(int argc, char *argv[],
		      real& eps,
		      bool& e_flag, real& e,
		      bool& m_flag, real& m,
		      bool& q_flag, real& q,
		      bool& r_flag, real& r)
{
    // Defaults:

    e_flag = false;
    m_flag = false;
    q_flag = false;
    r_flag = false;
    bool eps_flag = false;
    bool s_flag = false;

    m = 0, q = -1, e = 0, r = 0;
    eps = 0;

    extern char *poptarg;
    int c;
    char* param_string = "m:M:q:Q:e:E:r:R:sS";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'E': e_flag = true;
	    	      r_flag = false;
		      e = atof(poptarg);
		      break;
	    case 'e': eps_flag = true;
		      eps = atof(poptarg);
		      break;
	    case 'M':
	    case 'm': m_flag = true;
		      m = atof(poptarg);
		      break;
	    case 'Q':
	    case 'q': q_flag = true;
		      q = atof(poptarg);
		      break;
	    case 'R':
	    case 'r': r_flag = true;
	    	      e_flag = false;
	    	      r = atof(poptarg);
		      break;
	    case 'S':
	    case 's': s_flag = true;
		      m_flag = true;
		      r_flag = true;
		      q_flag = true;
		      m = 1;
	    	      q = 0.5;
		      r = 1;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      return false;
	}

    if (e_flag && e >= 0) warning("scale: specified energy >= 0");
    if (m_flag && m <= 0) warning("scale: specified mass <= 0");
    if (q_flag && q <  0) warning("scale: specified virial ratio < 0");
    if (r_flag && r <= 0) warning("scale: specified virial radius <= 0");

    return true;
}

#else

main(int argc, char ** argv)
{
    real m = 0, q = -1, e = 0, r = 0;
    real eps = 0;

    bool e_flag = false;
    bool m_flag = false;
    bool q_flag = false;
    bool r_flag = false;

    if (!parse_scale_main(argc, argv, eps,
			  e_flag, e, m_flag, m,
			  q_flag, q, r_flag, r)) {
	get_help();
	exit(1);
    }

    dyn *b = get_dyn(cin);
    b->log_history(argc, argv);

    scale(b, eps, e_flag, e, m_flag, m, q_flag, q, r_flag, r);

    put_node(cout, *b);
}

#endif
