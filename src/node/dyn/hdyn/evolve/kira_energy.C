
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions associated with calculating the energy of the system.
//
// Externally visible functions:
//
//	void print_recalculated_energies

#include "hdyn.h"

#define MAX_EKIN	1.0e15		// Assume a hardware error if ekin
					// is greater than this value
					// -- beware of runaway neutron
					// stars!

local void inc_pot(dyn *b, real dpot)	// allow dyn function to access set_pot
{
    hdyn *bb = (hdynptr)b;
    real pot = bb->get_pot();
    bb->set_pot(pot+dpot);
}

void calculate_energies_with_external(hdyn* b,
				      real& epot, real& ekin, real& etot,
				      bool cm,		// default = false
				      bool use_grape)	// default = true
{
    // Compute the total energy, including external terms; also compute
    // the "pot" class data.

    // First, get the internal energy (use GRAPE if available).

    if (!b->get_ignore_internal())
	kira_calculate_internal_energies(b, epot, ekin, etot, cm, use_grape);

    if (b->get_external_field() > 0) {

	// Add the external contribution to the total potential, and
	// add external terms to hdyn::pot of all top-level nodes.

	real dpot = get_external_pot(b, inc_pot);
	epot += dpot;
	etot += dpot;
    }
}

void print_recalculated_energies(hdyn* b,
				 bool print_dde,	// default = false
				 bool save_story)	// default = false
{
    static real de_prev = VERY_LARGE_NUMBER;

    real epot = 0;
    real ekin = 0;
    real etot = 0;

    calculate_energies_with_external(b, epot, ekin, etot);

    // Here, epot includes both external and internal potential energy
    // terms, and ekin includes the kinetic energy of the center of mass.

    int p = cerr.precision(INT_PRECISION);
    cerr << "Energies: " << epot << " " << ekin << " " << etot;

#if 0

    // Recompute on the front end to check the GRAPE calculation
    // (recall that the GRAPE-4 computes to ~single precision, and
    // effectively flattens all tree structures prior to determining
    // the energy, so close binaries may be very poorly handled).

    calculate_energies_with_external(b, epot, ekin, etot, false, false);
    cerr << endl << "Recomputed: " << epot << " " << ekin << " " << etot;

#endif

    cerr.precision(p);

    // Possible hardware problem?

    if (abs(ekin) > MAX_EKIN) {
	cerr << endl << "...retrying..." << endl;

	epot = ekin = etot = 0;
	calculate_energies_with_external(b, epot, ekin, etot);
	cerr << "Energies: " << epot << " " << ekin << " " << etot;
    }

    if (abs(ekin) > MAX_EKIN) print_dde = false;

    real de = 0;

    if (b->get_kira_counters()->initial_etot == 0) {

	b->get_kira_counters()->initial_etot = etot;
	b->get_kira_counters()->de_total = 0;

    } else {

	// real de = (etot - e_corr - initial_etot) / epot;
	// The above normalization fails if many stars have escaped
	// the system with HUGE kicks.

	de = (etot - b->get_kira_counters()->de_total
		- b->get_kira_counters()->initial_etot);
					//   / kira_stats.initial_etot;
					//     ^^^^^^^^ note! ^^^^^^^^
    }

    if (print_dde && de_prev != VERY_LARGE_NUMBER) {
	cerr << endl << "          de = " << de
	     << "  d(de) = " << de - de_prev;
	if (ekin > 0) cerr << "  (" << (de - de_prev) / ekin << ")";
	cerr << endl;
    } else
	cerr << "  de = " << de << endl;

    if (print_dde)
	de_prev = de;

    if (save_story) {

	// Use high precision here because results may be used for
	// detailed output elsewhere.

	putrq(b->get_dyn_story(), "energy_time", b->get_system_time(),
	      HIGH_PRECISION);
	putrq(b->get_dyn_story(), "potential_energy", epot,
	      HIGH_PRECISION);
	putrq(b->get_dyn_story(), "kinetic_energy", ekin,
	      HIGH_PRECISION);
	putrq(b->get_dyn_story(), "total_energy", ekin+epot,
	      HIGH_PRECISION);
    }
}
