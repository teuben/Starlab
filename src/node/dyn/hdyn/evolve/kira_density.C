
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// Functions associated with calculating densities associated with
// particles in the system.
//
// Externally visible function:
//
//	bool kira_calculate_densities
//
// Major reorganization to remove compile-time GRAPE selection.
//					      Steve 6/04

#include "hdyn.h"

// The following function is a switch between the various energy-
// computation means available -- currently, GRAPE or nothing!
//
// If other schemes are added (e.g. treecode), add them HERE.

bool kira_calculate_densities(hdyn* b, vec& cod_pos, vec& cod_vel)
{
    // Density computation (currently limited to GRAPE systems).
    // Called only by log_output() in kira_log.C.

    bool status = false;

    if (b->has_grape()) {

	cerr << "computing densities using GRAPE..." << endl;
	real cpu0 = cpu_time();

	// The second argument determines the squared radius at which
	// particles are deemed to have zero densities.  This allows
	// discrimination against low-density particles, and also limits
	// costly repeat GRAPE calls.

	if (b->has_grape6())
	    status = grape6_calculate_densities(b, 0.1);  // (densities are
							  //  saved in particle
							  //  dyn stories)
	else
	    status = grape4_calculate_densities(b, 0.1);

	real cpu1 = cpu_time();
	compute_mean_cod(b, cod_pos, cod_vel);
	real cpu2 = cpu_time();

	cerr << "CPU times:  density " << cpu1 - cpu0
	     << "  cod " << cpu2 - cpu1
	     << endl;

//	for_all_daughters(hdyn, b, bb) {
//	    real dens = getrq(bb->get_dyn_story(), "density");
//	    PRC(bb->format_label()); PRL(dens);
//	}

    } else {

	// Skip (too expensive if no GRAPE is available...).
	// Need a standalone (kd)tree module...

	cerr << "skipping density calculation..." << endl;
    }

    return status;
}
