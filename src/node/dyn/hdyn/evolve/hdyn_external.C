
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//  hdyn_external.C: functions related to external incluencef on the system.
//.............................................................................
//    version 1:  Jul 2001, Steve McMillan
//    version 2:
//.............................................................................
//
// Global functions:
//
//	void add_external
//	real de_external_pot
//
//.........................................................................

#include "hdyn.h"

// NOTE: The following two functions are the *only* places where the
//	 tidal field (the only external field known so far) is specified.
//
//	 Routines add_external() and de_external_pot() MUST be kept
//	 consistent when changes are made!

void add_external(hdyn * b,
		  bool pot_only)	// default = false
{
    // Add a tidal component to the aceleration of body b.
    // CANNOT neglect the contribution to the jerk.
    // Also include the tidal potential in pot.

    // ASSUME that b is a top-level node.

    real a1 = b->get_alpha1();
    if (a1 == 0) return;

    real a3 = b->get_alpha3();

    if (!pot_only) {

	vector da_tidal_centrifugal = -vector(a1*b->get_pred_pos()[0],
					      0.0,
					      a3*b->get_pred_pos()[2]);

	vector dj_tidal_centrifugal = -vector(a1*b->get_pred_vel()[0],
					      0.0,
					      a3*b->get_pred_vel()[2]);

	vector da_coriolis = 2 * b->get_omega()
      			       * vector(b->get_pred_vel()[1],
					-b->get_pred_vel()[0],
					0.0);

	// Must update acc BEFORE computing dj!

	b->inc_acc(da_tidal_centrifugal + da_coriolis);

	vector dj_coriolis = 2 * b->get_omega()
      			       * vector(b->get_acc()[1],
					-b->get_acc()[0],
					0.0);

	b->inc_jerk(dj_tidal_centrifugal + dj_coriolis);

	// PRL(da_tidal_centrifugal);
	// PRL(da_coriolis);
	// PRL(abs(da_tidal_centrifugal)/abs(b->get_acc()));

    }

    real x = b->get_pred_pos()[0];
    real z = b->get_pred_pos()[2];

    b->inc_pot(0.5*(a1*x*x + a3*z*z));
}

real de_external_pot(hdyn * b)
{
    // Determine tidal component to the potential of body b.
    // ASSUME that tidal_type is set and that b is a top-level node.

    // Add tidal and centrifugal terms for top-level nodes only.
    // (No potential term for the Coriolis force, note.)

    real a1 = b->get_alpha1();
    if (a1 == 0) return 0;

    real a3 = b->get_alpha3();

    real dpot = 0;
    for_all_daughters(hdyn, b, bb) {
	real x = bb->get_pos()[0];
	real z = bb->get_pos()[2];
	real dp = a1*x*x + a3*z*z;
	dpot += bb->get_mass() * dp;
    }

    dpot *= 0.5;
    return dpot;
}
