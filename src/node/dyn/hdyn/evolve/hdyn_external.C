
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//  hdyn_external.C: functions related to external influences on the system.
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

// NOTES:  1. This file is the *only* place where tidal and other
//            external fields are specified.
//
//	   2. Routines add_external() and de_external_pot() MUST be
//	      kept consistent when changes are made!
//
//	   3. We need functions add_xxx() and de_xxx_pot() for each
//	      external field xxx introduced.

local inline void add_tidal(hdyn * b,
			    bool pot_only)
{
    // Add a tidal component to the aceleration of body b.
    // CANNOT neglect the contribution to the jerk.
    // Also include the tidal potential in pot.

    // ASSUME that b is a top-level node.

    real a1 = b->get_alpha1();
    if (a1 == 0) return;

    real a3 = b->get_alpha3();

    vector dx = b->get_pred_pos() - b->get_tidal_center();

    if (!pot_only) {

	vector da_tidal_centrifugal = -vector(a1*dx[0], 0.0, a3*dx[2]);

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

    real x = dx[0];
    real z = dx[2];

    b->inc_pot(0.5*(a1*x*x + a3*z*z));
}

local inline real de_tidal_pot(hdyn * b)
{
    // Determine the tidal component of the potential of body b.
    // ASSUME that tidal_type is set and that b is a top-level node.

    // Add tidal and centrifugal terms for top-level nodes only.
    // (No potential term for the Coriolis force, note.)

    real a1 = b->get_alpha1();
    if (a1 == 0) return 0;

    real a3 = b->get_alpha3();

    real dpot = 0;
    for_all_daughters(hdyn, b, bb) {
	real x = bb->get_pos()[0] - b->get_tidal_center()[0];
	real z = bb->get_pos()[2] - b->get_tidal_center()[0];
	real dp = a1*x*x + a3*z*z;
	dpot += bb->get_mass() * dp;
    }

    return 0.5*dpot;
}

local inline void add_plummer(hdyn * b,
			      bool pot_only)
{
    // Add a Plummer field to the aceleration of body b.
    // CANNOT neglect the contribution to the jerk.
    // Also include the tidal potential in pot.

    // ASSUME that b is a top-level node.

    real M = b->get_p_mass();
    if (M == 0) return;

    real a2 = b->get_p_scale_sq();

    vector dx = b->get_pred_pos() - b->get_p_center();
    real r2 = square(dx) + a2;
    real r1 = sqrt(r2);

    if (!pot_only) {
	vector dv = b->get_pred_vel();
	b->inc_acc(-M*dx/(r1*r2));
	b->inc_jerk(-M*(3*dx*(dx*dv)/r2 - dv));
    }

    b->inc_pot(-M/r1);
}

local inline real de_plummer_pot(hdyn * b)
{
    // Determine the Plummer component of the potential of body b.
    // ASSUME that b is a top-level node.

    real M = b->get_p_mass();
    if (M == 0) return 0;

    real a2 = b->get_p_scale_sq();

    real dpot = 0;
    for_all_daughters(hdyn, b, bb) {
	vector dx = bb->get_pred_pos() - bb->get_p_center();
	real r2 = square(dx) + a2;
	dpot += bb->get_mass() / sqrt(r2);
    }

    return -M*dpot;
}

void add_external(hdyn * b,
		  bool pot_only)	// default = false
{
    // Add external components to the aceleration of body b.
    // Note: CANNOT neglect the contribution to the jerk.
    // Also include the external potential in pot.

    // ASSUME that b is a top-level node.

    int ext = b->get_external_field();

    // Loop through the known external fields.

    if (GETBIT(ext, 0)) add_tidal(b, pot_only);
    if (GETBIT(ext, 1)) add_plummer(b, pot_only);
    // if (GETBIT(ext, 2)) add_other(b, pot_only);
}

real de_external_pot(hdyn * b)
{
    // Determine the external component to the potential of body b.
    // ASSUME that b is a top-level node.

    // Add tidal and centrifugal terms for top-level nodes only.

    real dpot = 0;
    int ext = b->get_external_field();

    // Loop through the known external fields.

    if (GETBIT(ext, 0)) dpot += de_tidal_pot(b);
    if (GETBIT(ext, 1)) dpot += de_plummer_pot(b);
    // if (GETBIT(ext, 2)) dpot += de_other_pot(b);

    return dpot;
}
