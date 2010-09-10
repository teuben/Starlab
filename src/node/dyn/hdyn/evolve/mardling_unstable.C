
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

#include "hdyn.h"

// This is simply a wrapper for the Fortran-77 function (Steve, 8/09).

// *** Still need to develop code for use if we don't have F77    ***
// *** -- probably should force the Aarseth criterion to be used. ***

#define NSTAB_F77 F77_FUNC(nstab, NSTAB)

extern "C" int NSTAB_F77(real*, real*, real*, real*, real*, real*, real*);

int mardling_stable(hdyn *b, hdyn *sister, kepler& outerkep)
{
    // Check the stability of the triple consisting of binary b and
    // its sister, described by kepler structure outerkep.  In
    // application of the stability criterion, sister is always
    // regarded as a point mass.

    hdyn *od = b->get_oldest_daughter();
    if (!od) return 1;
    hdyn *yd = od->get_younger_sister();
    kepler innerkep;

    if (b->get_kepler() == NULL) {
	innerkep.set_time(b->get_time());
	innerkep.set_total_mass(b->get_mass());
	innerkep.set_rel_pos(od->get_pos() - yd->get_pos());
	innerkep.set_rel_vel(od->get_vel() - yd->get_vel());
	innerkep.initialize_from_pos_and_vel();
    } else
	innerkep = *(b->get_kepler());

    real sigma = outerkep.get_period()/innerkep.get_period();
    real ei0 = innerkep.get_eccentricity();
    real eo = outerkep.get_eccentricity();
    real relinc = acos(innerkep.get_normal_unit_vector()
		        * outerkep.get_normal_unit_vector());
    real m1 = od->get_mass();
    real m2 = yd->get_mass();
    real m3 = sister->get_mass();

    return 1 - NSTAB_F77(&sigma, &ei0, &eo, &relinc, &m1, &m2, &m3);
    // sigma      = period ratio (outer/inner) ** should be > 1 **
    // ei0        = initial inner eccentricity
    // eo         = outer eccentricity
    // relinc     = relative inclination (radians)
    // m1, m2, m3 = masses (any units; m3=outer body)

    // Note: A return value of 1 from NSTAB means that we can't treat
    // this object as stable.
}



