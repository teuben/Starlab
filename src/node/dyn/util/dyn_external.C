
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//  dyn_external.C: functions related to external influences on the system.
//.............................................................................
//    version 1:  Jul 2001, Steve McMillan
//    version 2:
//.............................................................................
//
// Global functions:
//
//	void get_external_acc		(single node)
//	real get_external_pot		(root node)
//	real get_tidal_pot		(root node)
//	real get_plummer_pot		(root node)
//	real get_external_virial	(root node)
//	void print_external
//
//.........................................................................

#include "dyn.h"

// NOTES:  1. This file should be the *only* place where tidal and other
//            external fields are specified.
//
//	   2. Must add consistent functions add_xxx(), de_xxx_pot(), and
//	      xxx_virial() for each new external field xxx introduced.
//
//	   3. pot_external uses current positions (get_pos)
//	      get_external uses positions and velocities passed as
//	      arguments (may be pos or pred_pos, depending on the
//	      application)
//

local inline void add_tidal(dyn * b,
			    vector pos,
			    vector vel,
			    real& pot,
			    vector& acc,
			    vector& jerk,
			    bool pot_only)
{
    // Compute the tidal components of the acceleration, jerk,
    // and pot of top-level node b.

    real a1 = b->get_alpha1();
    if (a1 == 0) return;

    real a3 = b->get_alpha3();

    vector dx = pos - b->get_tidal_center();

    if (!pot_only) {

	vector da_tidal_centrifugal = -vector(a1*dx[0], 0.0, a3*dx[2]);
	vector da_coriolis = 2 * b->get_omega()
	    		       * vector(vel[1], -vel[0], 0.0);

	// Must update acc BEFORE computing dj for velocity-dependent forces!

	acc += da_tidal_centrifugal + da_coriolis;

	vector dj_tidal_centrifugal = -vector(a1*vel[0], 0.0, a3*vel[2]);
	vector dj_coriolis = 2 * b->get_omega()
      			       * vector(acc[1], -acc[0], 0.0);

	jerk += dj_tidal_centrifugal + dj_coriolis;
    }

    real x = dx[0];
    real z = dx[2];

    pot += 0.5*(a1*x*x + a3*z*z);
}

local inline real tidal_pot(dyn * b)
{
    // Determine the tidal component of the potential of root node b.

    // Add tidal and centrifugal terms for top-level nodes only.
    // (No potential term for the Coriolis force, note.)

    real a1 = b->get_alpha1();
    if (a1 == 0) return 0;

    real a3 = b->get_alpha3();

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	real x = bb->get_pos()[0] - b->get_tidal_center()[0];
	real z = bb->get_pos()[2] - b->get_tidal_center()[0];
	real dp = a1*x*x + a3*z*z;
	dpot += bb->get_mass() * dp;
    }

    return 0.5*dpot;
}

local inline void add_plummer(dyn * b,
			      vector pos,
			      vector vel,
			      real& pot,
			      vector& acc,
			      vector& jerk,
			      bool pot_only)
{
    // Compute the Plummer-field components of the acceleration, jerk,
    // and pot of top-level node b.

    real M = b->get_p_mass();
    if (M == 0) return;

    real a2 = b->get_p_scale_sq();

    vector dx = pos - b->get_p_center();
    real r2 = square(dx) + a2;
    real r1 = sqrt(r2);

    if (!pot_only) {
	acc -= M*dx/(r1*r2);
	jerk -= M*(3*dx*(dx*vel)/r2 - vel);
    }

    pot -= M/r1;
}

local inline real plummer_pot(dyn * b)
{
    // Determine the Plummer-field component of the potential of
    // root node b.

    real M = b->get_p_mass();
    if (M == 0) return 0;

    real a2 = b->get_p_scale_sq();

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	vector dx = bb->get_pos() - bb->get_p_center();
	real r2 = square(dx) + a2;
	dpot += bb->get_mass() / sqrt(r2);
    }

    return -M*dpot;
}

local inline real plummer_virial(dyn * b)
{
    // Determine the Plummer-field component of the virial sum
    // of root node b.

    real M = b->get_p_mass();
    if (M == 0) return 0;

    real a2 = b->get_p_scale_sq();

    real vir = 0;
    for_all_daughters(dyn, b, bb) {
	vector dx = bb->get_pos() - bb->get_p_center();
	real r2 = square(dx);
	vir += bb->get_mass()*r2*pow(r2+a2, -1.5);
    }

    return -M*vir;
}

void get_external_acc(dyn * b,
		      vector pos,
		      vector vel,
		      real& pot,
		      vector& acc,
		      vector& jerk,
		      bool pot_only)	// default = false
{
    // Compute the external components of the acceleration, jerk,
    // and pot of top-level node b, using the pos and vel provided.
    // b is used mainly as a convenient means of passing flags.

    pot = 0;
    if (!pot_only) acc = jerk = 0;

    unsigned int ext = b->get_external_field();

    if (ext) {

	// Loop through the known external fields.  Must do the
	// velocity-dependent tidal field last (and note that we will
	// have to be more careful when other velocity-dependent fields
	// are added, as *all* accs should be computed before jerks are
	// updated).

	if (GETBIT(ext, 1))
	    add_plummer(b, pos, vel, pot, acc, jerk, pot_only);

	// if (GETBIT(ext, 2))
	//     add_other(b, pos, vel, pot, acc, jerk, pot_only);

	if (GETBIT(ext, 0))
	    add_tidal(b, pos, vel, pot, acc, jerk, pot_only);
    }
}

// Accessors: 

real get_tidal_pot(dyn *b) {return tidal_pot(b);}
real get_plummer_pot(dyn *b) {return plummer_pot(b);}

real get_external_pot(dyn * b,
		      void (*pot_func)(dyn *, real))	// default = NULL
{
    // Determine the external component of the potential of root
    // node b, using current positions.

    real pot = 0;
    unsigned int ext = b->get_external_field();

    if (ext) {

	real dpot = 0;

	// Loop through the known external fields.

	if (GETBIT(ext, 0)) dpot += tidal_pot(b);
	if (GETBIT(ext, 1)) dpot += plummer_pot(b);
	// if (GETBIT(ext, 2)) dpot += other_pot(b);

	if (pot_func) (*pot_func)(b, dpot);	// used by kira to set_pot()

	pot += dpot;
    }

    return pot;
}

real get_external_virial(dyn * b)
{
    // Determine the external component of the potential of root
    // node b, using current positions.

    real vir = 0;
    unsigned int ext = b->get_external_field();

    if (ext) {

	// Loop through the known external, non-tidal fields.

	if (GETBIT(ext, 1)) vir += plummer_virial(b);
	// if (GETBIT(ext, 2)) vir += other_virial(b);
    }

    return vir;
}

static bool sep = false;

local void print_ext(int n)
{
    if (sep) cerr << ",";

    switch(n) {
	case 0:		cerr << "TIDAL";
			break;
	case 1:		cerr << "PLUMMER";
			break;
	default:	cerr << "?";
			break;
    }

    sep = true;
}

void print_external(unsigned int ext,
		    bool shortbits)		// default = false
{
    if (!ext) return;

    if (shortbits)

	// Just print out ext as a binary number.

	printbits(ext);

    else {

	// Want to interpret the bits in ext.  Code follows printbits.C.

	int n = 31;
	while (n >= 0) {if (GETBIT(ext, n)) break; n--;}
	while (n >= 0) {if (GETBIT(ext, n)) print_ext(n); n--;}
	sep = false;
    }
}
