
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
// Member functions:
//
//	real dyn::get_external_scale_sq()
//	vector dyn::get_external_center();
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

//-------------------------------------------------------------------------
// Tidal field (quadrupole):
//-------------------------------------------------------------------------

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
    // Determine the tidal component of the potential energy
    // of root node b.

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

//-------------------------------------------------------------------------
// Plummer field:
//-------------------------------------------------------------------------

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
	real r3i = 1/(r1*r2);
	acc -= M*dx*r3i;
	jerk += M*(3*dx*(dx*vel)/r2 - vel)*r3i;
    }

    pot -= M/r1;
}

local inline real plummer_pot(dyn * b)
{
    // Determine the Plummer-field component of the potential energy
    // of root node b.

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

    // Don't make any assumptions about the locations of the
    // center of mass of the center of the Plummer field...

    vector com_pos, com_vel;
    compute_com(b, com_pos, com_vel);
    // PRL(com_pos);

    vector dR = com_pos - b->get_p_center();
    vector acc_com = dR * pow(square(dR)+a2, -1.5);
    // PRL(acc_com);

    // Don't actually need the acc_com term, as it should sum to zero
    // in the loop below...

    real vir = 0;
    for_all_daughters(dyn, b, bb) {
	vector dr = bb->get_pos() - com_pos;
	dR = bb->get_pos() - b->get_p_center();
	vector acc_ext = dR * pow(square(dR)+a2, -1.5);
	real dvir = bb->get_mass()*dr*(acc_ext - acc_com);
	// PRL(dvir);
	vir += dvir;
    }
    // PRL(vir);

    return -M*vir;
}

//-------------------------------------------------------------------------
// Power-law field, M(r) = A r^x (and note that exponent = 0 reduces to
// a Plummer field with mass = A):
//-------------------------------------------------------------------------




//-------------------------------------------------------------------------
// For now, handle here the bits and pieces related to dynamical friction.
// *** To be cleaned up!!

static real beta = 0;
void set_friction_beta(real b) {beta = b;}

static real Mfric = 0;
void set_friction_mass(real m) {Mfric = m;}

static vector Vfric = 0;
void set_friction_vel(vector v) {Vfric = v;}

local real density(dyn *b, real r)
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    return A*x*pow(r*r+a2, 0.5*(x-1))/(4*M_PI*r*r);
}

local real potential(dyn *b, real r)		// Assume x != 1 for now...
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    return -A*pow(r*r+a2, 0.5*x-1);
}

static vector Afric = 0;
void set_friction_acc(dyn *b, real r)
{
    Afric = -beta*Mfric*Vfric*density(b, r)/pow(-potential(b, r), 1.5);
}

//-------------------------------------------------------------------------




local inline void add_power_law(dyn * b,
				vector pos,
				vector vel,
				real& pot,
				vector& acc,
				vector& jerk,
				bool pot_only)
{
    // Compute the power-law-field components of the acceleration, jerk,
    // and pot of top-level node b.

    real A = b->get_pl_coeff();
    if (A == 0) return;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    vector dx = pos - b->get_pl_center();
    real r2 = square(dx) + a2;

    if (x == 1)					// special case
	pot += 0.5*A*log(r2);

    if (x != 1 || !pot_only) {			// ugly logic for efficiency

	real r1 = pow(r2, 0.5*(1-x));		// in analogy to Plummer case

	if (!pot_only) {
	    real r3i = 1/(r1*r2);
	    acc -= A*dx*r3i;
	    jerk += A*((3-x)*dx*(dx*vel)/r2 - vel)*r3i;
	}

	if (x != 1)
	    pot -= A/(r1*(1-x));
    }
}

local inline real power_law_pot(dyn * b)
{
    // Determine the power-law-field component of the potential energy
    // of root node b.

    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {
	vector dx = bb->get_pos() - bb->get_pl_center();
	real r2 = square(dx) + a2;
	if (x == 1)
	    dpot += bb->get_mass() * log(r2);
	else
	    dpot += bb->get_mass() * pow(r2, 0.5*(x-1));
    }

    if (x == 1)
	dpot *= 0.5;
    else
	dpot /= x-1;

    return A*dpot;
}

local inline real power_law_virial(dyn * b)
{
    // Determine the power-law-field component of the virial sum
    // of root node b.

    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    // Don't make any assumptions about the locations of the
    // center of mass of the center of the power-law field...

    vector com_pos, com_vel;
    compute_com(b, com_pos, com_vel);
    // PRL(com_pos);

    vector dR = com_pos - b->get_pl_center();
    vector acc_com = dR * pow(square(dR)+a2, 0.5*(x-3));
    // PRL(acc_com);

    // Don't actually need the acc_com term, as it should sum to zero
    // in the loop below...

    real vir = 0;
    for_all_daughters(dyn, b, bb) {
	vector dr = bb->get_pos() - com_pos;
	dR = bb->get_pos() - b->get_pl_center();
	vector acc_ext = dR * pow(square(dR)+a2, 0.5*(x-3));
	real dvir = bb->get_mass()*dr*(acc_ext - acc_com);
	// PRL(dvir);
	vir += dvir;
    }
    // PRL(vir);

    return -A*vir;
}

//-------------------------------------------------------------------------
// General "external" functions:
//-------------------------------------------------------------------------

// Member functions:

real dyn::get_external_scale_sq()
{
    // Just enumerate the possiblilties...

    if (!get_external_field())
	return 0;
    else if (get_tidal_field())
	return 0;
    else if (get_plummer())
	return p_scale_sq;
    else if (get_pl())
	return pl_scale_sq;
    else
	return 0;
}

vector dyn::get_external_center()
{
    // Just enumerate the possiblilties...

    if (!get_external_field())
	return 0;
    else if (get_tidal_field())
	return tidal_center;
    else if (get_plummer())
	return p_center;
    else if (get_pl())
	return pl_center;
    else
	return 0;
}

// Other functions:

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

	if (GETBIT(ext, 2))
	    add_power_law(b, pos, vel, pot, acc, jerk, pot_only);

	// if (GETBIT(ext, 3))
	//     add_other(b, pos, vel, pot, acc, jerk, pot_only);

	if (GETBIT(ext, 0))
	    add_tidal(b, pos, vel, pot, acc, jerk, pot_only);

	// Add dynamical friction term:

	acc += Afric;
    }
}

// Accessors: 

real get_tidal_pot(dyn *b) {return tidal_pot(b);}
real get_plummer_pot(dyn *b) {return plummer_pot(b);}
real get_power_law_pot(dyn *b) {return power_law_pot(b);}

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
	if (GETBIT(ext, 2)) dpot += power_law_pot(b);
	// if (GETBIT(ext, 3)) dpot += other_pot(b);

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
	if (GETBIT(ext, 2)) vir += power_law_virial(b);
	// if (GETBIT(ext, 3)) vir += other_virial(b);
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
	case 2:		cerr << "POWER-LAW";
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
