
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  dyn_external.C: functions related to external influences on the system.
//.............................................................................
//    version 1:  Jul 2001, Steve McMillan
//    version 2:  Sep 2001, Steve McMillan
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
//
// Notes from Steve (10/01):
//
// For now, handle here the bits and pieces related to dynamical friction.
// The expression given by Binney & Tremaine is:
//
//	Afric = -4 pi log(Lambda) beta Mfric Vcm rho
//				[erf(X) - 2Xexp(-X^2)/sqrt(pi)] / |Vcm|^3
//
// where the non-obvious terms are defined below and Lambda ~ N(<r).
// We obtain N(<r) from M(<r) assuming a mean mass of 1 Msun and using
// known scalings.  The quantity "1 Msun" is defined if physical units
// are enabled; if they aren't, then it is not clear what scaling we
// should use, or indeed what the meaning of dyamical friction is...
// Beta is a "fudge factor," = 1 according to Binney and Tremaine.
// We introduce some flexibility by allowing beta to be specified on
// the command line, and letting log Lambda default to 1 in the case
// where no physical scale is known.
//
//-------------------------------------------------------------------------

static real beta = 0;				// tunable parameter;
void set_friction_beta(real b) {beta = b;}	// BT say beta = 1

static real Mfric = 0;				// cluster effective mass
void set_friction_mass(real m) {Mfric = m;}

static vector Vcm = 0;				// cluster CM velocity
void set_friction_vel(vector v) {Vcm = v;}

local real density(dyn *b, real r)		// background density
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    return A*x*pow(r*r+a2, 0.5*(x-1))/(4*M_PI*r*r);
}

local real mass(dyn *b, real r)			// mass interior to r
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    return A*pow(r*r+a2, 0.5*x);
}

#define LAMBDA_FAC	1

local real logLambda(dyn *b, real r)
{
    real mass_unit = -1;
    if (b->get_starbase())
	    mass_unit = b->get_starbase()->conv_m_dyn_to_star(1);

    // Only case where this is meaningful is the power-law field.
    real LogLambda;
    if (beta <= 0 || !b->get_pl())
	LogLambda = 0;

    if (mass <= 0)				// no physical mass scale
	LogLambda = 1;
    else {

	// Use M(<r), assuming <m> = 1 Msun.
	LogLambda = 6.6; // Spinnato et al 2003

//	return log(LAMBDA_FAC*mass(b, r)*mass_unit);
    }

    return LogLambda;
}

local real potential(dyn *b, real r)		// background potential;
						// assume x != 1 for now...
{
    real A = b->get_pl_coeff();
    if (A == 0) return 0;

    real a2 = b->get_pl_scale_sq();
    real x = b->get_pl_exponent();

    //#ifdef OLD_ERROR_VERSION_5OCT2001
    //    return -A*pow(r*r+a2, 0.5*x-1);		// *** wrong! ***
    //#else
    return -A*pow(r*r+a2, 0.5*(x-1))/(x-1);
    //#endif
}

static vector Afric = 0;			// frictional acceleration
void set_friction_acc(dyn *b, real r)
{
    if (beta > 0) {

	real sigma2 = sqrt(-2*potential(b, r)/3);  // this is sqrt(2) * sigma;
						   // assume virial equilibrium
	real V = abs(Vcm);
	real X = sqrt(2-b->get_pl_exponent());
	// real X = V/sigma2;			// scaled velocity; BT p. 425

	//#ifndef NEW_4OCT2001
	//	real coeff = beta;		// discard after current runs
	//#else
	real coeff = 4*M_PI*beta*logLambda(b, r);
	//#endif

	if (X > 0.1)

	    Afric = -coeff * Mfric * Vcm * density(b, r) * pow(V, -3)
			   * (erf(X) - 2*X*exp(-X*X)/sqrt(M_PI));
	else

	    // Expand for small X:

	    Afric = -coeff * Mfric * Vcm * density(b, r)
			   * 4 / (3*sqrt(M_PI)*pow(sigma2, 3));

#if 1
	cerr << endl << "set_friction_acc: "; PRL(Afric);
	PRC(beta); PRC(coeff); PRC(Mfric); PRL(density(b, r));
	PRL(Vcm);
	PRC(r); PRC(sigma2); PRC(V); PRL(X);
	PRL(erf(X) - 2*X*exp(-X*X)/sqrt(M_PI));
#endif

    }
}

//-------------------------------------------------------------------------

bool acx_set = false;		// kludge...
static real acx = 0, acx1 = 0;

local inline void set_acx(real A, real c2, real x)
{
    if (c2 > 0) {
	acx = A*pow(c2, x/2);
	if (x != 1)
	    acx1 = acx*x/(sqrt(c2)*(x-1));
	else
	    acx1 = A*(1+0.5*log(c2));
    }
    acx_set = true;
}

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

    real c2 = b->get_pl_cutoff_sq();
    real M = b->get_pl_mass();
    real eps2 = b->get_pl_softening_sq();

#if 0
    PRC(A); PRC(a2); PRL(x);
    PRC(c2); PRC(M); PRL(eps2);
#endif

    if (!acx_set) set_acx(A, c2, x);

    vector dx = pos - b->get_pl_center();
    real dx2 = square(dx);
    real r2 = dx2 + a2;

    real dx1i;
    if (M > 0 || c2 > 0) dx1i = 1/sqrt(dx2+eps2);

    if (x == 1) {				// special case (acx = Ac)

	real dpot;
	if (dx2 <= c2)
	    dpot = -M*dx1i + acx1;
	else {
	    dpot = 0.5*A*log(r2);
	    if (M > 0 || c2 > 0) dpot -= (M-acx)*dx1i;
	}
	pot += dpot;
    }

    if (x != 1 || !pot_only) {			// ugly logic...

	real r1 = pow(r2, 0.5*(1-x));		// in analogy to Plummer case

	if (!pot_only) {

	    real dx3i;
	    if (M > 0 || c2 > 0) dx3i = dx1i/(dx2+eps2);

	    vector vr = dx*(dx*vel);

	    if (dx2 > c2) {

		real r3i = A/(r1*r2);
		acc -= dx*r3i;
		jerk += ((3-x)*vr/r2 - vel)*r3i;

		if (M > 0 || c2 > 0) {
		    dx3i *= M-acx;
		    acc -= dx*dx3i;
		    jerk += (3*vr/(dx2+eps2) - vel)*dx3i;
		}

	    } else if (M > 0) {

		dx3i *= M;
		acc -= dx*dx3i;
		jerk += (3*vr/(dx2+eps2) - vel)*dx3i;

	    }
	}

	// if (b->name_is("2481")) {
	// 	PRC(dx2); PRC(dx*vel); PRL(dx2*abs(acc));
	// }

	if (x != 1) {
	    real dpot;
	    if (dx2 <= c2)
		dpot = -M*dx1i + acx1;
	    else {
		real r1 = pow(r2, 0.5*(1-x));
		dpot = -A/(r1*(1-x));
		if (M > 0 || c2 > 0) dpot -= (M-acx)*dx1i;
	    }
	    pot += dpot;
	}
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

    real c2 = b->get_pl_cutoff_sq();
    real M = b->get_pl_mass();
    real eps2 = b->get_pl_softening_sq();

    if (!acx_set) set_acx(A, c2, x);	//////  *** not implemented yet ***

    real dpot = 0;
    for_all_daughters(dyn, b, bb) {

	vector dx = bb->get_pos() - bb->get_pl_center();
	real dx2 = square(dx);
	real r2 = dx2 + a2;

	real dx1i;
	if (M > 0 || c2 > 0) dx1i = 1/sqrt(dx2+eps2);

	real ddpot;
	if (x == 1) {				// same code as in add_power_law
	    if (dx2 <= c2)
		ddpot = -M*dx1i + acx1;
	    else {
		ddpot = 0.5*A*log(r2);
		if (M > 0 || c2 > 0) ddpot -= (M-acx)*dx1i;
	    }
	} else {
	    if (dx2 <= c2)
		ddpot = -M*dx1i + acx1;
	    else {
		real r1 = pow(r2, 0.5*(1-x));
		ddpot = -A/(r1*(1-x));
		if (M > 0 || c2 > 0) ddpot -= (M-acx)*dx1i;
	    }
	}

	dpot += bb->get_mass()*ddpot;
    }

    return dpot;
}

local inline real power_law_virial(dyn * b)
{
    // Determine the power-law-field component of the virial sum
    // of root node b.

    // *** New embedded mass/cutoff not yet implemented (Steve, 12/01). ***

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

    // We don't actually need the acc_com term, as it should sum
    // to zero in the loop below...

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
    // b is used only as a convenient means of passing and static
    // dyn data.

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

	// Add dynamical friction term to non-escapers only:

	if (getiq(b->get_dyn_story(), "esc") == 0)		// too slow?
	    acc += Afric;
    }
}

real vcirc(dyn *b, vector r)
{
    // Return the circular velocity at position r.  Node b is used only
    // as a convenient means of passing static dyn class data.

    vector acc, jerk;
    real pot;

    get_external_acc(b, r, vector(0), pot, acc, jerk);

    real vc2 = -r*acc;

    if (vc2 > 0)
	return sqrt(vc2);
    else
	return -sqrt(-vc2);	    // vcirc < 0 ==> no circular orbit exists
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
