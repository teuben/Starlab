
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// refine_mass2.C:  More careful determination of cluster mass and center,
//		    for use with an external non-tidal field.
//
// Externally visible functions:
//
//	void refine_cluster_mass2

#include "dyn.h"

local inline vector ext_acc(dyn *b, vector pos)
{
    vector v = 0, acc, j;
    real p;
    get_external_acc(b, pos, v, p, acc, j);		// (discard p and j)
    return acc;
}

local inline real ext_pot(dyn *b, vector pos)
{
    vector v = 0, a, j;
    real pot;
    get_external_acc(b, pos, v, pot, a, j, true);	// (discard a and j)
    return pot;
}

// Handy, for now...

static vector ext_center, Rhat, acc_CM;
static real R;

local real g1(dyn *b, real r)
{
    vector pos = ext_center + r*Rhat;
    return (ext_acc(b, pos) - acc_CM) * Rhat + pow(r-R, -2);
}

local real g2(dyn *b, real r)
{
    vector pos = ext_center + r*Rhat;
    return (ext_acc(b, pos) - acc_CM) * Rhat - pow(r-R, -2);
}

#define EPS_R 1.e-4

local inline real solve(real (*g)(dyn*, real), dyn *b, real r1, real r2)
{
    // Return the root of g(r) = 0 (if any) between r1 and r2.
    // Avoid evaluation of g at r1 or r2.  Start the search at r2,
    // proceed toward r1 (could have r1 > r2), and return the first
    // zero found.

    int fac = (r1 < r2);
    fac = 2*fac - 1;		// = +1 if r1 < r2, -1 otherwise
    r1 *= 1 + fac*EPS_R;
    r2 *= 1 - fac*EPS_R;

    real dr = 0.01 * (r2 - r1);

    // Search for a zero, starting at r2.

    real g2 = g(b, r2);
    real r = r2 - dr;
    while (fac*(r - r1) >= 0 && g(b, r)*g2 > 0) r -= dr;

    if (fac*(r - r1) < 0) return r;

    // Refine by bisection.

    r1 = r;
    r2 = r + dr;
    real g1 = g(b, r1);
    g2 = g(b, r2);

    while (abs(r2/r1 - 1) < EPS_R) {
	r = 0.5*(r1+r2);
	real gr = g(b, r);
	if (gr == 0) return r;
	if (gr*g1 > 0) {
	    r1 = r;
	    g1 = gr;
	} else {
	    r2 = r;
	    g2 = gr;
	}
    }

    // Final refinement by linear interpolation.

    r = r1 + (r2-r1)*(0-g1)/(g2-g1);

    return r;
}

local void get_rL(dyn *b,		// root node
		  real M,		// current mass estimate
		  vector center,	// current center estimate
		  real& r_L1,
		  real& r_L2)
{
    // Find the inner Lagrangian point in the external field.
    // Construction of the external field functions is such that
    // it is just as easy to work directly with 3-D vectors...

    ext_center = b->get_external_center();
    R = abs(center - ext_center);
    Rhat = (center - ext_center)/(R*M);	// factor of M avoids passing M to g
    acc_CM = ext_acc(b, center);

    r_L1 = solve(g1, b, sqrt(b->get_external_scale_sq()), R);
    r_L2 = solve(g2, b, R, 10*R);
}

local int bitcount(unsigned int i)
{
    // Count nonzero bits.  From K&R.

    int b;
    for (b = 0; i != 0; i >>= 1)
	if (i & 01) b++;
    return b;
}

#define M_TOL 1.e-4
#define M_ITER_MAX 20

void refine_cluster_mass2(dyn *b)
{
    unsigned int ext = b->get_external_field();

    if (!ext || b->get_tidal_field()) return;

    // Self-consistently determine the total mass within the outermost
    // closed zero-velocity surface under the specified external field(s).
    // Use a point-mass approximation for the cluster potential and iterate
    // until the actual mass within the surface agrees with the mass used
    // to generate the surface.
    //
    // Experimental code, implemented by Steve, 8/01.

    // Method only works for a single external, velocity-independent field...

    if (bitcount(ext) != 1) return;

    // Use the standard center as our starting point.  The center will be
    // redetermined self-consistently, along with the mass.

    vector center, vcenter;
    get_std_center(b, center, vcenter);

    real M_inside = total_mass(b);
    real M = -1;				// (to start the loop)
    int N_inside;
    vector cen = center, vcen = vcenter;

    cerr << endl << "  refine_cluster_mass2: getting mass by iteration"
	 << endl << "  initial total system mass = " << M_inside
	 << endl;
    PRI(2); PRL(center);

    int iter = 0;

    while (iter++ < M_ITER_MAX
	   && M_inside > 0
	   && abs(M_inside/M - 1) > M_TOL) {

	M = M_inside;
	center = cen;
	vcenter = vcen;

	// Reinitialize:

	M_inside = 0;
	N_inside = 0;
	cen = vcen = 0;

	// Determine the Lagrangian points of the (point-mass) cluster
	// in the external field, then count stars within the limiting
	// equipotential.

	real r_L1, r_L2;
	get_rL(b, M, center, r_L1, r_L2);
	cerr << endl;
	PRC(R); PRL(Rhat);
	PRC(r_L1); PRL(r_L2);

	// Limiting potential: 

	real phi_L1 = ext_pot(b, ext_center + r_L1*Rhat) - M/(R-r_L1);
	real phi_L2 = ext_pot(b, ext_center + r_L2*Rhat) - M/(r_L2-R);
	PRC(phi_L1); PRL(phi_L2);

	real phi_lim = max(phi_L1, phi_L2);	// maximize the cluster mass
	real r_L = max(R-r_L1, r_L2-R);
	PRC(phi_lim); PRL(r_L);

	for_all_daughters(dyn, b, bb) {
	    real r = abs(bb->get_pos() - center);
	    if (r < r_L) {
		if (r == 0 
		    || -M/r + ext_pot(b, bb->get_pos()) < phi_lim) {
		    N_inside++;
		    M_inside += bb->get_mass();
		    cen += bb->get_mass()*bb->get_pos();
		    vcen += bb->get_mass()*bb->get_vel();
		}
	    }
	}

	if (M_inside > 0) {
	    cen /= M_inside;
	    vcen /= M_inside;
	}

	PRC(iter); PRC(N_inside); PRL(M_inside);
	PRL(cen);
    }

    if (iter >= M_ITER_MAX)
	warning("refine_cluster_mass2: too many iterations");
}
