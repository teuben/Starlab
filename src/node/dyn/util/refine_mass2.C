
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

void refine_cluster_mass2(dyn *b,
			  int verbose)		// default = 0
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

    // Do nothing if all we want is to set the dyn story and the current
    // values are up to date.

    if (verbose == 0
	&& getrq(b->get_dyn_story(), "bound_center_time")
		== b->get_system_time())
	return;

    // Use the standard center as our starting point.  The center will be
    // redetermined self-consistently, along with the mass.

    vector center, vcenter;
    get_std_center(b, center, vcenter);

    // Choose the initial mass to include only the volume between the
    // standard center and the center of the external potential.

    real M_inside = 0;
    R = abs(center - b->get_external_center());	// global R, same definition

    for_all_daughters(dyn, b, bb)
	if (abs(bb->get_pos() - center) <= R)
	    M_inside += bb->get_mass();

    real M = -1;				// (to start the loop)
    int N_inside;
    vector cen = center, vcen = vcenter;

    if (verbose) {
	cerr << endl << "  refine_cluster_mass2: getting mass by iteration"
	     << endl << "  initial total system mass = " << M_inside
	     << endl << "  initial center = " << center
	     << endl;
    }

    int iter = 0;
    real r_L = 0, phi_lim = 0;

    // Iteration can (should!) run away to zero mass if the cluster
    // density is too low relative to the local tidal field strength.
    // Keep track of the "50% center of mass" during the iteration as
    // a convenient measure of the cluster center in this case.

    real M0 = M_inside;
    vector cen50 = center, vcen50 = vcenter;
    bool set50 = false;

    while (iter++ < M_ITER_MAX
	   && M_inside > 0
	   && abs(M_inside/M - 1) > M_TOL) {

	// Set current mass and center:

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

	if (verbose > 1) {
	    cerr << endl;
	    PRC(R); PRC(r_L1); PRC(r_L2);
	}

	// Limiting potential: 

	real phi_L1 = ext_pot(b, ext_center + r_L1*Rhat) - M/(R-r_L1);
	real phi_L2 = ext_pot(b, ext_center + r_L2*Rhat) - M/(r_L2-R);
	// PRC(phi_L1); PRC(phi_L2);

	phi_lim = max(phi_L1, phi_L2);		// maximize the cluster mass
	r_L = max(R-r_L1, r_L2-R);

	if (verbose > 1) PRL(r_L);

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

	// Linearly interpolate an estimate of the 50% center.

	if ((M > 0.5*M0 && M_inside <= 0.5*M0)
	    || (M < 0.5*M0 && M_inside >= 0.5*M0)) {
	    cen50 = center + (0.5*M0-M)*(cen-center)/(M_inside-M);
	    vcen50 = vcenter + (0.5*M0-M)*(vcen-vcenter)/(M_inside-M);
	    set50 = true;
	}

	if (verbose > 1) {
	    PRI(2); PRC(iter); PRC(N_inside); PRL(M_inside);
	    PRI(2); PRL(cen);
	}
    }

    if (iter >= M_ITER_MAX)
	warning("refine_cluster_mass2: too many iterations");

    if (verbose == 1) {
	PRI(2); PRC(iter); PRC(N_inside); PRL(M_inside);
	PRI(2); PRL(cen);
    }

    bool modify_center = false;

    if (iter >= M_ITER_MAX || M_inside < 0.01*M0) {

	// Looks like the cluster no longer exists.  Use 50% or ext center.

	if (set50) {
	    center = cen50;
	    vcenter = vcen50;
	} else {
	    center = ext_center;
	    vcenter = 0;
	}
	modify_center = true;

    } else {

	center = cen;
	vcenter = vcen;

    }

    // Now center and vcenter should be usable, even if M_inside and
    // N_inside aren't very meaningful.

    if (verbose && modify_center) {
	PRI(2); PRL(center);
    }

    // Write our best estimate of the center to the dyn story.

    putrq(b->get_dyn_story(), "bound_center_time", b->get_system_time());
    putvq(b->get_dyn_story(), "bound_center_pos", center);
    putvq(b->get_dyn_story(), "bound_center_vel", vcenter);

    // Repeat the inner loop above and flag stars as escapers or not.

    bool disrupted = (iter >= M_ITER_MAX || M_inside < 0.01*M0);

    for_all_daughters(dyn, b, bb) {
	bool escaper = true;
	if (!disrupted) {
	    real r = abs(bb->get_pos() - center);
	    if (r < r_L
		&& (r == 0 || -M/r + ext_pot(b, bb->get_pos()) < phi_lim))
		escaper = false;
	}
	putiq(bb->get_dyn_story(), "esc", escaper);
    }
}
