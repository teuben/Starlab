
//// radial_density:  Compute the 1-dimensional radial density profile
////                  of an N-body system.  The function fills in the density
////                  array corresponding to the provided radius array.
////                  The tool prints out radius and density in a form
////                  suitable for plotting.
////
////                  If a current density center is found in the root dyn
////                  story it is used as the system center for the density
////                  computation.  Otherwise, if a valid center of mass is
////                  found, it is used.  If neither center is found, the
////                  geometric center is used.
////
//// Options:     -n    specify number of radial bins [100]
////              -r    specify maximum radius [take from data]

//-----------------------------------------------------------------------------
//   version 1:  Apr 2003   Steve McMillan
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

void get_radial_densities(dyn * b, int n_zones,
			  vector cpos,
			  real r[], real rho[])
{
    if (n_zones < 2) return;

    // Determine the mass per bin.

    for (int i = 0; i < n_zones; i++) rho[i] = 0;

    for_all_daughters(dyn, b, bb) {

	real rr = abs(bb->get_pos()-cpos);	// would be more efficient
						// to use the square...

	// Don't assume linear binning, so locate the bin the hard way.
	// Could do much better here if we sorted the radii first.

	for (int i = 0; i < n_zones; i++)
	    if (rr <= r[i]) {
		rho[i] += bb->get_mass();
		break;
	    }
    }

    // Convert from mass to density.

    real v0 = 0;
    for (int i = 0; i < n_zones; i++) {
	real v1 = pow(r[i], 3);
	rho[i] /= (4*M_PI/3) * (v1 - v0);	// dM --> dM/dV
	v0 = v1;
    }
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  get_radial_densities() as a tool.
//-----------------------------------------------------------------------------

#define N_DEFAULT 100

main(int argc, char ** argv)
{
    int n_zones = N_DEFAULT;
    real r_max = 0;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "n:r:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'n': n_zones = atoi(poptarg);
		      break;

	    case 'r': r_max = atof(poptarg);
		      break;

            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            

    if (n_zones <= 1) n_zones = N_DEFAULT;
    PRC(r_max); PRL(n_zones);

    dyn *b;

    while (b = get_dyn()) {

	real r[n_zones], rho[n_zones];

	// Use the density center if known and up to date (preferred).
	// Otherwise, use modified center of mass, if known and up to date.

	vector cpos, cvel;
	get_std_center(b, cpos, cvel);

	cpos -= b->get_pos();			// std_center quantities
	cvel -= b->get_vel();			// include the root node

	// See if we need to compute r_max.

	if (r_max <= 0) {
	    r_max = 0;
	    for_all_daughters(dyn, b, bb) {
		real r2 = square(bb->get_pos()-cpos);
		if (r2 > r_max) r_max = r2;
	    }
	    r_max = sqrt(r_max);
	}

	// Set up the radial array.

	for (int i = 0; i < n_zones; i++)
	    r[i] = (i+1) * r_max / n_zones;	// r[i] = outer edge of zone i

	// Compute and print the density array.

	get_radial_densities(b, n_zones, cpos, r, rho);

	real r0 = 0;
	for (int i = 0; i < n_zones; i++) {
	    real r1 = r[i];
	    cout << i << " " << (r0+r1)/2 << " " << rho[i] << endl;
	    r0 = r1;
	}

	delete b;
    }
}

#endif
