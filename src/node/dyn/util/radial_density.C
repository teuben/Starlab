
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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

// Prototype function.  Could also be implemented using get_density_profile()
// with a mass selection function that accepts all stars.


// Sorting code is taken almost verbatim from lagrad.C.

typedef struct {
    real r_sq;
    real mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rm_pair_ptr) pi)->r_sq > ((rm_pair_ptr) pj)->r_sq)
        return +1;
    else if (((rm_pair_ptr)pi)->r_sq < ((rm_pair_ptr)pj)->r_sq)
        return -1;
    else
        return 0;
}

int get_radial_densities(dyn *b, vec cpos,
			 int n_zones, real r[], real rho[])
{
    if (n_zones < 2) return 1;

    // Set up an array of (r_sq, mass) pairs.

    int n = b->n_daughters();	// (NB implicit loop through the entire system)

    // Would be possible to determine n and set up the array simultaneously
    // using malloc and realloc.  Not so easy with new...
    
    rm_pair_ptr table = new rm_pair[n];

    if (!table) {
	cerr << "get_radial_densities: insufficient memory for table"
	     << endl;
	return 1;
    }

    int i = 0;
    for_all_daughters(dyn, b, bi) {
	table[i].r_sq = square(bi->get_pos() - cpos);
	table[i].mass = bi->get_mass();
	i++;
    }

    // Sort the array by radius (may repeat work done elsewhere...).

    qsort((void *)table, (size_t)i, sizeof(rm_pair), compare_radii);

    // Initialize the density array.

    int j;
    for (j = 0; j < n_zones; j++) rho[j] = 0;

    // Bin the (ordered) data.

    j = 0;
    real rj2 = r[j]*r[j];

    for (i = 0; i < n; i++) {
	real ri_sq = table[i].r_sq;
	while (ri_sq > rj2) {
	    j++;
	    if (j >= n_zones) break;
	    rj2 = r[j]*r[j];
	}
	if (j >= n_zones) break;
	rho[j] += table[i].mass;
    }

    // Convert from mass to density.

    real v0 = 0;	// assume that the first zone extends in to r = 0
    for (j = 0; j < n_zones; j++) {
	real v1 = pow(r[j], 3);
	rho[j] /= (4*M_PI/3) * (v1 - v0);	// dM --> dM/dV
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

	vec cpos, cvel;
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

	// Set up a linear radial array.

	for (int i = 0; i < n_zones; i++)
	    r[i] = (i+1) * r_max / n_zones;	// r[i] = outer edge of zone i

	// Compute and print the density array.

	get_radial_densities(b, cpos, n_zones, r, rho);

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
