
//// lagrad:  Compute Lagrangian (mass) radii for an N-body system.
////          The radii are stored in the system (root) dyn story.
////
////          If a current density center is found in the root dyn story,
////          the Lagrangian radii are calculated relative to it.
////          Otherwise, if a valid center of mass is found, it is used.
////          If neither center is found, the geometric center is used.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -n    specify number of Lagrangian zones (linear in mass) [4]
////              -s    use "special" nonlinear binning:
////                        0.005, 0.01, 0.02, 0.05,
////                        0.1, 0.25, 0.5, 0.75, 0.9
////              -t    same as -n 10

//-----------------------------------------------------------------------------
//   version 1:  Jul 1998   Simon Portegies Zwart
//.............................................................................
//   non-local functions: 
//      compute_gereral_luminosity_radii
//      compute_luminosity_radii_quartiles
//      compute_luminosity_radii_percentiles
//-----------------------------------------------------------------------------

//#include "dyn.h"
#include "sstar_to_dyn.h"
#include "star.h"

#ifndef TOOLBOX

typedef  struct
{
    real  radius;
    real  mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)
{
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return(1);
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return(-1);
    else
        return(0);
}

//-----------------------------------------------------------------------------
//  compute_general_luminosity_radii  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

static real nonlin_masses[9] = {0.005, 0.01, 0.02, 0.05, 0.1,
				0.25, 0.5, 0.75, 0.9};

real compute_projected_luminosity_radii(dyn * b, int axis, bool Llagr,
					int nzones,
					bool nonlin,
					boolfn bf) {

  real r1_r9 = 0;
    if (nzones < 2) return nzones;
    if (nonlin && nzones != 10) return nzones;		// Special case

    // Note: nzones specifies the number of radial zones to consider.
    //	     However, we only store radii corresponding to the nzones-1 
    //	     interior zone boundaries.

    real* mass_percent = new real[nzones-1];
    if (mass_percent == NULL) {
	cerr << "compute_general_luminosity_radii: "
	     << "not enough memory left for mass_percent\n";
	return r1_r9;
    }

    int n = 0;
    if (bf == NULL)
	n = b->n_leaves();
    else {
	for_all_leaves(dyn, b, bb)
	    if ((*bf)(bb)) n++;
    }

    rm_pair_ptr rm_table = new rm_pair[n];

    if (rm_table == NULL) {
	cerr << "compute_general_luminosity_radii: "
	     << "not enough memory left for rm_table\n";
	return r1_r9;
    }

    // Use the density center if known and up to date.
    // Otherwise, use center of mass if known and up to date.
    // Otherwise, use the geometric center.

    vec dc_pos = 0;
    bool try_com = false;

    if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {

	if (getrq(b->get_dyn_story(), "density_center_time")
		!= b->get_system_time()) {
	    warning("lagrad: neglecting out-of-date density center");
	    try_com = true;
	} else
	    dc_pos = getvq(b->get_dyn_story(), "density_center_pos");

    }

    if (try_com && find_qmatch(b->get_dyn_story(), "com_pos")) {

	if (getrq(b->get_dyn_story(), "com_time")
		!= b->get_system_time()) {
	    warning("lagrad: neglecting out-of-date center of mass");
	} else
	    dc_pos = getvq(b->get_dyn_story(), "com_pos");

    }
    if(axis>=0)
      dc_pos[axis] = 0;

    // Set up an array of (radius, mass) pairs.  Also find the total
    // mass of all nodes under consideration.

    real total_mass = 0;

    int i = 0;
    real lstar = 0;
    vec pos = 0;
    for_all_leaves(dyn, b, bi) {
	if (bf == NULL || (*bf)(bi)) {
	  if(Llagr) {
	    if (find_qmatch(bi->get_starbase()->get_star_story(), "L_eff")) {
	      lstar = getrq(bi->get_starbase()->get_star_story(), "L_eff");
	    }
	    else if (b->get_use_sstar()){
	      lstar = ((star*)bi->get_starbase())->get_luminosity();
	    }
	  }
	  else
	    lstar = bi->get_mass(); 

	    total_mass += lstar;
	    pos = bi->get_pos();
	    if(axis>=0)
	      pos[axis] = 0;
	    rm_table[i].radius = abs(pos - dc_pos);
	    rm_table[i].mass = lstar;
	    i++;
	}
    }

    // Sort the array by radius.

    qsort((void *)rm_table, (size_t)n, sizeof(rm_pair), compare_radii);

    // Determine Lagrangian radii.

    // cerr << "Determining Lagrangian radii 1" << endl << flush;

    int k;
    for (k = 0; k < nzones-1; k++) {
        if (!nonlin) 
	    mass_percent[k] = ((k + 1) / (real)nzones) * total_mass;
	else
	    mass_percent[k] = nonlin_masses[k] * total_mass;
    }

    real *rlagr = new real[nzones-1];
    real cumulative_mass = 0.0;
    i = 0;

    switch(axis) {
    case 0:
      cerr <<"     (projected along x-axis)\n";
      break;
    case 1:
      cerr <<"     (projected along y-axis)\n";
      break;
    case 2:
      cerr <<"     (projected along x-axis)\n";
      break;
    default:
      cerr <<"     (unprojected)\n";
    }
    
    cerr << "    r_lagr =";

    for (k = 0; k < nzones-1; k++) {

        while (cumulative_mass < mass_percent[k])
	    cumulative_mass += rm_table[i++].mass;

	rlagr[k] = rm_table[i-1].radius;
	
	cerr << "  " << rlagr[k];
    }
    cerr << endl;

    
    if (nzones==10 && rlagr[8]>0)
      r1_r9 = rlagr[0]/rlagr[8];

    delete [] mass_percent;
    delete [] rm_table;
    delete [] rlagr;

    return r1_r9;
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  compute_luminosity_radii() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    int n = 0;
    int axis = -1;
    bool  c_flag = false;      // if TRUE: a comment given on command line
    bool  t_flag = false;      // if TRUE: percentiles rather than quartiles
    bool  nonlin = false;
    bool  Llagr  = true;

    check_help();

    void compute_luminosity_radii(dyn *, int, bool, bool);

  extern char *poptarg;
    int c;
    const char *param_string = "a:c:mn:xyzst";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'a': axis = atoi(poptarg);
		      break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'x': axis = 0;
		      break;
	    case 'y': axis = 1;
		      break;
	    case 'z': axis = 2;
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 'm': Llagr = false;
		      break;
	    case 's': n = 10;
	      nonlin = true;  // fall through
	    case 't': t_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    dyn *b;

    while (b = get_dyn()) {

      //        addstar(b);
	
        if (c_flag == TRUE)
            b->log_comment(comment);

        b->log_history(argc, argv);

	if (Llagr) 
	  cerr << "Lagrangian mass radii"<<endl;
	else
	  cerr << "Lagrangian luminosity radii"<<endl;
	
	if (axis>=0) {
	  cerr << "\nprojected along the ";
	  switch(axis) {
	    case 0: cerr << "x";
	            break;
	    case 1: cerr << "y";
	            break;
	    case 2: cerr << "z";
	            break;
	  }
	  cerr << "-axis\n";
	}

	real r1_r9 = 0;
	if (t_flag)
	  real r1_r9 = compute_projected_luminosity_radii(b, axis, Llagr, 10);
	else {
	  if (n <= 1)
	    compute_projected_luminosity_radii(b, axis, Llagr, 4);
	  else
	    compute_projected_luminosity_radii(b, axis, Llagr, n, nonlin);
	}

	// Print out radii in case of quartiles only.

	if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

	    int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	    real *r_lagr = new real[n_lagr];
	    getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);
	    cerr << "r_lagr =";
	    for (int i = 0; i < n_lagr; i++) cerr << " " << r_lagr[i];
	    cerr << endl;

	}

	// put_dyn(b);
	delete b;
    }
}

#endif

