
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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
//   version 1:  May 1989   Piet Hut               email: piet@iassns.BITNET
//                           Institute for Advanced Study, Princeton, NJ, USA
//   version 2:  Dec 1992   Piet Hut  --  adapted to the new C++-based starlab
//   version 3:  Jul 1996   Steve McMillan & Jun Makino
//.............................................................................
//   non-local functions: 
//      compute_general_mass_radii
//      compute_mass_radii_quartiles
//      compute_mass_radii_percentiles
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

typedef  struct
{
    real  radius;
    real  mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return +1;
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return -1;
    else
        return 0;
}

//-----------------------------------------------------------------------------
//  compute_general_mass_radii  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

static real nonlin_masses[9] = {0.005, 0.01, 0.02, 0.05, 0.1,
				0.25, 0.5, 0.75, 0.9};

void compute_general_mass_radii(dyn * b, int nzones,
				bool nonlin,
				boolfn bf)
{
    if (nzones < 2) return;
    if (nonlin && nzones != 10) return;		// Special case

    // Note: nzones specifies the number of radial zones to consider.
    //	     However, we only store radii corresponding to the nzones-1 
    //	     interior zone boundaries.

    char lagr_string[64] = "geometric center";

    int n = 0;
    if (bf == NULL)
	n = b->n_daughters();
    else {
	for_all_daughters(dyn, b, bb)
	    if ((*bf)(bb)) n++;
    }

    // Use the density center if known and up to date (preferred).
    // Otherwise, use modified center of mass, if known and up to date.

    vec lagr_pos, lagr_vel;
    int which = get_std_center(b, lagr_pos, lagr_vel);
    if (which == 1)
	strcpy(lagr_string, "density center");
    else
	strcpy(lagr_string, "modified center of mass");

    rm_pair_ptr rm_table = new rm_pair[n];

    if (rm_table == NULL) {
	cerr << "compute_general_mass_radii: "
	     << "not enough memory left for rm_table\n";
	return;
    }

    // Set up an array of (radius, mass) pairs.  Also find the total
    // mass of all nodes under consideration.

    real total_mass = 0;
    int i = 0;

    vec dlagr_pos = lagr_pos - b->get_pos();

    for_all_daughters(dyn, b, bi) {
	if (bf == NULL || (*bf)(bi)) {
	    total_mass += bi->get_mass();
	    rm_table[i].radius = abs(bi->get_pos() - dlagr_pos);
	    rm_table[i].mass = bi->get_mass();
	    i++;
	}
    }

    // Sort the array by radius.  (Slightly wasteful, as get_std_center
    // will also sort the data if compute_mcom is called.)

    qsort((void *)rm_table, (size_t)i, sizeof(rm_pair), compare_radii);

    // Determine the Lagrangian radii.

    // cerr << "Determining Lagrangian radii 1" << endl << flush;

    real* mass_percent = new real[nzones-1];
    if (mass_percent == NULL) {
	cerr << "compute_general_mass_radii: "
	     << "not enough memory left for mass_percent\n";
	delete [] rm_table;
	return;
    }

    int k;
    for (k = 0; k < nzones-1; k++) {
        if (!nonlin) 
	    mass_percent[k] = ((k + 1) / (real)nzones) * total_mass;
	else
	    mass_percent[k] = nonlin_masses[k] * total_mass;
    }

    real *rlagr = new real[nzones-1];
    if (rlagr == NULL) {
	cerr << "compute_general_mass_radii: "
	     << "not enough memory left for r_lagr\n";
	delete [] rm_table;
	delete [] mass_percent;
	return;
    }
    real cumulative_mass = 0.0;
    i = 0;

    // cerr << "Determining Lagrangian radii 2" << endl << flush;

    for (k = 0; k < nzones-1; k++) {

        while (cumulative_mass < mass_percent[k])
	    cumulative_mass += rm_table[i++].mass;

	rlagr[k] = rm_table[i-1].radius;
    }

    // cerr << "writing stories" << endl << flush;

    // Place the data in the root dyn story.

    if (bf == NULL)
	putiq(b->get_dyn_story(), "boolfn", 0);
    else
	putiq(b->get_dyn_story(), "boolfn", 1);

    putiq(b->get_dyn_story(), "n_nodes", n);
    putrq(b->get_dyn_story(), "lagr_time", b->get_system_time());
    putvq(b->get_dyn_story(), "lagr_pos", lagr_pos);
    putvq(b->get_dyn_story(), "lagr_vel", lagr_vel);
    putsq(b->get_dyn_story(), "pos_type", lagr_string);
    putiq(b->get_dyn_story(), "n_lagr", nzones-1);
    putra(b->get_dyn_story(), "r_lagr", rlagr, nzones-1);

    delete [] mass_percent;
    delete [] rm_table;
    delete [] rlagr;
}

// Convenient synonyms:

void  compute_mass_radii_quartiles(dyn * b)
{
    compute_general_mass_radii(b, 4);
}

void  compute_mass_radii_percentiles(dyn * b)
{
    compute_general_mass_radii(b, 10);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  compute_mass_radii() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    int n = 0;
    bool  c_flag = false;      // if TRUE: a comment given on command line
    bool  t_flag = false;      // if TRUE: percentiles rather than quartiles
    bool  nonlin = false;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:n:st";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 's': n = 10;
                      nonlin = true;
	    case 't': t_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);

        b->log_history(argc, argv);

	if (t_flag)
	    compute_mass_radii_percentiles(b);
	else {
	    if (n <= 1)
		compute_mass_radii_quartiles(b);
	    else
		compute_general_mass_radii(b, n, nonlin);
	}

	// Print out radii in case of quartiles only.

	if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

	    int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	    real *r_lagr = new real[n_lagr];
	    getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);
	    cerr << "r_lagr =";
	    for (int i = 0; i < n_lagr; i++) cerr << " " << r_lagr[i];
	    cerr << endl;
	    delete [] r_lagr;
	}

	put_dyn(b);
	rmtree(b);
    }
}

#endif

// endof: lagrad.C
