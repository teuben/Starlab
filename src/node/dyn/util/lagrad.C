
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Compute Lagrangian (mass) radii for an N-body system.  The radii
//// are stored in the system (root) dyn story.  If a current density
//// center is found in the root dyn story, the Lagrangian radii are
//// calculated relative to it.  Otherwise, if a valid center of mass
//// is found, it is used.  If neither center is found, the geometric
//// center is used.
////
//// Usage: lagrad [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -n    specify number of Lagrangian zones (linear in mass) [4]
////              -s    use "special" nonlinear binning:
////                        0.005, 0.01, 0.02, 0.05,
////                        0.1, 0.25, 0.5, 0.75, 0.9
////              -t    same as -n 10
////
//// Written by Piet Hut, Steve McMillan, and Jun Makino.
////
//// Report bugs to starlab@sns.ias.edu.

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
//				    Return true iff the radii are usable.
//-----------------------------------------------------------------------------

static real nonlin_masses[9] = {0.005, 0.01, 0.02, 0.05, 0.1,
				0.25, 0.5, 0.75, 0.9};

bool compute_general_mass_radii(dyn * b, int nzones,
				bool nonlin,
				boolfn bf)
{
    if (nzones < 2) return false;
    if (nonlin && nzones != 10) return false;		// special case

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

    if (n <= 0) {				// possible in restricted cases
	cerr << endl
	     << "    compute_general_mass_radii: no stars in this category\n";
	return false;
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
	cerr << endl
	     << "    compute_general_mass_radii: "
	     << "not enough memory left for rm_table\n";
	return false;
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

    // cerr << "    determining Lagrangian radii 1" << endl << flush;

    real* mass_percent = new real[nzones-1];
    if (mass_percent == NULL) {
	cerr << endl
	     << "    compute_general_mass_radii: "
	     << "not enough memory left for mass_percent\n";
	delete [] rm_table;
	return false;
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
	cerr << endl
	     << "    compute_general_mass_radii: "
	     << "not enough memory left for r_lagr\n";
	delete [] rm_table;
	delete [] mass_percent;
	return false;
    }
    real cumulative_mass = 0.0;
    i = 0;

    // cerr << "    determining Lagrangian radii 2" << endl << flush;

    for (k = 0; k < nzones-1; k++) {

        while (cumulative_mass < mass_percent[k])
	    cumulative_mass += rm_table[i++].mass;

	rlagr[k] = rm_table[i-1].radius;
    }

    // cerr << "    writing stories" << endl << flush;

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

    return true;
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



// Print_lagrangian_radii and supporting functions.

local bool testnode(dyn * b)
{
    return (b->is_top_level_node());
}

local bool binary_fn(dyn * b)
{
    return (b->get_oldest_daughter());
}

local bool single_fn(dyn * b)
{
    return (getiq(b->get_log_story(), "mass_doubled") == 0);
}

local bool double_fn(dyn * b)
{
    return (getiq(b->get_log_story(), "mass_doubled") == 1);
}

#include <vector>
#include <algorithm>

static real cutoff_mass = 0;

void set_lagr_cutoff_mass(dyn *b, real f)
{
    // Determine a cutoff for the "most massive" single stars, as follows.
    // Start with the f-th percentile mass, then move up the list toward
    // the more massive end until an increase in mass occurs, and use
    // that new mass.  This will choose a mass very close to the f-th
    // percentile in case of a continuous mass function, and should pick
    // out the next group up if have discrete mass groups.

    vector<real> *m = new vector<real>;
    for_all_daughters (dyn, b, bb)
	if (bb->is_leaf()) m->push_back(bb->get_mass());

    sort(m->begin(), m->end());

    int nf = (int)(f*m->size()) - 1, i = nf;
    while (i < m->size() && (*m)[i] == (*m)[nf]) i++;
    if (i < m->size())
	cutoff_mass = (*m)[i];
    else
	cutoff_mass = VERY_LARGE_NUMBER;
}

real get_lagr_cutoff_mass()
{
    return cutoff_mass;
}

local bool massive_fn(dyn * b)
{
    return (b->get_mass() >= cutoff_mass);
}

real print_lagrangian_radii(dyn* b,
			    int which_lagr,	// default = 2 (nonlinear)
			    bool verbose,	// default = true
			    int which_star,	// default = 0 (all stars)
			    bool noprint)	// default = false (print!)
{
    bool nonlin = false;

    real rhalf = 0;
    int ihalf;

    int nl, indent;
    if (which_lagr == 0) {
	nl = 4;
	indent = 15;
	ihalf = 1;
    } else if (which_lagr == 1) {
	nl = 10;
	indent = 20;
	ihalf = 4;
    } else {
	nl = 10;
	indent = 26;
	nonlin = true;
	ihalf = 6;
    }

    bool status = false;

    if (which_star == 0)
	status = compute_general_mass_radii(b, nl, nonlin);
    else if (which_star == 1)
	status = compute_general_mass_radii(b, nl, nonlin, binary_fn);
    else if (which_star == 2)
	status = compute_general_mass_radii(b, nl, nonlin, single_fn);
    else if (which_star == 3)
	status = compute_general_mass_radii(b, nl, nonlin, double_fn);
    else if (which_star == 4)
	status = compute_general_mass_radii(b, nl, nonlin, massive_fn);

    if (status && find_qmatch(b->get_dyn_story(), "n_lagr")
	&& getrq(b->get_dyn_story(), "lagr_time") == b->get_system_time()) {

	// Assume that lagr_pos has been properly set if n_lagr is set and
	// lagr_time is current.  Redundant now that status is returned.

	vec lagr_pos = getvq(b->get_dyn_story(), "lagr_pos");

	if (verbose && !noprint) {
	    cerr << endl << "  Lagrangian radii relative to ("
		 << lagr_pos << "):" << endl;
	    if (find_qmatch(b->get_dyn_story(), "pos_type"))
		cerr << "                               ( = "
		     << getsq(b->get_dyn_story(), "pos_type")
		     << ")" << endl;
	    if (which_lagr == 0)
		cerr << "    quartiles: ";
	    else if (which_lagr == 1)
		cerr << "    10-percentiles: ";
	    else
		cerr << "    selected percentiles: ";
	}

        int n_lagr = getiq(b->get_dyn_story(), "n_lagr");  // should be nl - 1
	real *r_lagr = new real[n_lagr];

	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	for (int k = 0; k < n_lagr; k += 5) {
	    if (!noprint) {
		if (k > 0) {
		    cerr << endl;
		    for (int kk = 0; kk < indent; kk++) cerr << " ";
		}
		for (int i = k; i < k+5 && i < n_lagr; i++)
		    cerr << " " << r_lagr[i];
	    }
	    cerr << endl << flush;
	}

	rhalf = r_lagr[ihalf];
	delete [] r_lagr;
    }

    return rhalf;
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

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
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
