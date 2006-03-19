
//// set_radius:  set radii of specified particle(s) inthe  input snapshot
////              to the specified values. The function sets both member
////              data and R_eff in the dyn story.  Specified radii are
////              used unaltered -- i.e. no scaling to physical parameters
////              is performed.
////
////              *** Note new command-line arguments as of 3/06. ***
////
//// Options:    -A    set selected radius for all stars
////             -l    specify label of next particle to modify [no default]
////             -m    specify mass (or mass range) of particles to modify
////             -M    set radii for all stars using R = M^0.8
////             -r    specify new radius for particle [no default]
////
//// Usage:  set_radius -l l1 -r radius1 -l l2 -r radius2 ...
////    or:  set_radius -m m1 -r radius1 -m m2l m2u -r radius2 ...
////    or:  set_radius -A radius
////    or:  set_radius -M
////
//// Note that the label or mass range must be specified before the radius.

// Simon Portegies Zwart, MIT Oct 2000
// Ernest Mamikonyan, Drexel, Jul 2005
// Steve McMillan, Drexel, Mar 2006

#include "_dyn_.h"

#ifdef TOOLBOX

// leaf_with_label: Return a pointer to the leaf with the specified name.
//                  Function is very similar to node::node_with_name()...

local _dyn_* leaf_with_label(_dyn_* b, char* label)
{
    for_all_leaves(_dyn_, b, bi) {
	if (bi->get_name() != NULL) {
	    if (strcmp(bi->get_name(), label) == 0)
		return bi;
	} else if (bi->get_index() >= 0) {
	    char id[64];
	    sprintf(id, "%d", bi->get_index());
	    if (strcmp(id, label) == 0)
		return bi;
	}
    }
    return NULL;
}

local void set_radius_by_label(_dyn_* b, char* label, real radius)
{
    // Locate the particle to set.

    _dyn_* bi = leaf_with_label(b, label);

    if (bi == NULL)
	cerr << "Warning: particle \"" << label << "\" not found.\n";
    else
	bi->set_radius(radius, true);		// set radius and R_eff
}

local void set_radius_all(_dyn_* b, real radius)
{
    for_all_leaves(_dyn_, b, bi)
	bi->set_radius(radius, true);		// set radius and R_eff
}

local void set_radius_by_mass(_dyn_* b)
{
    real msf, lsf, tsf;
    get_physical_scales(b, msf, lsf, tsf);

    // Set the radius in units of the cluster virial radius.

    for_all_leaves(_dyn_, b, bi)
	bi->set_radius(pow(bi->get_mass(),0.8)*6.96e8/3.086e16, true);
}

local void set_radius_by_mass_range(_dyn_* b, real ml, real mu, real radius)
{
    for_all_leaves(_dyn_, b, bi) {
	real mi = bi->get_mass();
	if (mi >= ml && mi <= mu)
	    bi->set_radius(radius, true);
    }
}

int main(int argc, char ** argv)
{
    char label[64];
    label[0] = '\0';

    check_help();
    pgetopt(argc, argv, "", "$Revision$", _SRC_);

    _dyn_* b;
    b = get__dyn_();

    b->log_history(argc, argv);

    // Parse the command line by hand, modifying the system as we go.

    int i = 0;
    real ml = 0, mu = 0;
    int which_spec = 0;

    while (++i < argc)
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {

		case 'A': set_radius_all(b, (real)atof(argv[++i]));
			  break;

		case 'l': strcpy(label, argv[++i]);
			  which_spec = 1;
			  break;

		case 'm': ml = mu = atof(argv[++i]);
			  if (argv[i+1][0] != '-') mu = atof(argv[++i]);
			  which_spec = 2;
			  break;

		case 'M': set_radius_by_mass(b);
			  break;

		case 'r': if (which_spec == 1)
			      set_radius_by_label(b, label,
						  (real)atof(argv[++i]));
			  else if (which_spec == 2)
			      set_radius_by_mass_range(b, ml, mu,
						       (real)atof(argv[++i]));
			  else
			      cerr << "Must specify label or mass range"
				   << " before radius." << endl;
			  which_spec = 0;
			  break;

		default:  get_help();
	    }

    put__dyn_(b);
    return 0;
}

#endif
