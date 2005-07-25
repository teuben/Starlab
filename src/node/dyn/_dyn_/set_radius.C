
//// set_radius:  set radii of specified particle(s) in input snapshot
////           to specified values.
////
//// Usage:  set_radius -l l1 -m radius1 -l l2 -m radius2 ...
////    or:  set_radius -A radius
////    or:  set_radius -M
////
//// Options:    -A    set selected radius for all stars
////             -l    specify label of next particle to modify [no default]
////             -M    set radii for all stars using R = M^0.8
////             -m    specify new radius for particle [no default]

// Simon Portegies Zwart, MIT Oct 2000

#include "_dyn_.h"

#ifdef TOOLBOX

// leaf_with_label: Return a pointer to the leaf with the specified name.

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
    // Locate the particle to be split.

    _dyn_* bi = leaf_with_label(b, label);

    if (bi == NULL)
	cerr << "Warning: particle \"" << label << "\" not found.\n";
    else
	bi->set_radius(radius);

}

local void set_radius(_dyn_* b, real radius)
{
    // Locate the particle to be split.

    for_all_leaves(_dyn_, b, bi) {
      if (bi == NULL)
	cerr << "Warning: no particle found.\n";
      else
	bi->set_radius(radius);
    }

}

local void set_radius_by_mass(_dyn_* b)
{

    real msf, lsf, tsf;
    get_physical_scales(b, msf, lsf, tsf);

    for_all_leaves(_dyn_, b, bi) {
      if (bi == NULL)
	cerr << "Warning: no particle found.\n";
      else
	// set the radius in units of the cluster virial radius
	bi->set_radius(pow(bi->get_mass(),0.8)*6.96e8/3.086e16);
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
    while (++i < argc)
	if (argv[i][0] == '-')
	    switch (argv[i][1]) {

		case 'A': set_radius(b, (real)atof(argv[++i]));
			  break;

		case 'l': strcpy(label, argv[++i]);
			  break;

		case 'M': set_radius_by_mass(b);
			  break;

		case 'm': set_radius_by_label(b, label, (real)atof(argv[++i]));
			  break;

		default:  get_help();
	    }

    put__dyn_(b);
    return 0;
}

#endif
