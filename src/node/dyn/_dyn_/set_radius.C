
//// setradius:  set radiuses of specified particle(s) in input snapshot
////           to specified values.
////
//// Usage:  setradius -l l1 -m radius1 -l l2 -m radius2 ...
////    or:  setradius -A radius
////
//// Options:    -A    set selected radius for all stars
////             -l    specify label of next particle to modify [no default]
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

void main(int argc, char ** argv)
{
    char label[64];
    label[0] = '\0';

    check_help();

    _dyn_* b;
    b = get__dyn_(cin);

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

		case 'm': set_radius_by_label(b, label, (real)atof(argv[++i]));
			  break;

		default:  get_help();
	    }

    put__dyn_(cout, *b);
}

#endif
