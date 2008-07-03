
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                      //            _\|/_
//=======================================================//              /|\ ~

//// Rename all CM nodes bottom up to reflect the standard Starlab
//// naming scheme (a,b).  By default, we don't touch existing
//// numeric names, but any non-numeric name are overwritten and
//// all names are checked for uniqueness.
////
//// NB: New code is largely untested.
////
//// Options:
////
////          -c    don't check for uniqueness [check]
////          -p    don't preserve existing numeric names [preserve]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   Steve McMillan, Sep 2003, Sep 2005

#ifdef TOOLBOX

#include "node.h"

typedef  struct
{
    int index;
    node *b;
} ib_pair, *ib_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_indices  --  compare the indices of two particles
//-----------------------------------------------------------------------------

local int compare_indices(const void * pi, const void * pj)
{
    if (((ib_pair_ptr)pi)->index > ((ib_pair_ptr)pj)->index)
        return +1;
    else if (((ib_pair_ptr)pi)->index < ((ib_pair_ptr)pj)->index)
        return -1;
    else
        return 0;
}

local void check_unique_indices(node *b)
{
    const char *func = "check_unique_indices";

    // Check for uniqueness in node indices (potentially expensive),
    // and correct any duplicates.

    // Old-fashioned sort...

    int n = 0;
    for_all_leaves(node, b, bi) n++;

    ib_pair_ptr ib_table = new ib_pair[n];

    if (ib_table == NULL) {
	cerr << endl << "    " << func
		     << ": not enough memory for ib_table\n";
	return;
    }

    // Set up an array of (index, nodeptr) pairs.
    // Only look at nodes with index > 0.

    int i = 0, imax = 0;
    for_all_leaves(node, b, bi) {
	int index = bi->get_index();
	if (index > 0) {
	    ib_table[i].index = index;
	    ib_table[i++].b = bi;
	    if (imax < index) imax = index;
	}
    }

    // Sort the array by index.

    qsort((void *)ib_table, (size_t)i, sizeof(ib_pair), compare_indices);

    // Go down the list looking for duplicates.

    int lasti = ib_table[0].index, ndup = 0;

    for (int j = 1; j < i; j++) {
	if (ib_table[j].index == lasti) {

	    // Must replace this index.

	    (ib_table[i].b)->set_index(++imax);
	    ndup++;

	} else

	    lasti = ib_table[j].index;
    }

    cerr << func << ": " << ndup << " duplicate indices renamed" << endl;
    delete [] ib_table;
}

local void check_valid_indices(node *b)
{
    const char *func = "check_valid_indices";

    // Check for validity of indices and apply changes as needed.
    // Code may duplicate some of the operations of the previous function.

    int imax = 0;
    for_all_leaves(node, b, bi) {
	int index = bi->get_index();
	if (index > 0) {
	    if (imax < index) imax = index;
	}
    }

    int nmod = 0;
    for_all_leaves(node, b, bi) {
	int index = bi->get_index();
	if (index <= 0) {
	    bi->set_index(++imax);
	    nmod++;
	}
    }

    cerr << func << ": " << nmod << " invalid indices renamed" << endl;
}

// Should probably be moved to node and made global (Steve, 9/03).

local void name_parent_from_components(node *od, node *yd)
{
    char name1[256], name[256];
    strcpy(name1, od->format_label());
    sprintf(name, "(%s,%s)", name1, yd->format_label());
    od->get_parent()->set_name(name);
    od->get_parent()->set_index(-1);
}

int main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "cp";

    bool preserve = true, check = true;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'c': check = !check;
		      break;
	    case 'p': preserve = !preserve;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    node *b;
    while (b = get_node()) {
	b->log_history(argc, argv);

	if (preserve) {
	    if (check) check_unique_indices(b);
	    check_valid_indices(b);
	}

	int i = 0;
	for_all_leaves(node, b, bb) {
	    i++;
	    if (!preserve) bb->set_index(i);
	    if (bb->is_low_level_node()) {
		node *sis = bb->get_elder_sister();
		if (!sis || sis->is_parent()) {
		    if (!sis)
			name_parent_from_components(bb,
						    bb->get_younger_sister());
		    else
			name_parent_from_components(sis, bb);
		}
	    }
	}

	put_node(b);
	rmtree(b);
    }
}
#endif
