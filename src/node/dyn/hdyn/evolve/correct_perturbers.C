
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Remove escaping or destroyed nodes from all perturber lists.
//
// Externally visible functions:
//
//	void correct_perturber_lists
//	void expand_perturber_lists

#include "hdyn.h"

int counter = 0;

local void check_format_and_print(hdyn* bi)
{
    if (counter%5 == 0) {
	cerr << endl << "  ";
    }
    cerr << "  " << bi->format_label();
    counter++;
}

local void check_and_correct_perturbers(hdyn* bi, hdyn ** ilist, int n,
					hdyn * cm)
{
    if (bi->get_perturber_list() == NULL) {
	cerr << endl << "warning: check_and_correct_perturbers: "
	     << bi->format_label()
	     << " has valid_perturbers but no perturber_list" << endl;
	return;
    }

    hdynptr* plist = bi->get_perturber_list();		// convenient shorthand

    bool first = true;
    int null_count = 0;

    for (int i = 0; i < bi->get_n_perturbers(); i++) {
	hdyn* p = plist[i];
	for (int j = 0; j < n; j++) {
	    if (ilist[j] == p) {

		// Replace first plist member found to be on ilist by cm
		// (which may be NULL).  Set all others NULL.

		if (first) {
		    if (bi->get_kira_diag()->report_correct_perturber_list)
			check_format_and_print(bi);
		    plist[i] = cm;
		} else
		    plist[i] = NULL;

		first = false;
		if (plist[i] == NULL) null_count++;

		break;		
	    }
	}
    }

    // Do a second pass through the list to remove NULLs.

    if (null_count > 0) {
	int offset = 0;
	for (int i = 0; i < bi->get_n_perturbers(); i++) {
	    if (offset > 0) plist[i-offset] = plist[i];
	    if (plist[i] == NULL) offset++;
	}
	bi->set_n_perturbers(bi->get_n_perturbers()-null_count);
    }
}

// correct_perturber_lists:  Remove all nodes on the list from any perturber
//			     list that may contain them, replacing the first
//			     by cm, if non-NULL.

void correct_perturber_lists(hdyn * b, hdyn ** list, int n,
			     hdyn * cm) // default = NULL
{
    if (n <= 0) return;

    if (b->get_kira_diag()->report_correct_perturber_list) {
	cerr << endl << "correct_perturber_lists at time "
	     << b->get_system_time() << ": ";
	PRL(ALLOW_LOW_LEVEL_PERTURBERS);

	cerr << "removing nodes:";
	for (int i = 0; i < n; i++) cerr << "  " << list[i]->format_label();
	cerr << ",  cm = ";
	if (cm)
	    cerr << cm->format_label();
	else
	    cerr << "NULL";
	cerr << flush;

	counter = 0;
    }

    if (ALLOW_LOW_LEVEL_PERTURBERS) {

	for_all_nodes(hdyn, b, bi)
	    if (bi->get_valid_perturbers())
		check_and_correct_perturbers(bi, list, n, cm);

    } else {

	for_all_daughters(hdyn, b, bi)
	    if (bi->get_valid_perturbers())
		check_and_correct_perturbers(bi, list, n, cm);
    }

    if (b->get_kira_diag()->report_correct_perturber_list) {
	if (counter == 0) cerr << "  (no corrections)";
	cerr << endl;
    }

    // cerr << "Leaving correct_perturber_lists." << endl << flush;
}

local int expand_perturber_list(hdyn* bi, hdyn* bb, bool verbose)
{
    if (bi->get_perturber_list() == NULL) {
	cerr << endl << "warning: expand_perturber_list: "
	     << bi->format_label()
	     << " has valid_perturbers but no perturber_list" << endl;
	return 0;
    }

    hdynptr* plist = bi->get_perturber_list();		// convenient shorthand
    int np = bi->get_n_perturbers();

    for (int i = 0; i < np; i++) {
	hdyn* p = plist[i];
	if (plist[i] == bb) {

	    if (verbose) check_format_and_print(bi);

	    // Should only be expanding unperturbed binaries.  Flag
	    // expansion of perturbed binaries.

	    if (!bb->get_oldest_daughter()->get_kepler()) {
		cerr << "warning: expand_perturber_list for node "
		     << bi->format_label() << ": ";
		cerr << bb->format_label() << " is perturbed" << endl;
	    }

	    if (np >= MAX_PERTURBERS - 1) {

		// List will overflow if expanded.  Make it invalid.

		bi->set_valid_perturbers(false);

	    } else {

		// Remove bb from the list...

		for (int j = i+1; j < np; j++)
		    plist[j-1] = plist[j];
		np--;

		// ...then add its components to the end.

		for_all_leaves(hdyn, bb, bbb)
		    plist[np++] = bbb;

		bi->set_n_perturbers(np);
	    }

	    return 1;
	}
    }
    return 0;
}

// expand_perturber_lists:  Expand the specified node into components on all
//			    perturber lists that contain it.

void expand_perturber_lists(hdyn* b, hdyn* bb,	// bb is CM node
			    bool verbose)	// default = false
{
    if (bb->is_leaf()) return;

    if (verbose) {
	cerr << endl
	     << "expand_perturber_lists at time "
	     << b->get_system_time() << " to resolve "
	     << bb->format_label() << ":";
	counter = 0;
    }

    int np = 0;
    if (ALLOW_LOW_LEVEL_PERTURBERS) {
	for_all_nodes(hdyn, b, bi)
	    if (bi != bb && bi->get_valid_perturbers())
		np += expand_perturber_list(bi, bb, verbose);
    } else {
	for_all_daughters(hdyn, b, bi)
	    if (bi != bb && bi->get_valid_perturbers())
		np += expand_perturber_list(bi, bb, verbose);
    }

    if (verbose) {
	if (counter == 0) cerr << "  (no changes)";
	cerr << endl;
    }
}
