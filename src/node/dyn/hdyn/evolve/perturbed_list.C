
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions for maintaining the list of (top-level) perturbed binaries.
// "Maintenance" calls come from kira.C, kira_stellar.C, and kira_escape.C.
// Major use of the list is in hdyn_ev.C (correct_acc_and_jerk)
//
//						Steve McMillan, 3/99
//
// Externally visible functions:
//
//	bool hdyn::is_perturbed_cpt
//	bool hdyn::is_perturbed_cm
//	void hdyn::reconstruct_perturbed_list
//	void hdyn::dump_perturbed_list
//	int  hdyn::get_index_on_perturbed_list
//	bool hdyn::on_perturbed_list
//	void hdyn::check_perturbed_list
//	void hdyn::add_to_perturbed_list
//	void hdyn::remove_from_perturbed_list

#include "hdyn.h"

// Define "perturbed", for use everywhere perturbed_list is used.

bool hdyn::is_perturbed_cpt()
{
    return (!get_kepler() || !fully_unperturbed);
}

bool hdyn::is_perturbed_cm()
{
    return (is_parent() && get_oldest_daughter()->is_perturbed_cpt());
}

// Initial asumption was that, in use, the list would be checked
// frequently, but updated relatively rarely.  Thus, we expected to
// use an ordered list for convenience of searching, even though that
// means updates may be significantly more expensive -- monitor this!

// In fact, the only current use of these functions is to check entire
// list in correct_acc_and_jerk.  In this case, individual access time
// is unimportant, since we never do it, so we should perhaps optimize
// the list maintenance cost (Steve, 8/00).

//#define FAST_ACCESS true	// faster individual access, slower maintenance
#define FAST_ACCESS false

// Note from Steve (8/00):  Using the FAST_ACCESS code requires that the
// list be sorted.  Currently, we sort by address.  However, this means
// that the details of the perturber correction depend on where nodes lie
// in memory, making the results in general non-reproducible.  Differences
// crop up typically in the last bit because the order in which corrections
// are applied can vary from run to run.
//
// **** For this and the above reason, FAST_ACCESS is presently turned off.
// **** If it is restored, we must find another quantity to use for sorting
// **** purposes -- perhaps define a unique index for all nodes?

local int compare_ptr(const void * pi, const void * pj)
{
    hdynptr bi = *((hdynptr*)pi);
    hdynptr bj = *((hdynptr*)pj);

    if (bi > bj)
        return +1;
    else if (bi < bj)
        return -1;
    else
        return 0;
}

void hdyn::reconstruct_perturbed_list()
{
    if (!options->use_perturbed_list) return;

    if (perturbed_list) delete [] perturbed_list;

    perturbed_list = new hdynptr[get_root()->n_leaves()];
    n_perturbed = 0;

    for_all_daughters(hdyn, get_root(), bi)
	if (bi->is_perturbed_cm()) {

	    perturbed_list[n_perturbed++] = bi;

	    //cerr << "Added node " << bi << " " << bi->format_label()
	    //     << " as perturbed_list[" << n_perturbed-1
	    //	   << "] at time " << bi->get_system_time()
	    //	   << endl;

	}

    if (FAST_ACCESS) {

	// Sort the list -- O(N log N) ops.

	qsort(perturbed_list, n_perturbed, sizeof(hdynptr), compare_ptr);
    }

    cerr << endl
	 << "rebuilt perturbed_list; n_perturbed = " << n_perturbed
	 << endl;
}

void hdyn::dump_perturbed_list()
{
    if (n_perturbed <= 0) return;

    cerr << "perturbed_list (n_perturbed = " << n_perturbed << "):"
	 << endl;

    for (int i = 0; i < n_perturbed; i++)
	cerr << "    " << i << "  " << perturbed_list[i]
	     << " " << perturbed_list[i]->format_label() << endl;
}

int hdyn::get_index_on_perturbed_list(bool debug) // default = false
{
    // If 'this' is not on the list, but is within the range defined
    // by the list, return the element immediately *below* it.  If this
    // is below the bottom of the list, return -1; if this is above the
    // top element, return n_perturbed.
    //
    // Note from Steve, 8/00:  The dual nature of the return from this
    // function is never used, as we always test for < 0 and >= n_perturbed
    // in the calling function.  Thus, the convention used here is not wrong,
    // but it is redundant and slightly inefficient...

    if (!FAST_ACCESS) {

	// Simplest O(N) version:

	bool below = true;
	for (int i = 0; i < n_perturbed; i++) {
	    if (perturbed_list[i] == this)
		return i;
	    else if (below && perturbed_list[i] > this)
		below = false;
	}

	if (below)
	    return -1;
	else				// redundant -- return value of -1
	    return n_perturbed;		// in all out of range cases is OK

    } else {

	// Assume the list is ordered and do an O(log N) binary search.

	if (debug)
	    cerr << "In get_index_on_perturbed_list for " << format_label()
		 << " at time " << system_time << endl;

	if (this < perturbed_list[0])
	    return -1;
	else if (this > perturbed_list[n_perturbed-1])
	    return n_perturbed;

	int ilow = 0, ihigh = n_perturbed - 1;

	while (ilow < ihigh) {

	    int imid = (ilow + ihigh)/2;

	    // Logic can probably be improved...

	    if (imid == ilow) {
		if (perturbed_list[ihigh] == this)
		    return ihigh;
		else
		    return ilow;
	    }

	    if (perturbed_list[imid] == this)
		return imid;
	    else if (perturbed_list[imid] < this)
		ilow = imid;
	    else
		ihigh = imid;

	    if (debug) {
		PRC(this), PRC(ilow), PRL(ihigh);
		PRC(perturbed_list[ilow]), PRL(perturbed_list[ihigh]);
	    }
	}

	return ilow;			// redundant
    }
}

bool hdyn::on_perturbed_list()
{
    int i = get_index_on_perturbed_list();
    return (i >= 0 && i < n_perturbed && perturbed_list[i] == this);
}

void hdyn::check_perturbed_list()
{
    // Note: checking and reconstruction are O(N log N) operations
    //	     if FAST_ACCESS is true; otherwise checking is O(N^2),
    //	     reconstruction is O(N).

    if (!perturbed_list) {

	cerr << endl
	     << "check_perturbed_list:  perturbed_list = NULL!"
	     << endl;

    } else {

	int i;
	bool reconstruct_list = false;

	// Check that all perturbed top-level binaries are on the list.

	for_all_daughters(hdyn, get_root(), bi)
	    if (bi->is_perturbed_cm()) {
		if ((i = bi->get_index_on_perturbed_list()) < 0
		    || i >= n_perturbed
		    || perturbed_list[i] != bi) {
		    cerr << "check_perturbed_list:  " << bi->format_label()
			 << " not on perturbed_list" << endl;
		    reconstruct_list = true;
		}
	    }

	// Check that all list entries are perturbed top-level binaries.

	for (i = 0; i < n_perturbed; i++) {
	    if (!perturbed_list[i]->is_perturbed_cm()) {
		cerr << "check_perturbed_list:  entry " << i
		     << " " << perturbed_list[i]->format_label()
		     << " is not perturbed" << endl;
		 reconstruct_list = true;
	     }
	}

	// Check that the list is properly ordered.

	for (i = 0; i < n_perturbed-1; i++) {
	     if (perturbed_list[i] >= perturbed_list[i+1]) {
		 cerr << "perturbed_list corrupted at i = " << i << endl;
		 reconstruct_list = true;
	     }
	}

	// Fix any problems.

	if (reconstruct_list) {
	    // dump_perturbed_list();
	    reconstruct_perturbed_list();
	    cerr << "perturbed_list reconstructed" << endl;
	}
    }
}

void hdyn::add_to_perturbed_list(int id) // default = 0
{
    if (!options->use_perturbed_list) return;

    if (!perturbed_list) {
	perturbed_list = new hdynptr[get_root()->n_leaves()];
	n_perturbed = 0;
    }

    // This is an O(N) operation if we have to place a node in the middle
    // of the list.  Adding to the end is O(1), but the test for membership
    // is O(N) instead of O(log N) in the case FAST_ACCESS = false, so there
    // isn't really much advantage...

    int i = get_index_on_perturbed_list();

    if (i >= 0 && i < n_perturbed && perturbed_list[i] == this) {

	cerr << "add_to_perturbed_list: "
	     << format_label() << " already on list"
	     << " at location " << i << endl;

	check_perturbed_list();

    } else {

	if (FAST_ACCESS) {

	    // Redefine i to be the insertion point for 'this'.  Note
	    // the convention that, if 'this' lies above the range of the
	    // present list, i = n_perturbed.

	    if (i < n_perturbed) i++;

	    // Make room for the new node.

	    for (int j = n_perturbed; j > i; j--)
		perturbed_list[j] = perturbed_list[j-1];

	} else

	    i = n_perturbed;

	perturbed_list[i] = this;
	n_perturbed++;

	if (diag->report_adjust_perturbed_list) {

	    cerr << "added node " << this << " " << format_label()
		 << " to perturbed_list at time " << system_time
		 << endl
		 << "    index = " << i
		 << "  n_perturbed = " << n_perturbed;

	    if (id > 0) cerr << " (ID = " << id << ")";

	    cerr << endl;
	}
    }
}

void hdyn::remove_from_perturbed_list(int id) // default = 0
{
    if (!options->use_perturbed_list) return;

    // Same O(N) method, regardless of FAST_ACCESS.

    for_all_nodes(node, get_root(), bb)
	if (bb!=bb->get_starbase()->get_node()) {
	    PRC(bb->format_label());
	    PRC(bb);PRL(bb->get_starbase()->get_node());
	}

    if (!perturbed_list || n_perturbed <= 0) return;

    int i = get_index_on_perturbed_list();

    if (i < 0 || i >= n_perturbed || perturbed_list[i] != this) {

	cerr << "remove_from_perturbed_list(" << id
	     << "): " << format_label() << " not on list" << endl;
	if (i >= 0 && i < n_perturbed)
	  cerr << "    perturbed_list[i] = "
	       << perturbed_list[i]->format_label() << endl;

	check_perturbed_list();

    } else {

	// Remove 'this' and contract the list.

	for (int j = i; j < n_perturbed-1; j++)
	    perturbed_list[j] = perturbed_list[j+1];

	n_perturbed--;

	if (diag->report_adjust_perturbed_list) {

	    cerr << "removed node " << this << " " << format_label()
		 << " from perturbed_list at time " << system_time
		 << endl
		 << "    index = " << i
		 << "  n_perturbed = " << n_perturbed;

	    if (id > 0) cerr << " (ID = " << id << ")";

	    cerr << endl;
	}
    }

    for_all_nodes(node, get_root(), bb)
	if (bb != bb->get_starbase()->get_node()) {
	    PRC(bb->format_label());
	    PRC(bb);PRL(bb->get_starbase()->get_node());
	}
}
