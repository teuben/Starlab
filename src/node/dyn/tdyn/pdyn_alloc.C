
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// pdyn_alloc.C:  Memory management for the pdyn class (only!).
////
//// Options: none

// Externally visible functions:
//
//	new_pdyn
//	new_pdyn2
//	alloc_pdyn

#include "pdyn.h"

// Note: pdyn_mem is the base node of the pdyn management scheme; if it
// is non-NULL, we are doing our own memory management.

static int  n_pdyn_mem = 0;		// total number of available nodes

// Static arrays (pdyn, flag, pointer), length n_pdyn_mem, index i.
// Could be orgainzed as an array of structures, but this seems simpler.

static pdyn *pdyn_mem  = NULL;		// array of all nodes
static bool *pdyn_free = NULL;		// true iff node #i is available
static int  *pdyn_ptr  = NULL;		// pointer into pdyn_alloc

// Management array (pointer into the static array of nodes), index j.

static int  n_alloc = 0;		// number of allocated nodes
static int  *pdyn_alloc = NULL;		// nodes 0 to n_alloc-1 are allocated,
					// the remainder are free

// Pointers are such that	pdyn_alloc[pdyn_ptr[i]] = i
//				pdyn_ptr[pdyn_alloc[j]] = j

void dealloc_all()
{
    // Quick dealloc; don't clear stories or bases (careful!).

    n_alloc = 0;
    for (int i = 0; i < n_pdyn_mem; i++) {
	pdyn_free[i] = true;
	pdyn_ptr[i] = i;
	pdyn_alloc[i] = i;
    }
}

bool initialize_pdyn(int n,
		     hbpfp the_hbpfp,	// default = new_hydrobase
		     sbpfp the_sbpfp,	// default = new_starbase
		     bool use_stories)	// default = true
{
    // Create a block of n contiguous pdyn objects, with associated
    // pointers and management arrays.

    n_pdyn_mem = 0;
    pdyn_mem = new pdyn[n](the_hbpfp, the_sbpfp, use_stories);
    pdyn_free = new bool[n];
    pdyn_ptr = new int[n];
    pdyn_alloc = new int[n];

    if (pdyn_mem && pdyn_free && pdyn_ptr && pdyn_alloc) {
	n_pdyn_mem = n;
	dealloc_all();
	    
	cerr << "initialize_pdyn: created " << n
	     << " pdyn elements; pdyn_mem = " << pdyn_mem << endl;

	return true;

    } else
	return false;
}

// alloc_pdyn: replacement for new pdyn once static array has 
//	      been built (arguments from pdyn constructor)

pdyn *alloc_pdyn(hbpfp the_hbpfp,	// default = new_hydrobase
		 sbpfp the_sbpfp,	// default = new_starbase
		 bool use_stories,	// default = true
		 bool reinit)		// default = true
{
    if (pdyn_mem) {

	// Must create/delete *bases and stories as necessary.  Also, must
	// do all (node/dyn/...) initialization ourselves if reinit not set.

	// Get the next free node.  Should extend the array on overflow...

	if (n_alloc >= n_pdyn_mem) err_exit("pdyn array overflow");

	// Next available node is i = pdyn_alloc[n_alloc].  Allocate
	// it and adjust all pointers accordingly.

	int i = pdyn_alloc[n_alloc];
	pdyn_free[i] = false;
	pdyn_ptr[i] = n_alloc++;

	pdyn *p = pdyn_mem + i;

	if (reinit) {

	    // Initialize p using the same functions as the *dyn constructors.

	    p->pdyn_init();
	    p->dyn_init();
	    p->node_init();

	    // Optional (must delete if not wanted):

	    if (!use_stories)
		p->rmstory();
	    else
		p->node_set_stories(use_stories);

	    if (!the_hbpfp)
		p->rmhydrobase();
	    else
		p->node_set_hbase(the_hbpfp);

	    if (!the_sbpfp)
		p->rmstarbase();
	    else
		p->node_set_sbase(the_sbpfp);
	}

	return p;

    } else

	// The standard constructor:

	return new pdyn(the_hbpfp, the_sbpfp, use_stories);
}

// dealloc_pdyn: replacement for delete pdyn once static array has
//		 been built

void dealloc_pdyn(pdyn *p)
{
    if (p) {

	if (pdyn_mem) {

	    // Return p to the free list.  Clean up as though we were
	    // deleting the node, although this may not actually be necessary.

	    int i = p - pdyn_mem;		// easy computation of index!
	    pdyn_free[i] = true;
	    int j = pdyn_ptr[i];

	    // Location j in the management array is now free.
	    // Reorder the array to place the free slot at the end:
	    // swap elements j and n_alloc-1, then reduce n_alloc.

	    int itmp = pdyn_alloc[n_alloc-1];

	    pdyn_alloc[--n_alloc] = i;
	    pdyn_ptr[i] = n_alloc;

	    pdyn_alloc[j] = itmp;
	    pdyn_ptr[itmp] = j;

	    // cerr << "deleting " << p << " = " << p->format_label()
	    //      << ",  i = " << i << endl;

	    // pdyn destructor (no action):

	    // dyn destructor:

	    p->rmkepler();

	    // node destructor:

	    p->set_invalid();
	    p->clear_name();
	    if (p->is_root()) p->set_root(NULL);
	    p->rmstory();
	    p->rmhydrobase();
	    p->rmstarbase();

	} else

	    // The standard destructor:

	    delete p;
    }
}

// new_pdyn: used by the pdyn version of get_node()

node * new_pdyn(hbpfp the_hbpfp,
		sbpfp the_sbpfp,
		bool use_stories)	// default = true
{
    return (node *) new pdyn(the_hbpfp, the_sbpfp, use_stories);
}

// new_pdyn2: replacement for new_pdyn using our own memory management

node * new_pdyn2(hbpfp the_hbpfp,
		 sbpfp the_sbpfp,
		 bool use_stories)	// default = true
{
    return alloc_pdyn(the_hbpfp, the_sbpfp, use_stories);
}
