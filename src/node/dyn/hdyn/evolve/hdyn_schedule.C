
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  hdyn_schedule.C: scheduling routines
//.............................................................................
//    version 1:  Nov 1996   Jun Makino
//    version 2:  May 1999   Steve McMillan
//.............................................................................
//
//  Externally visible functions:
//
//	void print_sort_counters
//	void fast_get_nodes_to_move
//		- a single function which takes care of everything...
//	void dump_node_list
//	void clean_up_hdyn_schedule
//
//.............................................................................

#include "hdyn.h"

typedef struct {hdynptr b;
		xreal next_time;}	// Save next_time here because
node_time;				// get_next_time() can be EXTREMELY
					// slow -- avoiding it gains nearly
					// a factor of 5!  (Steve 5/99)

// Note that xreal comparisons likely take longer, so may need to
// optimize this... (Steve, 5/00)

//-----------------------------------------------------------------------------
//
// Static variables used to define and manage the scheduling:
//
//	static	node_time	*nodes			// list of node/time
//							// pairs for scheduling
//
//	static	int		work_size		// allocated size of
//							// the nodes array
//
//	static	int		nbody			// current number of
//							// nodes on the list
//
//	static	int		nprev			// number of nodes on
//							// previous "to_move"
//							// list
//
//	static	int		istack[NSTACK]		// used in quicksort_nr
//	static	node_time	arr_copy[MAX_COPY]	// used in merge_sort
//
//	static	int		counter[NC]		// used in statistics
//	static	int		total_count[NC]		// gathering
//
//-----------------------------------------------------------------------------

local inline void swap(node_time a[], int i, int j)
{
    // ASSUME assignment a = b is OK for structures.

    node_time tmp;

    tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
}

// Makino's quicksort -- recursive.

local void quicksort(node_time a[], int left, int right)
{
    // Sort nodes in array a, according to get_next_time().

    int i,j;
    xreal v;

    if (right > left) {
	bool is_sorted = true;
	for (i=left; i< right; i++) {
	    if (a[i].next_time > a[i+1].next_time) {
		is_sorted = false;
		i = right;
	    }
	}
	if (is_sorted) return;
	i = (left + right)/2;
	swap(a, i, right);
	v = a[right].next_time; i = left- 1; j = right;
	for (;;) {
	    while (a[++i].next_time < v && i < right) ;
	    while (a[--j].next_time > v && j > left) ;
	    if (i >= j) break;
	    swap(a, i, j);
	}
	swap(a, i, right);
	quicksort(a, left, i - 1);
	quicksort(a, i + 1, right);
    }
}

// Copy of Makino's quicksort.  Use of separate identical functions
// allows us to keep track of sorting for profiling purposes.

local void quicksort2(node_time a[], int left, int right)
{
    // Sort nodes in array a, according to get_next_time().

    int i,j;
    xreal v;

    if (right > left) {
	bool is_sorted = true;
	for (i=left; i< right; i++) {
	    if (a[i].next_time > a[i+1].next_time) {
		is_sorted = false;
		i = right;
	    }
	}
	if (is_sorted) return;
	i = (left + right)/2;
	swap(a, i, right);
	v = a[right].next_time; i = left- 1; j = right;
	for (;;) {
	    while (a[++i].next_time < v && i < right) ;
	    while (a[--j].next_time > v && j > left) ;
	    if (i >= j) break;
	    swap(a, i, j);
	}
	swap(a, i, right);
	quicksort2(a, left, i - 1);
	quicksort2(a, i + 1, right);
    }
}

// Numerical Recipes non-recursive quicksort.  VERY slow for nearly
// ordered data...  (Steve 5/99)

#define M	7
#define NSTACK	128

static int istack[NSTACK];

local void quicksort_nr(node_time arr[], int n)

// Note Numerical Recipes unit-offset arrays!  Array arr[] runs
// from 1 to n.  Handle in the calling sequence...

{
    int i, l = 1, ir = n, j, k;
    int jstack = 0;
    node_time a;

    for (;;) {
	if (ir-l < M) {

	    for (j = l+1; j <= ir; j++) {
		a = arr[j];
		for (i = j-1; i >= 1; i--) {
		    if (arr[i].next_time <= a.next_time) break;
		    arr[i+1] = arr[i];
		}
		arr[i+1] = a;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];

	} else {

	    k = (l+ir) >> 1;

	    swap(arr, k, l+1);

	    if (arr[l+1].next_time > arr[ir].next_time)
		swap(arr, l+1, ir);

	    if (arr[l].next_time > arr[ir].next_time)
		swap(arr, l, ir);

	    if (arr[l+1].next_time > arr[l].next_time)
		swap(arr, l+1, l);

	    i = l+1;
	    j = ir;
	    a = arr[l];
	    for (;;) {
		do i++; while (arr[i].next_time < a.next_time);
		do j--; while (arr[j].next_time > a.next_time);
		if (j < i) break;
		swap(arr, i, j);
	    }
	    arr[l] = arr[j];
	    arr[j] = a;
	    jstack += 2;

	    if (jstack > NSTACK)
		err_exit("NSTACK too small in sort_quicknr.");

	    if (ir-i+1 >= j-l) {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    } else {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
}

// Heapsort from Numerical Recipes.  Note Numerical Recipes unit-offset
// arrays!  Array arr[] runs from 1 to n...  Faster than NR quicksort,
// but still slower than Makino's version for our purposes.  (Steve 5/99)

local void heapsort(node_time arr[], int n)
{
    int i, ir, j, l;
    node_time a;

    if (n < 2) return;

    l = (n >> 1)+1;
    ir = n;

    for (;;) {
	if (l > 1) {
	    a = arr[--l];
	} else {
	    a = arr[ir];
	    arr[ir] = arr[1];
	    if (--ir == 1) {
		arr[1] = a;
		break;
	    }
	}
	i = l;
	j = l+l;
	while (j <= ir) {
	    if (j < ir && arr[j].next_time < arr[j+1].next_time)
		j++;
	    if (a.next_time < arr[j].next_time) {
		arr[i] = arr[j];
		i = j;
		j <<= 1;
	    } else j = ir+1;
	}
	arr[i] = a;
    }
}

// Insertion sort from Numerical Recipes, except that arrays start at 0.

inline void insertion_sort(node_time arr[], int n)	// sort from 0 to n-1
{
    int i, j;
    node_time a;

    for (j = 1; j < n; j++) {
	a = arr[j];
	i = j-1;
	while (i >= 0 && arr[i].next_time > a.next_time) {
	    arr[i+1] = arr[i];
	    i--;
	}
	arr[i+1] = a;
    }
}

// Quick and dirty merging of an array consisting of two ordered
// subarrays.  2-3 times faster than Makino's quicksort.  (Steve 5/99)

#define MAX_COPY	256
static node_time	arr_copy[MAX_COPY];

local void merge_sort(node_time arr[], int n, int np)
{
    // Ordered subarrays are 0 to np-1, np to n-1.  We expect np << n
    // in typical use, and arr[np-1] > arr[n-1] by construction..

    if (np > MAX_COPY) {
	quicksort2(arr, 0, n-1);
	return;
    }

    // Work with a copy of the first subarray (simplest to code).

    int i, j;

    for (i = 0; i < np; i++) arr_copy[i] = arr[i];

    // This could probably be accelerated using the block structure...

    i = 0, j = np;

    for (int k = 0; k < n; k++) {

	if (j >= n || arr_copy[i].next_time < arr[j].next_time)
	    arr[k] = arr_copy[i++];
	else
	    arr[k] = arr[j++];		// note that k < j always

    }
}

//------------------------------------------------------------------------

#define NC 12
static int counter[NC] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static int total_count[NC] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

local void update_sort_counters(int n)
{
    int nc = 1, ic = 0;
    while (nc < n) {
	nc *= 2;
	ic++;
    }
    if (ic >= NC) ic = NC - 1;
    counter[ic]++;
    total_count[ic]++;
}

void print_sort_counters()
{
    cerr << endl << "Sort counters:" << endl << "    ";
    for (int i = 0; i < NC; i++) cerr << counter[i] << " ";
    cerr << endl << "    ";
    for (int i = 0; i < NC; i++) cerr << total_count[i] << " ";
    cerr << endl << endl;
    for (int i = 0; i < NC; i++) counter[i] = 0;
}

local void print_sort_compare(node_time nodes[], int iprev, int icurr,
			      char* s)
{
    // iprev was the last element in the previously sorted list.
    // icurr  is the last element in the new sorted list.
    // The list has not yet been resorted.

    // Count how many elements of the "prev" list have moved into
    // the rest of the list.  On entry, the list is ordered from 0
    // to iprev, and from iprev+1 to icurr.

    int n_moved = 0;
    xreal tp = nodes[iprev+1].next_time;

    for (int i = iprev; i >= 0; i--) {
	if (nodes[i].next_time <= tp) break;
	n_moved++;
    }	

    fprintf(stderr, "quicksort: %20.12f %8d %8d %8d %s\n",
	    (real)nodes[0].b->get_system_time(), iprev+1, icurr+1, n_moved, s);
}

//------------------------------------------------------------------------

#define WHICH_SORT	4		// 4 is preferred (Steve 5/99)

#define N_ISORT		10

// Convenient to separate sorts by length for debugging/profiling...

inline local void QS(node_time a[], int n, int np = 0)
{
    if (n <= N_ISORT && (WHICH_SORT != 4 || np <= 0)) {
	insertion_sort(a, n);
	return;
    }

    if (WHICH_SORT == 1) {

	// Makino's quicksort:

	if (np < 0)
	    quicksort(a, 0, n-1);
	else
	    quicksort2(a, 0, n-1);

    } else if (WHICH_SORT == 2) {

	// NR quicksort:

	quicksort_nr(a-1, n);			// note NR offset!

    } else if (WHICH_SORT == 3) {

	// NR heapsort:

	heapsort(a-1, n);			// note NR offset!

    } else if (WHICH_SORT == 4) {

	// Steve's merge_sort:

	if (np > 0)

	    merge_sort(a, n, np);		// np is the start of the
						// second ordered subarray
	else {

	    if (np < 0)
		quicksort(a, 0, n-1);
	    else
		quicksort2(a, 0, n-1);

	}
    }
}

local void sort_node_array(node_time nodes[], int nprev, int nbody)
{

#if 0
    cerr << "before sort\n";
    for (int i = 0; i < nprev; i++) {
	PRC(i); nodes[i].b->print_label(cerr);
	PRL(nodes[i].next_time);
    }
#endif

    // Update times in the node_time structure.

    for (int i = 0; i < nprev; i++)
	nodes[i].next_time = nodes[i].b->get_next_time();

    // Sort the previous list (bodies last stepped forward).
    // (No need to sort if nprev = 1.)

    int p = cerr.precision(13);		// nonstandard precision

    if (nprev > 1) {

	// Sort nodes array from 0 to nprev-1.

	update_sort_counters(nprev);
	QS(nodes, nprev, -1);
    }

#if 0
    for (int i = 0; i < nprev; i++) {
	PRC(i); nodes[i].b->print_label(cerr);
	PRL(nodes[i].next_time);
    }
#endif

    // Perhaps some check is necessary for unperturbed motion,
    // whose timestep is integer multiple of the basic step...

    if (nprev < nbody) {

	// The list from 0 to nprev-1 is ordered, as is the list from
	// nprev to nbody-1.  See if the two parts need to be merged.

	// int p = cerr.precision(HIGH_PRECISION);
	// PRC(nodes[nprev-1].next_time); PRL(nodes[nprev].next_time);
	// cerr.precision(p);

	if (nodes[nprev-1].next_time > nodes[nprev].next_time) {

	    // Subarrays overlap.  Need to sort or merge.
	    // These searches could probably be accelerated using
	    // bisection (see NR, p. 117).

	    int ilow = nprev-1;
	    for (int i = nprev-2; i >= 0; i--) {
		if (nodes[i].next_time <= nodes[nprev].next_time) break;
		ilow = i;
	    }
	
	    int ihigh = nprev;
	    for (int i = nprev+1; i < nbody; i++) {
		if (nodes[i].next_time > nodes[nprev-1].next_time) break;
		ihigh = i;
	    }

	    // Out-of-order region runs from ilow to ihigh, inclusive.

#if 0
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << "sort_node_array: correct sequence from ";
	    cerr << ilow  << " ("; nodes[ilow].b->print_label(cerr);
	    cerr << ") to ";
	    cerr << ihigh << " ("; nodes[ihigh].b->print_label(cerr);
	    cerr << ")" << endl << "                 at time "
		 << nodes[nprev-1].next_time << "  nprev = " << nprev
		 << endl;

	    for (int i = 0; i < 10; i++)
	    	cerr << i << " " << nodes[i].b->format_label()
	    	     << " " << nodes[i].next_time << endl;

	    print_sort_compare(nodes, nprev-1, ihigh, "");
	    cerr.precision(p);
#endif

	    // Sort from ilow to ihigh, possibly using the additional
	    // information that the subarrays below and above nprev are
	    // already ordered.

	    update_sort_counters(ihigh-ilow+1);
	    QS(nodes+ilow, ihigh-ilow+1, nprev-ilow);

	}
    }
    cerr.precision(p);
}

// Static data:

static int work_size = 0;
static node_time * nodes = NULL;

static int nbody = 0;
static int nprev = 0;

// Allow possibility of cleaning up if necessary:

void clean_up_hdyn_schedule() {if (nodes) delete [] nodes;}

// ***** NOTE that the scheduling may become badly corrupted if bodies not
// ***** on the list have been integrated (e.g. by synchronize_tree()).
// ***** If this is done, then a reset *must* be forced.

void fast_get_nodes_to_move(hdyn * b,
			    hdyn * list[],
			    int &nlist,
			    xreal & tmin,
			    bool & reset)
{
    if (nbody == 0) reset = true;

    if (reset) {

	// Initialize the array of node pointers.

	int n = 0;
	for_all_nodes(hdyn, b, bb) n++;

	if (work_size < n) {
	    work_size = n*2;
	    if (nodes) delete [] nodes;
	    nodes = new node_time[work_size];
	}

	n = 0;
	for_all_nodes(hdyn, b, bbb) {
	    if (!bbb->is_root()
		&& (bbb->is_top_level_node()
		    || bbb->get_elder_sister() == NULL)) {
		nodes[n].b = bbb;
		nodes[n].next_time = bbb->get_next_time();
		n++ ;
	    }
	}
	nbody = nprev = n;
    }

    reset = false;

    bool debug = false;
    // debug = (b->get_system_time() > 100);

    sort_node_array(nodes, nprev, nbody);

    nlist = 0;
    tmin = nodes[0].next_time;
    for (int i = 0; i < nbody; i++) {
	if (nodes[i].next_time < tmin) {

	    // Fatal scheduling error...

	    cerr << endl << "Error in fast_get_nodes_to_move():" << endl;
	    cerr.precision(HIGH_PRECISION);

	    PRC(tmin), PRL((real)tmin);
	    PRC(i); cerr << "node " << nodes[i].b->format_label() << endl;
	    PRL(nodes[i].next_time);

	    for (int j = 0; j <= i; j++) {
		PRC(j);
		cerr << nodes[j].b->format_label() << "  "
		     << nodes[j].b->get_time() << "  "
		     << nodes[j].next_time << endl;
	    }

	    err_exit("internal error: sort failed\n");

	} else if (nodes[i].next_time > tmin) {

	    i = nbody;	// aka break!

	} else {

	    list[i] = nodes[i].b;
	    nlist = i + 1;
	}
    }
    nprev = nlist;

    if (debug) {
      cerr << endl;

      int p = cerr.precision(HIGH_PRECISION);
      PRC(nbody); PRC(nlist); PRC(tmin);
      cerr.precision(p);

      cerr << list[0]->format_label() << endl;
    }
}


void dump_node_list(int n) // default = 1000000000
{
    int p = cerr.precision(HIGH_PRECISION);
    cerr << endl << "Current node list (nbody = " << nbody << "):\n";
    for (int i = 0; i < nbody && i < n; i++) {
        cerr << i << "  " << nodes[i].b->format_label()
	     << "  " << nodes[i].next_time << endl;
    }
    cerr << endl << flush;
    cerr.precision(p);
}


void dump_node_list_for(char *s)
{
    int p = cerr.precision(HIGH_PRECISION);
    cerr << "Current node list entry for " << s << endl;
    for (int i = 0; i < nbody; i++) {
	if (streq(nodes[i].b->format_label(), s)) {
	    if (i > 0)
		cerr << i-1 << "  " << nodes[i-1].b->format_label()
		     << "  " << nodes[i-1].next_time << endl;
	    cerr << i << "  " << nodes[i].b->format_label()
		 << "  " << nodes[i].next_time << endl;
	    if (i < nbody-1)
		cerr << i+1 << "  " << nodes[i+1].b->format_label()
		     << "  " << nodes[i+1].next_time << endl;
	}
    }
    cerr << endl << flush;
    cerr.precision(p);
}
