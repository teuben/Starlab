
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// worldline:  Functions to manage and manipulate 4-D trees
////             Defines worldbundle, worldline, and segment
////             classes and member functions.
////
////             First argument is input data file.
////             No options.
//
//----------------------------------------------------------------------
//
// Worldbundle member functions:
//
//	worldbundle::worldbundle(tdyn *b)
//	void worldbundle::print()
//	void worldbundle::dump()
//	int worldbundle::find_index(real id)
//	int worldbundle::find_index(char *name)
//	int worldbundle::find_index(pdyn *b)
//	worldline *worldbundle::find_worldline(real id)
//	worldline *worldbundle::find_worldline(char *name)
//	worldline *worldbundle::find_worldline(pdyn *b)
//	void worldbundle::attach(tdyn *bn, int verbose)
//	void worldbundle::print_worldline(char *name, real dt)
//
// Worldline member function:
//
//	void worldline::dump(int offset)
//
// Segment member function:
//
//	void segment::print(char *label)
//
// Other externally visible functions:
//
//	real unique_id(char *name)
//	real unique_id(pdyn *b)
//	void print_id(void *id, char *label)
//	worldbundle *read_bundle(istream &s, int verbose)
//	tdyn *find_event(tdyn *bn, real t)
//	void print_event(tdyn *bn, real t)
//	vector interpolate_pos(tdyn *p, real t, tdyn *bn)
//	vector get_pos(tdyn *b, tdyn *bn, real t)
//	pdyn *create_interpolated_tree(worldbundle *wb, real t)
//
//	real physical_mass_scale()
//
//----------------------------------------------------------------------

#include "worldline.h"

#define NEW 0

#ifndef TOOLBOX

// Functions to define unique identification for particles.

#define FAC	(1.0/8192)		// power of 2...

// Overloaded function...

real unique_id(int index)
{
    return 1 + FAC*index;
}

real unique_id(char *name)
{
    // Associate a unique identifier with the specified name.
    // This ID will be used to track the particle throughout.
    // A leaf has 1 < id < 2, a binary CM has 2 < id < 3, etc.

    // First see if name can be interpreted as an integer.

    int index = atoi(name);
    if (index > 0) return unique_id(index);

    // This is presumably a center of mass node.  Turn its
    // name into an identifying number.

    // Look for commas and parentheses separating integers.
    // Combine ((a,b),c) --> a + FAC*(b + FAC*c), etc.
    // Factors are ~arbitrary.

    // Really need a better way of generating unique IDs, but
    // this should do for now.

    char tname[1024];
    strcpy(tname, name);	// work with a local copy

    // Parse the name.

    real id = 0, fac = 1;
    int num = 0, n = 0, l = strlen(tname);

    for (int i = 0; i < l; i++) {
	if (tname[i] >= '0' && tname[i] <= '9') {
	    if (num == 0) {
		fac *= FAC;
		num = i+1;
	    }
	} else {
	    if (num > 0) {
		tname[i] = '\0';
		id += fac*atoi(tname+num-1);
		n++;
		num = 0;
	    }
	}
    }

    id += n;
    return id;
}

real unique_id(pdyn *b)
{
    // Return a unique non-negative real number identifying this node.
    // Use the index if available (assume unique!), otherwise, construct
    // an integer from the node name.  If the node has no name or index,
    // ignore it (id = -1).

    // The root node should be called "root" or "?", and will therefore
    // be assigned an ID of 0.

    if (b->get_index() > 0)
	return unique_id(b->get_index());

    else if (b->get_name())
	return unique_id(b->format_label());

    else
        return -1;
}

void print_id(void *id, char *s = NULL, int offset = 0)
{
    real *r = (real *)id;
    unsigned long long *u = (unsigned long long *)id;

    real rr = *r;
    unsigned long long uu = *u;

    if (offset > 0) PRI(offset);
    if (s)
	cerr << s << ": ";
    else
	cerr << "    ";

    fprintf(stderr, "%.16lf = %qu = %qx\n", rr, uu, uu);
}

void worldline::dump(int offset)	// default = 0
{
    // Print a worldline, start to finish.

    segment *s = get_first_segment();
    int is = 0;
    while (s) {
	real t_start = s->get_t_start();
	real t_end = s->get_t_end();
	PRI(offset+4); cerr << "segment " << is << ": ";
		PRC(t_start); PRL(t_end);
	print_id(&id, "id", offset+8);
	tdyn *b = s->get_first_event();
	int ie = 0;
	while (b) {
	    real t = b->get_time();
	    PRI(offset+8);
	    cerr << "event " << ie << ": name = " << b->format_label() << ", ";
	    PRL(t); 
	    b = b->get_next();
	    ie++;
	}
	if (ie == 1) {
	    PRI(offset+8); cerr << "*** single event ***" << endl;
	}
	s = s->get_next();
	is++;
    }
}

void segment::print(char *label)
{
    if (label == NULL) label = "    ";
    cerr << label << "name " << first_event->format_label()
	 << ", id = " << id
	 << ", t_start = " << t_start
	 << ", t_end = " << t_end
	 << endl;
}

//======================================================================

local int wlcompare(const void *a, const void *b)	// a and b are actually
							// of type worldlineptr*
{
    worldlineptr wa = *((worldlineptr*)a);		// ugly!
    worldlineptr wb = *((worldlineptr*)b);

    int result = 0;
//    cerr << "comparing " << wa->get_id() << " and " << wb->get_id() << endl;

    if (wa->get_id() < wb->get_id())
	result = -1;
    else if (wa->get_id() > wb->get_id())
	result = 1;

//    PRL(result);
    return result;
}

worldbundle::worldbundle(tdyn *b)
{
    // Create a new worldbundle from the nodes below b.

    nw_max = 0;
    for_all_nodes(tdyn, b, bb) nw_max++;

    nw_max *= 2;				// conservative

    bundle = new worldlineptr[nw_max];
    nw = 0;

    t_min = VERY_LARGE_NUMBER;
    t_max = -VERY_LARGE_NUMBER;

    // Add nodes to the list.

    for_all_nodes(tdyn, b, bb) {

	// *** Flag recomputation of acc and jerk. ***

	bb->clear_jerk();

	real t = bb->get_time();
	worldline *w = new worldline(bb);

	bundle[nw++] = w;

	t_min = min(t_min, t);
	t_max = max(t_max, t);
    }

    // Sort the list by ID.

    qsort((void*)bundle, (size_t)nw, sizeof(worldlineptr), wlcompare);

//    for (int i = 0; i < nw; i++) {
//	PRC(i); PRC(bundle[i]->get_id());
//	cerr << bundle[i]->get_first_segment()
//			 ->get_first_event()->format_label() << endl;
//    }
}

void worldbundle::print()
{
    // Print out basic information about this worldbundle.

    for (int i = 0; i < nw; i++) {

	worldline *w = bundle[i];
	segment *s = w->get_first_segment();
	int segnum = 0;

	cerr << i << " (" << s->get_first_event()->format_label();
	int p = cerr.precision(HIGH_PRECISION);
	cerr << ", id = " << w->get_id() << ")";
	cerr.precision(p);
	cerr << "  t_start = " << w->get_t_start()
	     << ", t_end = " << w->get_t_end()
	     << endl;

	while (s) {

	    segnum++;
	    int nn = 1;

	    tdyn *b = s->get_first_event();
	    while (b->get_next()) {
		b = b->get_next();
		nn++;
	    }

	    cerr << "        segment " << segnum << ":  "
		 << s->get_t_start() << " (" << nn << ") " << s->get_t_end()
		 << endl;

	    s = s->get_next();
	}
    }
}

void worldbundle::dump(int offset)	// default = 0
{
    for (int i = 0; i < get_nw(); i++) {
	PRI(offset); cerr << "worldline " << i << ":" << endl;
	get_worldline(i)->dump(offset);
    }
}
//======================================================================

int worldbundle::find_index(real id)
{
    // Find the index of the worldline corresponding to id.

    // Return values:	0 - nw-1	index of id
    //			-1		id out of range below
    //			-nw-1		id out of range above
    //			-nw -- -2	-(first index above id) - 1
    //
    // (so -find_index - 1 is where a new worldline should go in the list).

    // The list is ordered in id; use bisection to search it.

    int loc;

    if (nw <= 0 || id < bundle[0]->get_id())
	loc = -1;

    else if (id > bundle[nw-1]->get_id())
	loc = -nw-1;

    else {

	int low = 0, high = nw-1;
	loc = (low+high)/2;

	if (bundle[low]->get_id() == id)
	    loc = low;

	else if (bundle[high]->get_id() == id)
	    loc = high;

	else if (bundle[loc]->get_id() != id) {

	    bool found = false;

	    while (low < high-1) {

		if (bundle[loc]->get_id() < id)
		    low = loc;
		else if (bundle[loc]->get_id() > id)
		    high = loc;
		else {
		    found = true;
		    break;
		}

		loc = (low+high)/2;
	    }

	    if (!found) loc = -high - 1;
	}
    }

    return loc;
}

int worldbundle::find_index(char *name)	{return find_index(unique_id(name));}
int worldbundle::find_index(pdyn *b)	{return find_index(unique_id(b));}

worldline *worldbundle::find_worldline(real id)
{
    int i = find_index(id);
    if (i >= 0 && i < nw)
	return bundle[i];
    else
	return NULL;
}

worldline *worldbundle::find_worldline(char *name)
{return find_worldline(unique_id(name));}

worldline *worldbundle::find_worldline(pdyn *b)
{return find_worldline(unique_id(b));}

//======================================================================

// Update the world bundle.

void worldbundle::attach(tdyn *bn,
			 int verbose)		// default = 0
{
    // Locate and attach all components of bn in the 4tree hierarchy.

    if (bn) {

	for_all_nodes(tdyn, bn, b) {

	    // *** Flag recomputation of acc and jerk. ***

	    b->clear_jerk();

	    // Find the start of b's current worldline segment, or
	    // add b to the list.  There are three possibilities:
	    //
	    //	(1) b's ID is on the list and b is valid
	    //				==> extend the current segment
	    //				    or begin a new one if the
	    //				    previous b was defunct
	    //
	    //	(2) b's ID is on the list and b is defunct
	    //				==> end the current segment
	    //
	    //	(3) b's ID is not on the list
	    //				==> extend the list and open
	    //				    a new segment
	    //
	    // Cases (1) and (2) are handled at the same time.

	    // Find b on the worldline list or add it to the list
	    // (maintaining the list ordering).

	    real id = unique_id(b);

	    if (id >= 0) {

		real t = b->get_time();
		int loc = find_index(id);

		if (loc >= 0) {

		    // Extend an existing worldline.

		    worldline *w = bundle[loc];
		    segment *s = w->get_last_segment();	// "last" = "current"

		    // Check case (2) for previous instance of b.

		    if (s->get_last_event()->is_defunct()) {

			if (verbose)
			    cerr << "starting new segment for "
				 << b->format_label()
				 << " (id = " << id
				 << ") at t = " << t
				 << endl;

			// Create a new segment starting at b.

			s = new segment(b);
			w->add_segment(s);

		    } else {

			// Case (1): extend the current segment.

			s->add_event(b, true);	  // "true" ==> skip id check

			// Note from Steve (8/01):  For unknown reasons, the id
			// check in add_event may cause errors when optimization
			// is turned on.  Problems arise in multiple systems (so
			// we are working near the last significant digit), where
			// the ids may return unequal even though the bits are
			// identical.  Didn't occur in the July 29, 2001 version
			// of the code.  Does occur in the August 2 version!
			//
			// Workaround for now: skip the check, which is redundant
			// here anyway.  Longer-term: improve the handling of
			// unique_id.

		    }

		    w->set_t_end(t);

		} else {

		    // Case (3): create a new worldline starting at b.
		    // Maintain the ordering of the array by placing the
		    // new pointer at location -loc - 1.

		    if (nw >= nw_max) {

			// Extend the storage array -- use new and
			// delete to mimic realloc().

			nw_max *= 5;
			nw_max /= 4;

			cerr << "extending bundle array for worldbundle "
			     << this << ":  nw_max = "
			     << nw_max << endl;

			worldlineptr *tmp = new worldlineptr[nw_max];
			for (int i = 0; i < nw; i++)
			    tmp[i] = bundle[i];
			delete [] bundle;
			bundle = tmp;
		    }

		    // Where to insert the new worldline:

		    loc = -loc - 1;

		    for (int i = nw-1; i >= loc; i--)
			bundle[i+1] = bundle[i];

		    bundle[loc] = new worldline(b);
		    nw++;

		    if (verbose)
			cerr << "adding " << b->format_label()
			     << " (id = " << id
			     << ") at t = " << b->get_time()
			     << endl;
		}

		t_max = max(t_max, t);
	    }
	}
    }
}

//======================================================================

static real physical_mass = 0;

real mass_scale_factor()
{
    return physical_mass;
}

// From Steve (10/00):
//
// Root dumps are done by kira twice per synchronization interval, so each
// worldbundle is self-contained (this is the *same* as the way in which
// tree changes are handled...).

worldbundle *read_bundle(istream &s,
			 int verbose)			// default = 0
{
    // Read a bundle of worldlines, from root dump to root dump.

    // First read the base root node.  Assume without checking that
    // the next portion of input data starts with a complete system.

    tdyn *b = get_tdyn(s, NULL, NULL, false);		// "stripped tdyn"
    if (!b || !streq(b->format_label(), "root")) return NULL;

    // Need to keep track of the physical initial mass, but this only
    // appears in the first Log story and is not saved in the usual
    // story structure, to save space.  Grab it here and make it
    // globally accessible.

    if (b->get_log_story())
	physical_mass = getrq(b->get_log_story(), "physical_initial_mass");

    // Create the initial list of node IDs.

    worldbundle *wb = new worldbundle(b);

    if (verbose)
	cerr << endl
	     << "new worldbundle: created initial list, nw = "
	     << wb->get_nw() << endl;

    // Now read in individual nodes and attach them to the 4tree.
    // Stop when another tree is read in. 

    while (b = get_tdyn(s, NULL, NULL, false)) {	// "stripped tdyn"

	// A new tree fragment will be read in every time the tree
	// changes, so pointers should be correctly maintained.
	// Links are found at the start of each worldline segment.

	// Connect b to the 4tree.

	wb->attach(b, verbose);

	// Stop once we have read in another complete system.

	if (b->get_oldest_daughter()
	    && streq(b->format_label(), "root")) {
	    if (verbose)
		cerr << "break at t = " << b->get_time() << endl;
	    break;
	}
    }

    if (verbose) {

	cerr << "after last input:  nw = " << wb->get_nw()
	     << endl;

	// Print some global diagnostics:

	for (int i = 0; i < wb->get_nw(); i++) {
	    worldline *w = wb->get_bundle()[i];

	    // First check segments.

	    int seg = 0, nseg = 0;
	    segment *s = w->get_first_segment();
	    tdyn *b = s->get_first_event();
	    real t = s->get_t_start();

	    while (s) {
		nseg++;
		s = s->get_next();
	    }

	    s = w->get_first_segment();
	    int n_jump = 0;

	    while (s) {

		if (s->get_t_start() >= s->get_t_end()) {
		    int p = cerr.precision(HIGH_PRECISION);
		    cerr << "worldline " << i << " (" << b->format_label()
			 << "), id = " << w->get_id() << endl;
		    cerr.precision(p);
		    PRI(10); cerr << "zero-length segment at t = "
				  << s->get_t_start()
				  << ";  segment " << seg << " of " << nseg
				  << endl;
		}

		if (s->get_t_start() > t) {
		    n_jump++;
		    if (verbose > 1) {
			int p = cerr.precision(HIGH_PRECISION);
			cerr << "worldline " << i << " (" << b->format_label()
			    << "), id = " << w->get_id() << endl;
			cerr.precision(p);
			PRI(10); cerr << "jump from " << t
				      << " to " << s->get_t_start()
				      << " after segment " << seg
				      << " of " << nseg
				      << endl;
		    }
		}

		t = s->get_t_end();
		seg++;
		s = s->get_next();
	    }

	    if (verbose == 1 && n_jump > 0) {
		cerr << "worldline " << i << " (" << b->format_label()
		     << "), id = " << w->get_id()
		     << " has " << n_jump << " jump";
		if (n_jump > 1) cerr << "s";
		cerr << endl;
	    }

	    // Then check events within each segment.

	    seg = 0;
	    s = w->get_first_segment();
	    while (s) {

		b = s->get_first_event();

		if (b->get_time() != s->get_t_start()) {
		    int p = cerr.precision(HIGH_PRECISION);
		    cerr << "worldline " << i << " (" << b->format_label()
			 << "), id = " << w->get_id() << endl;
		    cerr.precision(p);
		    PRI(10); cerr << "t_start =  " << s->get_t_start()
				  << ", t = "
				  << b->get_time() << " in segment " << seg
				  << " of " << nseg << endl;
		}

		while (b && b->get_next()) b = b->get_next();

		if (b->get_time() != s->get_t_end()) {
		    int p = cerr.precision(HIGH_PRECISION);
		    cerr << "worldline " << i << " (" << b->format_label()
			 << "), id = " << w->get_id() << endl;
		    cerr.precision(p);
		    PRI(10); cerr << "t_end =  " << s->get_t_end() << ", t = "
				  << b->get_time() << " in segment " << seg
				  << " of " << nseg << endl;
		}

		seg++;
		s = s->get_next();
	    }
	}
    }

    return wb;
}

int count_segments(worldbundle *wb)
{
    int ns = 0;
    for (int i = 0; i < wb->get_nw(); i++) {
	worldline *w = wb->get_worldline(i);
	segment *s = w->get_first_segment();
	while (s) {
	    ns++;
	    s = s->get_next();
	}
    }
    return ns;
}

int count_events(worldbundle *wb)
{
    int ne = 0;
    for (int i = 0; i < wb->get_nw(); i++) {
	worldline *w = wb->get_worldline(i);
	segment *s = w->get_first_segment();
	while (s) {
	    tdyn *b = s->get_first_event();
	    while (b) {
		ne++;
		b = b->get_next();
	    }
	    s = s->get_next();
	}
    }
    return ne;
}

//======================================================================

// Navigation of the 4tree data and diagnostic output.

tdyn *find_event(tdyn *bn, real t)
{
    // Find time t along the worldline segment starting at bn.
    // Return a pointer to the event immediately preceding t.

    tdyn *b = bn;

    while (b && b->get_next() && b->get_next()->get_time() < t)
	b = b->get_next();

    // The portion of the trajectory spanning t runs from b to b->next.

    return b;
}

void print_event(tdyn *bn, real t)
{
    // Print info on the portion of the worldline portion
    // starting at bn that spans t.

    tdyn *b = find_event(bn, t);

    if (b) {
	PRC(t); PRC(b->get_time()); PRL(b->get_next());
	if (b->get_next()) PRL(b->get_next()->get_time());
    }
}

//----------------------------------------------------------------------
#include "print_local.C"	// avoid repeating local print functions
//----------------------------------------------------------------------

// Interpolation of part of the 4tree to a specific time.

vector interpolate_pos(tdyn *p, real t,
		       tdyn *bn)		// default = NULL; specifies
						// actual base node
{
    // The range (p to p->next) includes time t.
    // Interpolate and return the pos of this particle.

    // Check...

    if (p->get_time() > t) {
	cerr << "interpolate_pos: error 1: specified time below range: ";
	PRC(p->format_label()); PRC(t); PRL(p->get_time());
	// print_details(wb, p, t); // wb not available here...

	return p->get_pos();
    }

    // Special case:

    if (p->get_time() == t) return p->get_pos();

    tdyn *n = p->get_next();

    if (!n) {
	cerr << "interpolate_pos: error 2: next node not found: ";
	PRC(p->format_label()); PRL(t);
	PRI(26); PRC(bn); PRL(bn->get_time());
	PRI(26); PRC(p); PRL(p->get_time());
	PRI(26); PRL(p->get_time()-t);
	if (bn) {
	    PRI(26); PRL(bn->format_label());
	}
	// print_details(wb, p, t);

	return p->get_pos();
    }

    if (n->get_time() < t) {
	cerr << "interpolate_pos: error 3: specified time above range: ";
	PRC(p->format_label()); PRC(t); PRL(p->get_time());
	// print_details(wb, p, t);

	return p->get_pos();
    }

    // Time t is included in the range.

    // Interpolate using pos and vel for now...
    // Note that we overwrite acc and jerk by equivalent
    // vectors that guarantee continuity of pos and vel.

    real tp = p->get_time();
    real dt = n->get_time() - tp;

    if (dt > 0) {

	// Recompute acc/2 and jerk/6 to fit pos and vel.

	// *** Flag this to prevent recalculation by setting ***
	// *** jerk = 0 as the tree is constructed.          ***

	if (p->get_jerk()[0] == 0) {
	    real dti = 1/dt;

	    p->set_acc((3*(n->get_pos()-p->get_pos())
			- dt*(2*p->get_vel()+n->get_vel()))*dti*dti);
	    p->set_jerk((p->get_vel()+n->get_vel()
			 - 2*dti*(n->get_pos()-p->get_pos()))*dti*dti);
	}

	dt = t - tp;

	return p->get_pos() + dt * (p->get_vel()
				    + dt * (p->get_acc()
					    + dt * p->get_jerk()));
    } else

	return p->get_pos();
}

vector interpolate_vel(tdyn *p, real t,
		       tdyn *bn)		// default = NULL; specifies
						// actual base node
{
    // Interpolate and return the vel of this particle.

    // Only called after a call to interpolate_pos, so skip checks
    // and recomputation of acc and jerk.

    // Special case:

    if (p->get_time() == t) return p->get_vel();

    tdyn *n = p->get_next();
    real tp = p->get_time();
    real dt = n->get_time() - tp;

    if (dt > 0) {

	dt = t - tp;
	return p->get_vel() + dt * (2*p->get_acc()
				    + dt * 3*p->get_jerk());
    } else

	return p->get_vel();
}

// Old versions of interpolate_pos() and interpolate_vel() are
// (for now) saved in old_interpolate.C

vector get_pos(tdyn *b, tdyn *bn,
	       real t)			// default = -VERY_LARGE_NUMBER
{
    // Return the current absolute position of node b, properly
    // corrected for its location in the tree.

    // The node of interest is b; the base node for the worldline
    // segment (containing all parental information) is bn.

    if (t == -VERY_LARGE_NUMBER) t = b->get_time();
    vector pos = interpolate_pos(b, t, bn);

    // The tree structure associated with the base node should extend
    // to the top level, but will not contain the root node unless bn
    // is part of a full dump.

    // Ascend the base tree, including parent contributions (at time t)
    // as needed.  Exclude root, but note that is_root() will fail, as
    // it simply checks for a null parent.  For now, check the name
    // explicitly...

    tdyn *bp = bn->get_parent();

    while (bp) {

	// Find the portion of the parent worldline spanning t.

	tdyn *p = find_event(bp, t);

	pos += interpolate_pos(p, t, bn);
	bp = bp->get_parent();
    }

    return pos;
}

void worldbundle::print_worldline(char *name,
				  real dt)	// default = 0
{
    // Print out the entire worldline (all segments) of particle 'name'.
    // Take steps of length dt (0 ==> just use the stored times).

    int loc = find_index(name);			// find the worldline
						// for particle 'name'
    if (loc >= 0) {

	worldline *w = bundle[loc];
	segment *s = w->get_first_segment();

	if (dt == 0) {

	    while (s) {				    // loop over segments

		tdyn *bn = s->get_first_event();    // starting point
		tdyn *b = bn;
		while (b) {			    // loop over events
		    cout << "    " << b->get_time()
			 << " " << get_pos(b, bn) << endl << flush;
		    b = b->get_next();
		}

		s = s->get_next();
	    }

	} else {

	    real t = get_t_min();
	    while (t < get_t_max() - 0.5*dt) {

		// Code here assumes that we are moving sequentially
		// through the data, as we start with the previous s.

		while (s && s->get_t_end() < t) s = s->get_next();

		if (s) {
		    tdyn *bn = s->get_first_event();
		    tdyn *b = find_event(bn, t);
		    cout << "    " << t
			 << " " << get_pos(b, bn) << endl << flush;
		}

		t += dt;
	    }
	}

    } else
	cerr << name << " not found" << endl;
}

//======================================================================

// Creation of an entire interpolated tree at some specific time.

//----------------------------------------------------------------------
#include "print_local2.C"	// avoid repeating local print functions
//----------------------------------------------------------------------

//----------------------------------------------------------------------
#include "attach_new_node.C"	// avoid repeating local print functions
//----------------------------------------------------------------------

//----------------------------------------------------------------------
#include "update_node.C"	// avoid repeating local functions
//----------------------------------------------------------------------

local inline void add_to_interpolated_tree(worldbundle *wb,
					   worldline *w, segment *s,
					   pdyn *root, real t, bool vel,
					   bool debug)
{
    // Current worldline is w, segment is s.  Find the base node
    // (containing all relevant tree information) for this particle,
    // and the current event.

    // *** By construction, on entry, w and s represent a leaf. ***

    tdyn *bn = s->get_first_event();

    if (debug) {
	cerr << "s = " << s << endl
	     << "first event = " << bn << " "
	     << bn->format_label() << " "
	     << bn->get_time() << endl;

	print_details(wb, bn, t);
    }

    // Create the entire portion of the tree corresponding to the
    // top-level node of bn.  Note that top will probably be bn most
    // of the time.

    // Need to be careful with top-level nodes, as they  generally
    // won't have parent nodes (no root), so standard functions like
    // "is_top_level_node()" and "get_top_level_node()" may fail.
    // (Choice of root in calling function should correct this...)

    tdyn *top = bn;
    while (top->get_parent()
	   && unique_id(top->get_parent()) > 0)
	top = top->get_parent();

    if (debug) {
	PRC(top); PRL(top->format_label());
	// put_node(cerr, *top, false);
    }

    // Copy the entire tree below top to the new tree.

    for_all_nodes(tdyn, top, bb) {

	// Find the worldline of bb.

	worldline *ww;

	if (bb == bn)
	    ww = w;
	else {
	    ww = wb->find_worldline(bb);
	    if (!ww)
		cerr << "create_interpolated_tree: error:"
		     << "can't find worldline of subtree member." << endl;
	}

	if (ww) {

	    // Create and attach a new node.

	    attach_new_node(wb, ww, root, top, bb, debug);

	    // Update the new tree entry: copy or interpolate all
	    // relevant quantities from bb to curr.

	    update_node(wb, ww, t, bb, top, vel, debug);
	}
    }
}

#define EPS 1.e-12

pdyn *create_interpolated_tree(worldbundle *wb, real t,
			       bool vel)		// default = false
{
    // All leaves on the bundle list are still current, by construction.
    // Find the base segment corresponding to each and use it to build
    // a new pdyn tree interpolated to the current time.

    bool debug = false;

    // Try to take care of rounding error in t.  Management of worldbundles
    // should be the responsibility of the calling program.  (This would not
    // be an issue if timesteps were constrained to be powers of 2...)

    real dt = t - wb->get_t_max();
    if (dt > EPS)
	return NULL;
    else if (dt > 0)
	t = wb->get_t_max();

    dt = t - wb->get_t_min();
    if (dt < -EPS)
	return NULL;
    else if (dt < 0)
	t = wb->get_t_min();

    // Find the array of worldlines (worldbundle), and initialize
    // the tree_node pointers.

    worldlineptr *bundle = wb->get_bundle();

    for (int i = 0; i < wb->get_nw(); i++)
	bundle[i]->clear_tree_node();

    // Create a root node.

    pdyn *root = new pdyn(NULL, NULL, false);
    root->set_system_time(t);
    root->set_pos(0);

    // Establish root as the global root node (should clean up problems
    // with "top_level_node" functions...).  Stored in location 0 of the
    // worldline array.

    root->set_root(root);
    bundle[0]->set_tree_node(root);

    // Loop through remaining worldlines and take action for leaves only.
    // Logic: we are really dealing with top-level nodes, but we don't know
    // which component we will encounter first.  When we find a component
    // of a clump of particles, we update the entire clump, and flag all
    // components accordingly.

    for (int i = 1; i < wb->get_nw(); i++) {

	worldline *w = bundle[i];

	if (!w->get_tree_node()) {

	    // Tree entry for this node has not yet been created/updated.

	    real id = w->get_id();

	    if (id >= 1 && id < 2) {

		// Worldline w represents a leaf.  Locate its current segment.

		if (debug) {
		    cerr.precision(16);
		    cerr << endl
			 << "following worldline " << i << ", id = "
			 << w->get_id() << endl;
		}

		segment *s = w->get_first_segment();
		while (s && s->get_t_end() < t) s = s->get_next();

		if (debug) {
		    PRC(w); cerr << "segment "; PRC(s);
		    print_events(s, t);
		}

		// If s is NULL, it should mean that w refers to a CM node
		// that does not exist at time t.  Particle worldlines
		// should be continuous -- could check in debugging mode...

		if (s) add_to_interpolated_tree(wb, w, s,
						root, t, vel, debug);
	    }
	}
    }

    return root;
}

#else

main(int argc, char *argv[])
{
    // Fixed format on the command line:  filename  time.

    if (argc <= 2) exit(1);

    ifstream s(argv[1]);
    if (!s) exit(2);

    real t = atof(argv[2]);

    worldbundle *wb = read_bundle(s);

    pdyn *root = create_interpolated_tree(wb, t, true);
    put_node(cout, *root, false);
    rmtree(root);

#if 0
    // Some statistics:

    wb->print();

    cerr << wb->get_nw() << " worldlines, "
	 << count_segments(wb) << " segments, "
	 << count_events(wb) << " events, t = "
	 << wb->get_t_min() << " to " << wb->get_t_max()
	 << endl;
#endif

#if 0
    // Locate a specific time along the worldline of some particle.

    char name[128];
    real t;
    while (1) {
	cerr << endl << "time, name: "; cin >> t >> name;
	if (cin.eof()) {
	    cerr << endl;
	    break;
	}
	int loc = wb->find_index(name);
	if (loc >= 0) {
	    PRC(unique_id(name)); PRL(loc);
	    worldline *w = wb->get_bundle()[loc];
	    segment *s = w->get_first_segment();
	    while (s && s->get_t_start() > t) s = s->get_next();
	    print_event(s->get_first_event(), t);
	} else
	    cerr << "not found" << endl;
    }
#endif

#if 0
    // Print out the entire worldline of a specified particle,
    // taking steps of the specified length.

    real dt;
    char name[128];
    while (1) {

	cerr << endl << "dt, name: " << flush; cin >> dt >> name;
	if (cin.eof()) {
	    cerr << endl;
	    break;
	}

	wb->print_worldline(name, dt);
    }
#endif

#if 0
    real t = 0.8;
    while (t < 0.85) {
	pdyn *root = create_interpolated_tree(wb, t);
	put_node(cout, *root, false);
	rmtree(root);
	t += 0.01;
    }
#endif

}

#endif
