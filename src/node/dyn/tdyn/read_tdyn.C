
#include "tdyn.h"

// Functions to define unique identification for particles.

#define FAC	(1.0/8192)		// power of 2...

// Overloaded function...

local inline real unique_id(int index)
{
    return 1 + FAC*index;
}

local inline real unique_id(char *name)
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

local inline real unique_id(tdyn *b)
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

//======================================================================

// A worldline is an indexed pointer to the start of a linked list of
// worldline segments.  Each segment consists of a series of events
// (tdyns) along a particle trajectory.  Tree changes result in new
// worldline segments for all particles involved.  The full worldline
// is the entirety of all such segments.
//
// There is presently considerable redundancy in the data stored.
// We will refine the data structures as the package evolves...

class segment {

    private:

	real id;			// global identifier
	tdyn *first_event;		// first event on the list
	tdyn *last_event;		// last event on the list
	real t_start;			// start time of this segment
	real t_end;			// end time of this segment
	segment *next;			// pointer to the next segment

    public:

	segment() {			// empty segment
	    id = -1;
	    first_event = last_event = NULL;
	    t_start = t_end = 0;
	    next = NULL;
	}

	segment(tdyn *b) {		// segment containing a single event
	    id = unique_id(b);
	    first_event = last_event = b;
	    t_start = t_end = b->get_time();
	    next = NULL;
	}

	real get_id()			{return id;}
	tdyn *get_first_event()		{return first_event;}
	tdyn *get_last_event()		{return last_event;}
	real get_t_start()		{return t_start;}
	real get_t_end()		{return t_end;}

	segment *get_next()		{return next;}
	void set_next(segment *s)	{next = s;}

	void add_event(tdyn *b) {
	    if (unique_id(b) == id) {
		last_event->set_next(b);
		b->set_prev(last_event);
		last_event = b;
		t_end = b->get_time();
	    }
	}

	void print(char *label = NULL);
};

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

class worldline {

    private:

	real id;			// global identifier
	segment *first_segment;		// first segment
	segment *last_segment;		// last segment
	real t_start;			// start time
	real t_end;			// end time

	// Management of tree traversal:

	real t_curr;			// current time
	tdyn *current;			// current node
	tdyn *current_base;		// current base node
	tdyn *tree_node;		// pointer to the corresponding node
					// the interpolated tree at time t

    public:

	worldline() {			// empty worldline
	    id = -1;
	    first_segment = last_segment = NULL;
	    t_start = t_end = 0;

	    t_curr = 0;
	    current = current_base = tree_node = NULL;
	}

	worldline(segment *s) {		// worldline of a single segment
	    id = s->get_id();
	    first_segment = last_segment = s;
	    t_start = s->get_t_start();
	    t_end = s->get_t_end();

	    t_curr = 0;
	    current = current_base = tree_node = NULL;
	}

	worldline(tdyn *b) {		// worldline of a single event
	    segment *s = new segment(b);
	    id = s->get_id();
	    first_segment = last_segment = s;
	    t_start = s->get_t_start();
	    t_end = s->get_t_end();

	    t_curr = 0;
	    current = current_base = tree_node = NULL;
	}

	real get_id()			{return id;}
	segment *get_first_segment()	{return first_segment;}
	segment *get_last_segment()	{return last_segment;}
	real get_t_start()		{return t_start;}

	real get_t_end()		{return t_end;}
	void set_t_end(real t)		{t_end = t;}

	real get_t_curr()		{return t_curr;}
	void set_t_curr(real t)		{t_curr = t;}

	tdyn *get_current()		{return current;}
	void set_current(tdyn* b)	{current = b;}

	tdyn *get_current_base()	{return current_base;}
	void set_current_base(tdyn* b)	{current_base = b;}

	tdyn *get_tree_node()		{return tree_node;}
	void set_tree_node(tdyn* b)	{tree_node = b;}
	void clear_tree_node()		{tree_node = NULL;}

	void add_segment(segment *s) {
	    if (s->get_id() == id) {
		last_segment->set_next(s);
		last_segment = s;
		t_end = s->get_t_end();
	    }
	}
};

typedef worldline *worldlineptr;	// convenient...

local int wlcompare(const void *a, const void *b)	// a and b are actually
							// of type worldlineptr
{
    if (((worldlineptr)a)->get_id() < ((worldlineptr)b)->get_id())
	return -1;
    else if (((worldlineptr)a)->get_id() < ((worldlineptr)b)->get_id())
	return 1;
    else
	return 0;
}

//======================================================================

// A worldbundle is a group of worldlines, including data structures
// for management purposes.

class worldbundle {

    private:

	worldlineptr *bundle;		// array of worldline pointers
	int	     nw;		// length of the array
	int	     nw_max;		// maximum length of the array
	real	     t_min;		// minimum time
	real	     t_max;		// maximum time

    public:

	worldbundle() {
	    bundle = NULL;
	    nw = nw_max = 0;
	    t_min = VERY_LARGE_NUMBER;
	    t_max = -VERY_LARGE_NUMBER;
	}

	worldbundle(tdyn *b);

	worldlineptr *get_bundle()	{return bundle;}

	int get_nw()			{return nw;}
	real get_t_min()		{return t_min;}
	real get_t_max()		{return t_max;}

	int find_id(real id);
	int find_id(char *name);
	int find_id(tdyn *b);

	worldline *find_worldline(real id);
	worldline *find_worldline(char *name);
	worldline *find_worldline(tdyn *b);

	void attach(tdyn *b);

	void print();
	void print_worldline(char *name, real dt = 0);
};

worldbundle::worldbundle(tdyn *b)
{
    // Create a new worldbundle from the nodes below b.

    nw_max = 0;
    for_all_nodes(tdyn, b, bb) nw_max++;

    nw_max *= 2;				// conservative

    bundle = new worldlineptr[nw_max];
    nw = 0;

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

//======================================================================

int worldbundle::find_id(real id)
{
    // Find the index of the worldline corresponding to id.

    // Return values:	0 - nw-1	index of id
    //			-1		id out of range below
    //			-nw-1		id out of range above
    //			-nw -- -2	-(first index above id) - 1
    //
    // (so -find_id - 1 is where a new worldline should go in the list).

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

int worldbundle::find_id(char *name)	{return find_id(unique_id(name));}
int worldbundle::find_id(tdyn *b)	{return find_id(unique_id(b));}

worldline *worldbundle::find_worldline(real id)
{
    int i = find_id(id);
    if (i >= 0 && i < nw)
	return bundle[i];
    else
	return NULL;
}

worldline *worldbundle::find_worldline(char *name)
{return find_worldline(unique_id(name));}

worldline *worldbundle::find_worldline(tdyn *b)
{return find_worldline(unique_id(b));}

//======================================================================

// Update the world bundle.

void worldbundle::attach(tdyn *bn)
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
		int loc = find_id(id);

		if (loc >= 0) {

		    // Extend an existing worldline.

		    worldline *w = bundle[loc];
		    segment *s = w->get_last_segment();	// "last" = "current"

		    // Check case (2) for previous instance of b.

		    if (s->get_last_event()->is_defunct()) {

			cerr << "starting new segment for "
			     << b->format_label()
			     << " (id = " << id
			     << ") at t = " << t
			     << endl;

			// Create a new segment starting at b.

			s = new segment(b);
			w->add_segment(s);

		    } else

			// Case (1): extend the current segment.

			s->add_event(b);

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

			cerr << "extending bundle array:  nw_max = "
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

local worldbundle *read_bundle(istream &s)
{
    // Read a bundle of worldlines, from root dump to root dump.

    // First read the base root node.  Assume without checking that
    // the input data start with a complete system.

    tdyn *b = get_tdyn(s);
    if (!b || !streq(b->format_label(), "root")) return NULL;

    // Create the initial list of node IDs.

    worldbundle *wb = new worldbundle(b);
    cerr << "created initial list:  nw = " << wb->get_nw() << endl;

    // Now read in individual nodes and attach them to the 4tree.
    // Stop when another tree is read in. 

    while (b = get_tdyn(s)) {

	// A new tree fragment will be read in every time the tree
	// changes, so pointers should be correctly maintained.
	// Links are found at the start of each worldline segment.

	// Connect b to the 4tree.

	wb->attach(b);

	// Stop once we have read in another complete system.

	if (b->get_oldest_daughter()
	    && streq(b->format_label(), "root")) {
	    cerr << "break at t = " << b->get_time() << endl;
	    break;
	}
    }

    cerr << "after last input:  nw = " << wb->get_nw() << endl;
    return wb;
}

//======================================================================

// Navigation of the 4tree data.

local tdyn *find_event(tdyn *bn, real t)
{
    // Find time t along the worldline segment starting at bn.
    // Return a pointer to the event immediately preceeding t.

    tdyn *b = bn;

    while (b && b->get_next() && b->get_next()->get_time() < t)
	b = b->get_next();

    // The portion of the trajectory spanning t runs from b to b->next.

    return b;
}

local void print_event(tdyn *bn, real t)
{
    // Print info on the portion of the worldline portion
    // starting at bn that spans t.

    tdyn *b = find_event(bn, t);

    if (b) {
	PRC(t); PRC(b->get_time()); PRL(b->get_next());
	if (b->get_next()) PRL(b->get_next()->get_time());
    }
}

local vector interpolate_pos(tdyn *p, real t)
{
    // The range (p to p->next) includes time t.
    // Interpolate and return pos.

    // Check...

    if (p->get_time() > t) {
	cerr << "interpolate_pos: error 1: ";
	PRC(p->format_label()); PRC(t); PRL(p->get_time());
	return vector(0);
    }

    // Special case:

    if (p->get_time() == t) return p->get_pos();

    tdyn *n = p->get_next();

    if (!n) {
	cerr << "interpolate_pos: error 2: ";
	PRC(p->format_label()); PRC(t); PRL(p->get_time());
	return vector(0);
    }

    if (n->get_time() < t) {
	cerr << "interpolate_pos: error 3: ";
	PRC(p->format_label()); PRC(t); PRL(p->get_time());
	return vector(0);
    }

    // Time t is included in the range.

    // Interpolate using pos and vel for now...
    // Note that we overwrite acc and jerk by equivalent
    // vectors that guarantee continuity of pos and vel.

    real tp = p->get_time();
    real dt = n->get_time() - tp;

    if (dt > 0) {

	real dti = 1/dt;

	// Recompute acc/2 and jerk/6 to fit pos and vel.

	// *** Flag this to prevent recalculation by setting ***
	// *** jerk = 0 as the tree is constructed.          ***

	if (p->get_jerk()[0] == 0) {
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

local vector get_pos(tdyn *b, tdyn *bn, real t = -VERY_LARGE_NUMBER)
{
    // Return the current position of node b, properly corrected
    // for its location in the tree.

    // The node of interest is b; the base node for the worldline
    // segment (containing all parental information) is bn.

    if (t == -VERY_LARGE_NUMBER) t = b->get_time();
    vector pos = interpolate_pos(b, t);

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

	pos += interpolate_pos(p, t);
	bp = bp->get_parent();
    }

    return pos;
}

void worldbundle::print_worldline(char *name,
				  real dt)	// default = 0
{
    // Print out the entire worldline (all segments) of particle 'name'.
    // Take steps of length dt (0 ==> just use the stored times).

    int loc = find_id(name);			// find the worldline
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

local tdyn *create_interpolated_tree(worldbundle *wb, real t)
{
    // All leaves on the bundle list are still current, by construction.
    // Find the base segment corresponding to each and use it to build
    // a new tree interpolated to the current time.

    worldlineptr *bundle = wb->get_bundle();

    for (int i = 0; i < wb->get_nw(); i++)
	bundle[i]->clear_tree_node();

    tdyn *root = new tdyn(NULL, NULL, false);	// stripped version
    root->set_time(t);
    root->set_system_time(t);
    root->set_pos(0);

    bundle[0]->set_tree_node(root);

    for (int i = 0; i < wb->get_nw(); i++) {

	worldline *w = bundle[i];

	if (!w->get_tree_node()) {

	    // Tree entry for this node has not yet been created.

	    real id = w->get_id();

	    if (id >= 1 && id < 2) {

		// Worldline w represents a leaf.  Locate its current segment.

		segment *s = w->get_first_segment();
		while (s && s->get_t_end() < t) s = s->get_next();

		if (s) {

		    // Current segment is s.  Find the base node (containing
		    // all relevant tree information) and the current event.

		    tdyn *bn = s->get_first_event();

		    // Create the entire portion of the tree corresponding
		    // to the top-level node of bn.  Note that top will
		    // likely be bn most of the time.

		    // Need to be careful with top-level nodes, as they
		    // generally won't have parent nodes (no root), so
		    // standard funtions like "is_top_level_node()" and
		    // "get_top_level_node()" will likely fail.

		    tdyn *top = bn;
		    while (top->get_parent()
			   && unique_id(top->get_parent()) > 0)
			top = top->get_parent();

		    // Copy the entire tree below top to the new tree.

		    for_all_nodes(tdyn, top, bb) {

			// Find the worldline of bb.

			worldline *ww;

			if (bb == bn)
			    ww = w;
			else {
			    ww = wb->find_worldline(bb);
			    if (!ww) cerr << "error 1" << endl;
			}

			if (ww) {

			    tdyn *curr = new tdyn(NULL, NULL, false);

			    // The traversal of the tree is such that
			    // parents are always seen before children
			    // and elder sisters are always seen before
			    // younger sisters.

			    // Don't use add_node, as it will create a
			    // tree with nodes in the reverse order!

			    if (bb == top) {

				// Add bb to the end of the top-level list.

				tdyn *n = root->get_oldest_daughter();

				if (!n)
				    add_node(*curr, *root);

				else {
				    while (n->get_younger_sister())
					n = n->get_younger_sister();
				    n->set_younger_sister(curr);
				    curr->set_elder_sister(n);
				    curr->set_parent(root);
				}

			    } else {

				// Parent and elder sister pointers
				// are derived (awkwardly!) from the
				// existing tree structures.

				tdyn *par = wb->find_worldline(bb->get_parent())
				    	      ->get_tree_node();
				curr->set_parent(par);

				tdyn * bb_sis = bb->get_elder_sister();

				if (!bb_sis)
				    par->set_oldest_daughter(curr);
				else {
				    tdyn *sis = wb->find_worldline(bb_sis)
						  ->get_tree_node();
				    curr->set_elder_sister(sis);
				    sis->set_younger_sister(curr);
				}
			    }

			    ww->set_tree_node(curr);

			    // Update the new tree entry.

			    tdyn *b = find_event(bb, t);

			    if (b->get_name())
				curr->set_name(b->get_name());
			    else if (b->get_index() >= 0)
				curr->set_index(b->get_index());

			    curr->set_mass(b->get_mass());
			    curr->set_time(t);
			    curr->set_pos(interpolate_pos(b, t));

			    // cerr << "added " << curr->format_label() << endl;

			}
		    }
		}
	    }
	}
    }

    return root;
}


main(int argc, char *argv[])
{
    if (argc <= 1) exit(1);
    ifstream s(argv[1]);
    if (!s) exit(2);

    worldbundle *wb = read_bundle(s);

    // Some statistics:

    cerr << "final system:" << endl;
    wb->print();

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
	int loc = find_id(wb, name);
	if (loc >= 0) {
	    PRC(unique_id(name)); PRL(loc);
	    worldline *w = wb->bundle[loc];
	    segment *s = w->first_segment;
	    while (s && s->t_start > t) s = s->next;
	    print_event(s->first_event, t);
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

#if 1

    real t = 0.8;
    while (t < 0.85) {
	tdyn *root = create_interpolated_tree(wb, t);
	put_node(cout, *root);
	rmtree(root);
	t += 0.01;
    }

#endif

}
