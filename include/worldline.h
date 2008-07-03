
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file worldline.h   Classes and definitions for 4-D trees.
//
//  version 1:  Oct 2000   Steve McMillan
//  version 2:
//
//  This file includes:
//  1) definition of classes segment, worldline, and worldbundle

#ifndef  STARLAB_WORLD_H
#  define  STARLAB_WORLD_H

#include "tdyn.h"

// Define id data type:

#define NEW_UNIQUE_ID_T			// remove to restore old version

#ifndef NEW_UNIQUE_ID_T
    typedef real unique_id_t;		// old
#else
    typedef int unique_id_t;		// new
#endif

unique_id_t unique_id(int index);
unique_id_t unique_id(char *name);
unique_id_t unique_id(node *b);

// Class definitions:

// A worldline is an indexed pointer to the start of a linked list of
// worldline segments.  Each segment consists of a series of events
// (tdyns) along a particle trajectory.  Tree changes result in new
// worldline segments for all particles involved.  The full worldline
// is the entirety of all such segments.
//
// There is presently considerable redundancy in the data stored.
// We will refine the data structures as the package evolves...

///\a segment: A series of events (tdyns) along a particle trajectory.

class segment {

    private:

	unique_id_t id;			///< global identifier
	tdyn *first_event;		///< first event on the list
	tdyn *last_event;		///< last event on the list
	real t_start;			///< start time of this segment
	real t_end;			///< end time of this segment
	segment *next;			///< pointer to the next segment

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

	unique_id_t get_id()		{return id;}
	tdyn *get_first_event()		{return first_event;}
	tdyn *get_last_event()		{return last_event;}
	real get_t_start()		{return t_start;}
	real get_t_end()		{return t_end;}

	segment *get_next()		{return next;}
	void set_next(segment *s)	{next = s;}

	void add_event(tdyn *b, bool accept = false) {

	    if (accept				// id check is often redundant
		|| unique_id(b) == id) {	// and may cause problems if
						// optimization is turned on
		if (last_event)
		    last_event->set_next(b);
		b->set_prev(last_event);

		last_event = b;
		t_end = b->get_time();
	    }
	}

	void print(const char *label = NULL);
};

/// \a worldline: Indexed pointer to the start of a linked list of segments.

/// A worldline is an indexed pointer to the start of a linked list of
/// worldline segments.  Each segment consists of a series of events
/// (tdyns) along a particle trajectory.  Tree changes result in new
/// worldline segments for all particles involved.  The full worldline
/// is the entirety of all such segments.

class worldline {

    private:

	unique_id_t id;			///< global identifier
	segment *first_segment;		///< first segment
	segment *last_segment;		///< last segment
	real t_start;			///< start time
	real t_end;			///< end time

	int start_esc_flag;		///< initial escaper flag; set by kira
	int end_esc_flag;		///< final escaper flag; set by kira
	real t_esc;			///< time when flag changed (to come)

	// Management of tree traversal:

	real t_curr;			///< current time
	tdyn *current_event;		///< current event
	segment *current_segment;	///< current segment

	pdyn *tree_node;		///< pointer to the corresponding node
					// in the interpolated tree at time t

    public:

	worldline() {			// empty worldline
	    id = -1;
	    first_segment = last_segment = NULL;
	    t_start = t_end = 0;

	    t_curr = 0;
	    current_event = NULL;
	    current_segment = NULL;

	    tree_node = NULL;

	    start_esc_flag = end_esc_flag = 0;
	    t_esc = VERY_LARGE_NUMBER;
	}

	worldline(segment *s) {		// worldline of a single segment
	    id = s->get_id();
	    first_segment = last_segment = s;
	    t_start = s->get_t_start();
	    t_end = s->get_t_end();

	    t_curr = 0;
	    current_event = NULL;
	    current_segment = NULL;

	    tree_node = NULL;

	    start_esc_flag = end_esc_flag = 0;
	    t_esc = VERY_LARGE_NUMBER;
	}

	worldline(tdyn *b) {		// worldline of a single event
	    segment *s = new segment(b);
	    id = s->get_id();
	    first_segment = last_segment = s;
	    t_start = s->get_t_start();
	    t_end = s->get_t_end();

	    t_curr = -VERY_LARGE_NUMBER;
	    current_event = NULL;
	    current_segment = NULL;

	    tree_node = NULL;

	    // Note:  We use the prev and next pointers in tdyn to carry
	    // temporary information from scan_dyn_story() about cluster
	    // membership, to avoid adding more data to the class.  At
	    // this stage, all prev and next pointers should be NULL.
	    // Restore those settings here and transfer the membership
	    // information to the worldline, where it belongs.
	    //
	    // NOTE:  The default (no esc_flag set via prev) is esc_flag = 0.

	    start_esc_flag = end_esc_flag = 0;
	    t_esc = VERY_LARGE_NUMBER;

	    if (b->get_prev()) {		// escape flag is set
		b->set_prev(NULL);
		start_esc_flag = end_esc_flag = 1;
		t_esc = b->get_time();		// should be overwritten
	    }

	    if (b->get_next()) {		// t_esc is specified
		real *x = (real*) b->get_next();
		b->set_next(NULL);
		t_esc = *x;
		delete x;
	    }
	}

	unique_id_t get_id()		{return id;}
	segment *get_first_segment()	{return first_segment;}
	segment *get_last_segment()	{return last_segment;}
	real get_t_start()		{return t_start;}

	real get_t_end()		{return t_end;}
	void set_t_end(real t)		{t_end = t;}

	inline void set_start_esc_flag(int f)	{start_esc_flag = f;}
	inline int get_start_esc_flag()	{return start_esc_flag;}
	inline bool is_member(real t)	{
	    if (start_esc_flag == 1)
		return false;
	    else if (end_esc_flag == 0)
		return true;
	    else
		return (t < t_esc);
	}

	inline void set_end_esc_flag(int f)	{end_esc_flag = f;}
	inline int get_end_esc_flag()	{return end_esc_flag;}

	real get_t_esc()		{return t_esc;}
	void set_t_esc(real t)		{t_esc = t;}

	real get_t_curr()		{return t_curr;}
	void set_t_curr(real t)		{t_curr = t;}

	tdyn *get_current_event()	{return current_event;}
	void set_current_event(tdyn* b)	{current_event = b;}

	segment *get_current_segment()	{return current_segment;}
	void set_current_segment(segment* s)	{current_segment = s;}

	pdyn *get_tree_node()		{return tree_node;}
	void set_tree_node(pdyn* b)	{tree_node = b;}
	void clear_tree_node()		{tree_node = NULL;}

	void add_segment(segment *s, bool accept = false) {
	    if (accept				// id check is often redundant
		|| s->get_id() == id) {		// and may cause problems if
						// optimization is turned on
		last_segment->set_next(s);
		last_segment = s;
		t_end = s->get_t_end();
	    }
	}

	void dump(int offset = 0, real t1 = 0, real t2 = VERY_LARGE_NUMBER);
	void check(int i = -1);
};

typedef worldline *worldlineptr;	// convenient...

/// \a worldbundle:  A group of worldlines representing an entire N-body system.
// Also includes structures for data management purposes.

class worldbundle {

    private:

	worldlineptr *bundle;		///< array of worldline pointers
	int	     nw;		///< length of the array
	int	     nw_max;		///< maximum length of the array
	real	     t_min;		///< minimum time
	real	     t_max;		///< maximum time
	real	     t_int;		///< time of most recent interpolation
	pdyn	     *root;		///< root node of interpolated tree

    public:

	worldbundle() {
	    bundle = NULL;
	    nw = nw_max = 0;
	    t_min = VERY_LARGE_NUMBER;
	    t_max = -VERY_LARGE_NUMBER;
	    t_int = -VERY_LARGE_NUMBER;
	    root = NULL;
	}

	worldbundle(tdyn *b);

	worldlineptr *get_bundle()	{return bundle;}
	worldline *get_worldline(int i)	{return bundle[i];}

	int get_nw()			{return nw;}
	real get_t_min()		{return t_min;}
	real get_t_max()		{return t_max;}

	void set_t_int(real t)		{t_int = t;}
	real get_t_int()		{return t_int;}

	void set_pdyn_root(pdyn *p)	{root = p;}
	pdyn *get_pdyn_root()		{return root;}

	int find_index(unique_id_t id);
	int find_index(char *name);
	int find_index(_pdyn_ *b);

	worldline *find_worldline(unique_id_t id);
	worldline *find_worldline(char *name);
	worldline *find_worldline(_pdyn_ *b);

	void attach(tdyn *b, int verbose = 0);

	void print();
	void print_worldline(char *name, real dt = 0);
	void dump(int offset = 0);
	void check();
	int get_n_daughters();
};

typedef worldbundle *worldbundleptr;

// Globally visible functions:

/// Print unique ID corresponding to id.

void print_id(void *id, const char *label = NULL);

/// Read a bundle of worldlines describing an N-body system for a time interval.

worldbundle *read_bundle(istream &s, int verbose = 0);

/// Read worldbundles describing an N-body system over an extended time.

void read_bundles(istream &s, worldbundleptr wh[], int& nh,
		  int verbose = 0);

/// Count the number of segments in a worldbundle.

int count_segments(worldbundle *wb);

/// Count the number of events in a worldbundle.

int count_events(worldbundle *wb);

/// Find time t along the segment of worldline w starting at base node bn.

tdyn *find_event(worldline *w, tdyn *bn, real t);

/// Print info on the portion of the worldline starting at bn that spans time t.

void print_event(worldline *w, tdyn *bn, real t);

// Old:

/// Return the interpolated position at time t in the segment starting at p.

vec interpolate_pos(tdyn *p, real t, tdyn *bn = NULL);

// New:

/// Return the interpolated position at time t in the segment starting at p.

void interpolate_pos(tdyn *p, real t, vec& pos, bool inc, tdyn *bn);

/// Set pos to the interpolated value.

void set_interpolated_pos(tdyn *p, real t, vec& pos,
			  tdyn *bn = NULL);
/// Set curr->pos to the interpolated value.

void set_interpolated_pos(tdyn *p, real t, pdyn *curr,
			  tdyn *bn = NULL);

/// Increment pos by the interpolated value.

void inc_interpolated_pos(tdyn *p, real t, vec& pos,
			  tdyn *bn = NULL);

/// Return the interpolated velocity at time t in the segment starting at p.

vec interpolate_vel(tdyn *p, real t, tdyn *bn = NULL);

/// Return the physical mass scale.

real mass_scale_factor();

/// Create an interpolated tree at time t and return a pointer to the root node.

pdyn *create_interpolated_tree(worldbundle *wb, real t, bool vel = false);

/// Create an interpolated tree at time t and return a pointer to the root node.

pdyn *create_interpolated_tree2(worldbundle *wb, real t, bool vel = false);

/// Force memory allocation for each bundle by creating an interpolated tree.

void preload_pdyn(worldbundleptr wh[], int nh,
		  bool verbose = false, bool vel = false);

/// Set all root nodes to use the specified center for center tracking.

const char *set_center(worldbundleptr wh[], int nh, int center_number,
		       bool verbose = false);

/// Identify which type of center we are using for center tracking.

int get_n_center();

/// Return the center we are currently using for center tracking (1 or 2).

int get_center();

/// Current possibilities are (1) "standard-center" and (2) "bound-center".

const char *get_center_id(int center_number = -1);

/// Position of the current center.

vec get_center_pos();

/// Velocity of the current center.

vec get_center_vel();

/// Is p a member of the specified worldbundle?

bool is_member(worldbundle *wb, pdyn *p);

/// Nuber of particles in the clump containing id.

int id_n_clump(unique_id_t id);

#endif
