
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// worldline2:  New functions to manage and manipulate 4-D trees
////
////              First argument is input data file.
////              No options.
//
//----------------------------------------------------------------------
//
// Other externally visible functions:
//
//	pdyn *create_interpolated_tree2(worldbundle *wb, real t)
//	char *set_center(worldbundleptr wh[], int nh,
//			 int center_number, bool verbose)
//	int get_center()
//	char *get_center_id(int center_number)
//	vector get_center_pos()
//	vector get_center_vel()
//
//----------------------------------------------------------------------

#include "worldline.h"
#include "inline.h"

#define NEW 1

#ifndef TOOLBOX

//----------------------------------------------------------------------
#include "print_local.C"	// avoid repeating local print functions
//----------------------------------------------------------------------

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

#if 0
local void dealloc_tree(pdyn *b,
			bool delete_b = true)
{
    // This is rmtree, except that it uses our own memory management,
    // if present (dealloc_pdyn).

    pdyn* d = b->get_oldest_daughter();
    while (d) {
	pdyn* tmp = d->get_younger_sister();
	dealloc_tree(d);
	d = tmp;
    }
    if (delete_b) dealloc_pdyn(b);  // optionally leave node itself untouched
}
#endif

local INLINE void clean_up_subtree(worldbundle *wb, pdyn *old, bool debug)
{
    if (old) {
	old = old->get_top_level_node();	    // should be OK, as
	    					    // pdyn root is defined

	if (debug)
	    cerr << "cleanup below " << old->format_label() << endl;

	// Zero out all tree_node() pointers, since it is likely
	// that not all will be reset during this function call.

	int ndel = 0;

	for_all_nodes(pdyn, old, o) {

	    if (debug)
		cerr << "deleting " << o->format_label() << endl;
	    ndel++;

	    int wi = o->get_worldline_index();

	    if (wi <= 0)
		cerr << "clean_up_subtree: error: "
		     << "deleting node with no worldline index." << endl;
	    else {
		wb->get_bundle()[wi]->clear_tree_node();
		wb->get_bundle()[wi]->set_t_curr(-VERY_LARGE_NUMBER);
	    }
	}		
	    
	detach_node_from_general_tree(*old);

	// Function rmtree() actually deletes the old nodes.  
	// May not be necessary to go that far.

	rmtree(old);

	// Replace by dealloc_tree() if we do our own memory management.

	if (debug)
	    cerr << "deleted " << ndel << " nodes" << endl;
    }
}

local INLINE void update_interpolated_tree(worldbundle *wb,
					   worldline *w, segment *s,
					   pdyn *root, real t, real t_int,
					   bool vel, bool debug)
{
    // Current worldline is w, segment is s, and w has not yet been
    // updated, according to t_curr.

    // *** By construction, on entry, w and s represent a leaf. ***

    // If times t and t_int lie in different segments, we must
    // rebuild the portion of the tree containing this particle.

    // Note that the "=" here are necessary.  If t_int was at the end
    // of the previous power-of-two segment, or at the start of the next,
    // we must still rebuild.

    // PROBLEM: This will cause an unnecessary rebuild if t_int happens
    // to be at the start or end of the proper segment (e.g. the start
    // or end of the worldbundle time range), which can slow the code
    // significantly.
    //
    // FIX: Save and compare the previous s rather than t_int.

    // bool rebuild = (t_int <= s->get_t_start() || t_int >= s->get_t_end());

    // Probably only need the last part of this test...

    bool rebuild = (t_int < s->get_t_start() || t_int > s->get_t_end()
		    || w->get_current_segment() != s);

    // Rebuild/update the current subtree.  Start by finding the base
    // node (containing all relevant tree information) for this particle,
    // and the current event.

    // Need to be careful with top-level tdyn nodes, as they
    // generally won't have parent nodes (no root), so standard
    // functions like is_top_level_node() and get_top_level_node()
    // may fail.

    tdyn *bn = s->get_first_event();
    tdyn *top = bn;

    while (top->get_parent()
	   && unique_id(top->get_parent()) > 0)
	top = top->get_parent();

    if (debug) {
	cerr << "s = " << s << endl
	     << "first event = " << bn << " "
	     << bn->format_label() << " "
	     << bn->get_time() << endl;

	print_details(wb, bn, t);

	PRC(top); PRL(top->format_label());
	// put_node(cerr, *top, false);
    }

    // Create the portion of the tree corresponding to the top-level
    // node of bn.  Note that top will probably be bn most of the time.

    // Copy the entire tree below top into the new tree.

    for_all_nodes(tdyn, top, bb) {

	// Find the worldline of bb.

	worldline *ww;

	if (bb == bn)
	    ww = w;
	else {
	    ww = wb->find_worldline(bb);
	    if (!ww)
		cerr << "update_interpolated_tree: error:"
		     << "can't find worldline of subtree member."
		     << endl;
	}

	if (ww) {

	    // Create and attach a new node, if necessary.

	    if (rebuild) {

		// Delete the old subtree.  The tree_node() pointer
		// still points into the old tree if not NULL, so use
		// this to locate and delete the out-of-date structure.
		// Note that the old "subtree" may consist of two
		// separate pieces (e.g. in a binary-binary interaction).

		// For efficiency, binary components are now handled
		// together, rather than separately.  Attach the younger
		// component of a binary too when we create the elder
		// one, because update_node will expect it to exist.

		if (bb == top || !bb->get_elder_sister()) {

		    // Function clean_up_subtree will check whether or
		    // not there is any work to be done.

		    clean_up_subtree(wb, ww->get_tree_node(), debug);
		    attach_new_node(wb, ww, root, top, bb, debug);

		    if (debug)
			cerr << "attached " << bb->format_label()
			     << " at time " << t << endl;

		    if (bb != top) {

			// Clean up and attach the binary sister too.

			tdyn *bbsis = bb->get_younger_sister();
			worldline *wwsis = wb->find_worldline(bbsis);

			clean_up_subtree(wb, wwsis->get_tree_node(), debug);
			attach_new_node(wb, wwsis, root, top, bbsis, debug);

			if (debug)
			    cerr << "attached " << bbsis->format_label()
				 << " at time " << t << endl;
		    }
		}
	    }

	    // Update the tree entry: copy or interpolate all relevant
	    // quantities from bb to curr.

	    if (ww->get_tree_node())
		update_node(wb, ww, t, bb, top, vel, debug);
	    else {
		cerr << "update_interpolated_tree: no tree_node"
		     << " for " << bb->format_label() << " at time " << t
		     << " t_int = " << t_int
		     << endl;
		PRL(rebuild);
	    }
	}
    }

    w->set_current_segment(s);		// note that only w has segment set;
					// other subtree members are unchanged
}

// For center tracking.  The position and velocity of the current
// center are computed in create_interpolated_tree2.  To change the
// center type, just swap in the new root pos and vels from the stories
// and set jerk = 0.  The rest will be handled automatically by the
// interpolation routines.

static int n_center = 2;
static int which_center = 0;	// 0 to n_center-1 are legal
static int which_std = 0;	// 1 or 2 is legal

static char *center_id[] = {"standard-center", "bound-center"};
static char *std_center_id[] = {"density-center", "modified-com"};

int get_n_center() {return n_center;}
int get_center() {return which_center;}

char *get_center_id(int center_number)		// default = -1
{
    // Don't attempt to tie the descriptive string to the strings used
    // in the root dyn stories.

    if (center_number < 0)
	return get_center_id(which_center);	// print the current center
    else if (center_number >= n_center)
	return NULL;
    else if (center_number == 0) {
	if (which_std == 1)
	    return std_center_id[0];
	else if (which_std == 2)
	    return std_center_id[1];
	else
	    return center_id[center_number];
    } else
	return center_id[center_number];
}

local bool scan_root_nodes(worldbundleptr wh[],
			   int nh,
			   bool verbose,
			   char *pos_id,
			   char *vel_id,
			   bool check)
{
    for (int ih = 0; ih < nh; ih++) {

	worldbundle *wb = wh[ih];
	worldline *w = wb->get_bundle()[0];	// root worldline

	// Each root worldline should contain only one segment, made
	// up of two events, but do this more generally and check.

	segment *s = w->get_first_segment();
	int is = 0;

	while (s) {

	    tdyn *b = s->get_first_event();
	    int ie = 0;

	    while (b) {

		if (check) {

		    // Look for pos and vel entries in the dyn story.

		    if (!find_qmatch(b->get_dyn_story(), pos_id)
			|| !find_qmatch(b->get_dyn_story(), vel_id))
			return false;
		} else {

		    // Story entries exist.  Set the root pos and vels.
		    // (Actually, only pos is currently needed, but...)

		    b->set_pos(getvq(b->get_dyn_story(), pos_id));
		    b->set_vel(getvq(b->get_dyn_story(), vel_id));
		    b->clear_acc();
		    b->clear_jerk();

		    if (which_center != 0)
			which_std = 0;
		    else if (which_std == 0)
			which_std = getiq(b->get_dyn_story(), "center_type");
		}

		b = b->get_next();
		ie++;
	    }

	    if (verbose & ie == 1)
		cerr << "warning: set_next_center: segment "
		     << is << "is a single event" << endl;

	    s = s->get_next();
	    is++;
	}

	if (verbose & is != 1)
	    cerr << "warning: set_next_center: root worldline  "
		<< ih << "contains more than one segment" << endl;
    }

    return true;
}

char *set_center(worldbundleptr wh[],	// entire worldbundle array
		 int nh,
		 int new_center,
		 bool verbose)		// default = false
{
    // Set all root nodes to use the specified center, and return
    // a string describing that center.

    if (new_center < 0 || new_center >= n_center) {
	if (verbose)
	    cerr << "set_next_center: invalid center number" << endl;
	return NULL;
    }

    int old_center = which_center;
    which_center = new_center;

    char center_pos_id[128], center_vel_id[128];
    if (which_center == 0) {
	strcpy(center_pos_id, "center_pos");	// default in read_tdyn
	strcpy(center_vel_id, "center_vel");
    } else if (which_center == 1) {
	strcpy(center_pos_id, "bound_center_pos");
	strcpy(center_vel_id, "bound_center_vel");
    }

    // Only make the change if the specified id exists in all root entries.
    // Scan the entire worldbundle array; also make some other basic checks.

    bool change = scan_root_nodes(wh, nh, verbose,
				  center_pos_id, center_vel_id, true);
    if (change) {
	scan_root_nodes(wh, nh, false, center_pos_id, center_vel_id, false);
	if (verbose)
	    cerr << "set_next_center: new center is "
		 << get_center_id(which_center) << endl;
    } else {
	if (verbose)
	    cerr << "set_next_center: unable to change center to "
		 << get_center_id(which_center) << endl;
	which_center = old_center;
    }

    return get_center_id(which_center);
}

static vector center_pos = 0;
vector get_center_pos()	{return center_pos;}

static vector center_vel = 0;
vector get_center_vel()	{return center_vel;}

// Membership determination:

bool is_member(worldbundle *wb, pdyn *p)
{
    // A node is a member if any of its children is.
    // Use knowledge of t_esc to determine the precise time of escape.

    for_all_leaves(pdyn, p, pp) {
	worldline *w = wb->find_worldline(pp);
	if (w && w->is_member(p->get_system_time())) return true;
    }
    return false;
}

local void interpolate_tree(worldbundle *wb, real t, real t_int,
			    pdynptr& root, bool vel, bool debug)
{
    // Compute the new tree for worldbundle wb at time t, with root
    // as the root node.

    // Find the array of worldlines (worldbundle), and initialize
    // the tree_node pointers.

    worldlineptr *bundle = wb->get_bundle();

    // All leaves on the bundle list are still current, by construction.
    // Find the base segment corresponding to each and use it to update
    // the pdyn tree interpolated to the current time.

    // Create a root node, if necessary.

    if (!root) {
	root = new pdyn(NULL, NULL, false);
	// if (!root) root = alloc_pdyn(NULL, NULL, false, true, debug);
	root->set_name("root");
	root->set_worldline_index(wb->find_index(root));
    }

    root->set_system_time(t);

    root->set_pos(0);			// unnecessary, but not wrong...
    root->set_vel(0);

    // Establish root as the global root node (should clean up problems
    // with "top_level_node" functions...).  Stored in location 0 of the
    // worldline array.

    root->set_root(root);
    bundle[0]->set_tree_node(root);

    // bundle[0]->set_t_curr(t);	// will be set in the loop below

    // Loop through remaining worldlines and take action for leaves only.
    // Logic: we are really dealing with top-level nodes, but we don't know
    // which component we will encounter first.  When we find a component
    // of a clump of particles, we update the entire clump, and flag all
    // components accordingly.

    // Handle the interpolation of the root node here too...
    // Root acc and jerk have been set up on input, so just interpolate
    // them as worldline 0 in the following loop.

    for (int i = 0; i < wb->get_nw(); i++) {

	worldline *w = bundle[i];

	if (w->get_t_curr() != t) {

	    // Worldline needs to be updated.

	    real id = w->get_id();

	    if (i == 0 || (id >= 1 && id < 2)) {

		// Worldline w represents a leaf.  Locate its current segment.

		if (debug) {
		    cerr.precision(16);
		    cerr << endl
			 << "updating worldline " << i << ", id = "
			 << w->get_id() << endl;
		}

		// If possible, start at the previous segment visited in
		// this worldline (may speed things up in normal use).
		// -- hard to see much improvememt in speed...

		segment *s;
#if 1
		s = w->get_current_segment();
		if (!s || s->get_t_start() > t)
#endif
		    s = w->get_first_segment();

		while (s && s->get_t_end() < t) s = s->get_next();

		if (debug) {
		    PRC(w); cerr << "segment "; PRC(s);
		    print_events(s, t);
		}

		// If s is NULL, it should mean that w refers to a CM node
		// that does not exist at time t.  Particle worldlines
		// should be continuous -- check in debugging mode...

		if (s) {
		    if (i == 0) {

			// Root node.

			tdyn *b = s->get_first_event();
			update_node(wb, w, t, b, b, vel, debug);

			center_pos = root->get_pos();
			center_vel = root->get_vel();

			// PRC(b->get_pos());
			// PRL(center_pos);

		    } else

			// Update entire clump of nodes, if necessary.

			update_interpolated_tree(wb, w, s,
						 root, t, t_int,
						 vel, debug);
		} else if (debug)
		    cerr << "NULL segment for leaf..." << endl;
	    }
	}
    }
}

#define EPS  1.e-12
#define EPS1 1.e-12				// kludge -- should be 0

local bool trim(worldbundle *wb, real& t)
{
    // Try to take care of rounding error in t.  Management of worldbundles
    // should be the responsibility of the calling program.  (This would not
    // be an issue if timesteps were constrained to be powers of 2...)

    // TEMPORARY BUG WORKAROUND: don't allow t to be exactly t_min or t_max!!!

    real dt = t - wb->get_t_max();
    if (dt > EPS)
	return false;
    else if (dt > -EPS1)
	t = wb->get_t_max() - EPS1;

    dt = t - wb->get_t_min();
    if (dt < -EPS)
	return false;
    else if (dt < EPS1)
	t = wb->get_t_min() + EPS1;

    return true;
}

pdyn *create_interpolated_tree2(worldbundle *wb, real t,
				bool vel)		// default = false
{
    static worldbundle *wb_last = NULL;		// last worldbundle handled
    static real t_int = -VERY_LARGE_NUMBER;	// last interpolation time

    static pdyn *root = NULL;			// root node of the
    						// interpolated tree

    if (!wb_last) set_kepler_fast_flag();	// use the "fast" kepler solver

    // Try to take care of rounding error in t.

    if (!trim(wb, t)) return NULL;

    if (wb != wb_last) {
	if (wb->get_pdyn_root()) {
	    root = wb->get_pdyn_root();		// restore an existing tree
	    t_int = wb->get_t_int();
	} else {
	    root = NULL;			// build a new tree
	    t_int = -VERY_LARGE_NUMBER;
	}
    }

    if (t == t_int && root) return root;

    // We no longer create the tree from scratch each time.  Rather,
    // we modify all or part of an existing tree (starting at root)
    // if possible.

    // All leaves on the bundle list are still current, by construction.
    // Find the base segment corresponding to each and use it to update
    // the pdyn tree interpolated to the current time.

    bool debug = false;

    // Build the new tree.

    interpolate_tree(wb, t, t_int, root, vel, debug);

    wb->set_t_int(t);
    if (wb != wb_last) wb->set_pdyn_root(root);

    wb_last = wb;
    t_int = t;

    return root;
}

void preload_pdyn(worldbundleptr wh[], int nh,
		  bool verbose)			// default = false
{
    for (int i = 0; i < nh; i++) {
	create_interpolated_tree2(wh[i], wh[i]->get_t_min());
	if (verbose)
	    cerr << "allocated memory for worldbundle " << i
		 << ",  t_min = " << wh[i]->get_t_min() << endl;
    }
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

    pdyn *root = create_interpolated_tree2(wb, t, true);
    put_node(cout, *root, false);

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
	    print_event(w, s->get_first_event(), t);
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
	pdyn *root = create_interpolated_tree2(wb, t);
	put_node(cout, *root, false);
	t += 0.01;
    }
#endif

}

#endif
