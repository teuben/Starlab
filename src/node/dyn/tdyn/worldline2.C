
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// worldline2:  Functions to manage and manipulate 4-D trees
////             Defines worldbundle, worldline, and segment
////             classes and member functions.
////
////             First argument is input data file.
////             No options.
//
//----------------------------------------------------------------------
//
// Other externally visible functions:
//
//	pdyn *create_interpolated_tree2(worldbundle *wb, real t)
//
//----------------------------------------------------------------------

#include "worldline.h"

#ifndef TOOLBOX

local void print_events(segment *s, real t)
{
    PRI(4); cerr << "events for s = " << s << ":" << endl;

    tdyn *bn = s->get_first_event();
    PRI(4); cerr << "base node "; PRC(bn); PRC(bn->get_time());
    PRL(bn->format_label());

    tdyn *b = bn;
    while (b->get_time() < t) {
	PRI(8); PRC(b); PRL(b->get_time());
	if (b->get_next()) b = b->get_next();
    }
    PRI(8); PRC(b); PRC(b->get_time());
    if (b->get_next()) {
	PRL(b->get_next()->get_time());
    } else {
	cerr << "next = NULL" << endl << endl;
    }
}

local void print_details(worldbundle *wb, tdyn *p, real t)
{
    // Re-locate node p at time t and print out relevent information
    // on the local worldline/segment/event structure.

    cerr << "details..." << endl;

    worldline *w = wb->find_worldline(p);
    segment *s = w->get_first_segment();
    segment *sprev = NULL;

    PRI(4); PRL(p->format_label());
    PRI(4); PRC(w); PRL(s);

    while (s->get_t_end() < t) {
	sprev = s;
	s = s->get_next();
    }

    PRI(4); PRL(s);

    PRI(4); PRC(t); PRC(s->get_t_start()); PRL(s->get_t_end());
    if (sprev) {
	PRI(4); PRC(sprev->get_t_start()); PRL(sprev->get_t_end());
	PRI(4); PRL(sprev->get_first_event());
    }

    print_events(s, t);
}

//======================================================================

// Creation of an entire interpolated tree at some specific time.

local void print_binary_diagnostics(real t, tdyn *b, tdyn *bb, pdyn *curr)
{
    PRC(t); PRL(b->format_label());
    PRL(bb->get_time());
    PRL(b->get_time());

    tdyn *n = b->get_next();
    if (n) PRL(n->get_time());

    real fac = 1 + bb->get_mass()
		/ bb->get_binary_sister()->get_mass();
    real M = bb->get_parent()->get_mass();

    real rbb = abs(fac*bb->get_pos());
    real vbb = square(fac*bb->get_vel());
    real ebb = 0.5*vbb - M / rbb;

    real rb = abs(fac*b->get_pos());
    real vb = square(fac*b->get_vel());
    real eb = 0.5*vb - M / rb;

    real rc = abs(fac*curr->get_pos());
    real vc = square(fac*curr->get_vel());
    real ec = 0.5*vc - M / rc;

    PRI(4); PRC(rbb); PRC(vbb); PRL(ebb);
    PRI(4); PRC(rb); PRC(vbb); PRL(eb);
    PRI(4); PRC(rc); PRC(vc); PRL(ec);

    if (n) {
	real rn = abs(fac*n->get_pos());
	real vn = square(fac*n->get_vel());
	real en = 0.5*vn - M / rn;
	PRI(4); PRC(rn); PRC(vn); PRL(en);
    }

    PRI(4); PRL(-M/(2*ebb));

    PRL(interpolate_pos(b, b->get_time(), bb));
    PRL(interpolate_pos(b, t, bb));
    if (n)
	PRL(interpolate_pos(b, n->get_time(),
			    bb));

    real dt = t - (real)b->get_time();
    real dtn = 0;
    if (n) dtn = n->get_time() - b->get_time();
    cerr << "interpolation..." << endl;
    PRC(dt); PRL(dtn);
    PRL(b->get_pos());
    PRL(dt*b->get_vel());
    PRL(dt*dt*b->get_acc());
    PRL(dt*dt*dt*b->get_jerk());
}

local inline pdyn *attach_new_node(worldbundle *wb, worldline *ww,
				   pdyn *root, tdyn *top, tdyn *bb,
				   bool debug)
{
    // Attach a new node corresponding to worldline ww.  The root node
    // for the interpolated tree is root; bb is the base node for worldline
    // ww; top is its top-level base node (at the start of the segment).

    // The traversal of the tree is such that parents are always seen
    // before children and elder sisters are always seen before younger
    // sisters.

    if (debug)
	cerr << "base node " << bb << " "
	     << bb->get_time() << " "
	     << bb->format_label() << endl;

    pdyn *curr = new pdyn(NULL, NULL, false);
    cerr << "created new node for " << bb->format_label() << endl;

    // Don't use add_node, as it will create a tree with nodes in the
    // reverse order!

    if (bb == top) {

	// Attach the top-level node.
	// Add curr to the end of the top-level list.

	pdyn *n = root->get_oldest_daughter();

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

	// Attach a low-level node.

	// Parent and elder sister pointers are derived
	// (awkwardly!) from the existing tree structures.

	pdyn *par = wb->find_worldline(bb->get_parent())
	    	      ->get_tree_node();
	curr->set_parent(par);

	tdyn * bb_sis = bb->get_elder_sister();

	if (!bb_sis)
	    par->set_oldest_daughter(curr);
	else {
	    pdyn *sis = wb->find_worldline(bb_sis)
			  ->get_tree_node();
	    curr->set_elder_sister(sis);
	    sis->set_younger_sister(curr);
	}
    }

    ww->set_tree_node(curr);
    return curr;
}

local inline void update_node(worldbundle *wb,
			      worldline *ww, real t,
			      tdyn *bb, tdyn *top, bool vel,
			      bool debug)
{
    // Copy or interpolate all relevant quantities from the base node
    // bb to the node curr in the interpolated tree.

    tdyn *b = find_event(bb, t);
    pdyn *curr = ww->get_tree_node();

    if (debug)
	cerr << "current node " << b << " "
	     << b->get_time() << " "
	     << b->format_label() << endl;

    // Very helpful to attach the worldline ID as an index to the pdyn,
    // but then we lose the connection between the index seen by the
    // display program and the index in the original data.  Save both!

    // Index is the original index; worldline ID is worldline_index.

    if (b->get_name()) {
	curr->set_name(b->get_name());
	curr->set_index(atoi(b->get_name()));
    }
    if (b->get_index() >= 0)
	curr->set_index(b->get_index());

    // Cleanest way to get the worldline index:

    curr->set_worldline_index(wb->find_index(b));

    curr->set_mass(b->get_mass());

    // Note that index and mass (and some others) do not need to
    // be updated -- should be set when the node is created.

    if (!b->get_kepler())
	curr->set_pos(interpolate_pos(b, t, bb));
    else
	curr->set_pos(b->get_pos());

    if (bb != top || vel) {

	// Need velocity information on low-level nodes.
	// Very inefficient!

	if (!b->get_kepler())
	    curr->set_vel(interpolate_vel(b, t, bb));
	else
	    curr->set_vel(b->get_vel());

#if 0
	// Diagnostic output, invoked for the elder sister of a
	// perturbed binary.

	if (!b->get_elder_sister()
	    && !b->get_kepler())
	    print_binary_diagnostics(t, b, bb, curr);
#endif

    }

    if (b->get_kepler())
	curr->set_kepler((kepler*)1);

    // **** NEW for pdyn data (Steve, 5/01). ****

    curr->set_stellar_type(b->get_stellar_type());
    curr->set_temperature(b->get_temperature());
    curr->set_luminosity(b->get_luminosity());

    ww->set_t_curr(t);

    // cerr << "updated " << curr->format_label()
    //      << endl;
}

local inline void update_interpolated_tree(worldbundle *wb,
					   worldline *w, segment *s,
					   pdyn *root, real t, real t_int,
					   bool vel, bool debug)
{
    // Current worldline is w, segment is s, and w has not yet been
    // updated, according to t_curr.

    // *** By construction, on entry, w and s represent a leaf. ***

    bool rebuild = false;

    if (t_int < s->get_t_start() || t_int > s->get_t_end()) {

	// Times t and t_int lie in different segments, so we must
	// rebuild the portion of the tree containing this particle.

	// Delete the old subtree.  The tree_node() pointer still
	// points into the old tree, so use this to locate and delete
	// the out-of-date structure.

	pdyn *old = w->get_tree_node();

	if (old) {
	    old = old->get_top_level_node();	    // should be OK, as
	    					    // pdyn root is defined

	    // Zero out all tree_node() pointers, since it is likely
	    // that not all will be reset during this function call.

	    int ndel = 0;

	    for_all_nodes(pdyn, old, o) {

		cerr << "deleting " << o->format_label() << endl;
		ndel++;

		int wi = o->get_worldline_index();
		PRL(wi);
		if (wi <= 0)
		    cerr << "update_interpolated_tree: error: "
			 << "deleting node with no worldline index." << endl;
		else
		    wb->get_bundle()[wi]->clear_tree_node();
	    }		
	    
	    detach_node_from_general_tree(*old);
	    rmtree(old);

	    cerr << "deleted " << ndel << " nodes" << endl;
	}

	rebuild = true;
    }

    // Rebuild/update the current subtree.  Start by finding the base node
    // (containing all relevant tree information) for this particle, and
    // the current event.

    tdyn *bn = s->get_first_event();

    if (debug) {
	cerr << "s = " << s << endl
	     << "first event = " << bn << " "
	     << bn->format_label() << " "
	     << bn->get_time() << endl;

	print_details(wb, bn, t);
    }

    // Create the portion of the tree corresponding to the top-level
    // node of bn.  Note that top will probably be bn most of the time.

    // Need to be careful with top-level tdyn nodes, as they
    // generally won't have parent nodes (no root), so standard
    // functions like is_top_level_node() and get_top_level_node()
    // will fail.

    tdyn *top = bn;
    while (top->get_parent()
	   && unique_id(top->get_parent()) > 0)
	top = top->get_parent();

    if (debug) {
	PRC(top); PRL(top->format_label());
	// put_node(cerr, *top, false);
    }

    // Copy the entire tree below top into the new tree.

    for_all_nodes(tdyn, top, bb) {

	// Find the worldline of bb.

	worldline *ww;

	if (bb == bn)
	    ww = w;
	else {
	    ww = wb->find_worldline(bb);
	    if (!ww)
		cerr << "create_interpolated_tree2: error:"
		     << "can't find worldline of subtree member."
		     << endl;
	}

	if (ww) {

	    // Create and attach a new node, if necessary.

	    if (rebuild) attach_new_node(wb, ww, root, top, bb, debug);

	    // Update the tree entry: copy or interpolate all relevant
	    // quantities from bb to curr.

	    if (ww->get_tree_node())
		update_node(wb, ww, t, bb, top, vel, debug);
	    else
		cerr << "update_interpolated_tree: error: no tree_node"
		     << " for " << bb->format_label()
		     << endl;
	}
    }
}

#define EPS 1.e-12

pdyn *create_interpolated_tree2(worldbundle *wb, real t,
				bool vel)		// default = false
{
    static pdyn *root = NULL;			// root node of the
    						// interpolated tree

    static real t_int = -VERY_LARGE_NUMBER;	// time last the interpolated
						// tree was updated

    if (t == t_int && root) return root;

    // We no longer create the tree from scratch.  Rather, we modify
    // all or part of an existing tree (starting at root) if possible.

    // All leaves on the bundle list are still current, by construction.
    // Find the base segment corresponding to each and use it to update
    // the pdyn tree interpolated to the current time.

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

    // Create a root node.

    if (!root) root = new pdyn(NULL, NULL, false);
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

	if (w->get_t_curr() != t) {

	    // Worldline needs to be updated.

	    real id = w->get_id();

	    if (id >= 1 && id < 2) {

		// Worldline w represents a leaf.  Locate its current segment.

		if (debug) {
		    cerr.precision(16);
		    cerr << endl
			 << "updating worldline " << i << ", id = "
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
		// should be continuous -- check in debugging mode...

		if (s)
		    update_interpolated_tree(wb, w, s,
					     root, t, t_int,
					     vel, debug);
		else if (debug)
		    cerr << "NULL segment for leaf..." << endl;
	    }
	}
    }

    t_int = t;
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

    pdyn *root = create_interpolated_tree2(wb, t, true);
    put_node(cout, *root, false);
    //rmtree(root);

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
	pdyn *root = create_interpolated_tree2(wb, t);
	put_node(cout, *root, false);
	//rmtree(root);
	t += 0.01;
    }
#endif

}

#endif
