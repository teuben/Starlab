
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
//	real physical_mass_scale()
//
//----------------------------------------------------------------------

#include "worldline.h"

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

local inline void clean_up_subtree(worldbundle *wb, pdyn *old, bool debug)
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

	if (debug)
	    cerr << "deleted " << ndel << " nodes" << endl;
    }
}

local inline void update_interpolated_tree(worldbundle *wb,
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

    bool rebuild = (t_int <= s->get_t_start() || t_int >= s->get_t_end());

    // Rebuild/update the current subtree.  Start by finding the base
    // node (containing all relevant tree information) for this particle,
    // and the current event.

    // Need to be careful with top-level tdyn nodes, as they
    // generally won't have parent nodes (no root), so standard
    // functions like is_top_level_node() and get_top_level_node()
    // will fail.

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
}

#define EPS 1.e-12

pdyn *create_interpolated_tree2(worldbundle *wb, real t,
				bool vel)		// default = false
{
    static worldbundle *wb_last = NULL;		// last worldbundle to be
    						// handled

    static pdyn *root = NULL;			// root node of the
    						// interpolated tree

    static real t_int = -VERY_LARGE_NUMBER;	// time last the interpolated
						// tree was updated

    // Use the "fast" kepler solver...

    if (!root) set_kepler_fast_flag();

    if (wb != wb_last) {
	if (root) {
	    rmtree(root);
	    root = NULL;
	}
	for (int i = 0; i < wb->get_nw(); i++) {
	    wb->get_bundle()[i]->clear_tree_node();
	    wb->get_bundle()[i]->set_t_curr(-VERY_LARGE_NUMBER);
	}

	t_int = -VERY_LARGE_NUMBER;
    }

    if (t == t_int && root) return root;

    // We no longer create the tree from scratch each time.  Rather,
    // we modify all or part of an existing tree (starting at root)
    // if possible.

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
    bundle[0]->set_t_curr(t);

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

    wb_last = wb;
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
