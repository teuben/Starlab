local inline void attach_new_node(worldbundle *wb, worldline *ww,
				  pdyn *root, tdyn *top, tdyn *bb,
				  bool debug)
{
    // Create and attach a new node corresponding to worldline ww.
    // The root node for the interpolated tree is root; bb is the base
    // node for worldline ww; top is its top-level base node (at the
    // start of the segment).

    // The traversal of the tree is such that parents are always seen
    // before children and elder sisters are always seen before younger
    // sisters.

    if (debug)
	cerr << "base node " << bb << " "
	     << bb->get_time() << " "
	     << bb->format_label() << endl;

    pdyn *curr = new pdyn(NULL, NULL, false);

    // Attach curr to the tree.  Don't use add_node, as it will
    // create a tree with nodes in the reverse order!

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

    if (NEW == 1) {

	//-----------------------------------------------------------------
	// The following "static" quantities should be constant within a
	// segment, and thus need only be set when node curr is created.
	// Should not be necessary to copy these data at each update.
	//
	//		name/index
	//		mass
	//-----------------------------------------------------------------

	// Moved here from update_node.C (Steve, 5/30/01):

	// Very helpful to attach the worldline ID as an index to the pdyn,
	// but then we lose the connection between the index seen by the
	// display program and the index in the original data.  Save both!

	// Index is the original index; worldline ID is worldline_index.

	// Moved here from update_node.C (Steve, 5/30/01):

	if (bb->get_name()) {
	    curr->set_name(bb->get_name());
	    curr->set_index(atoi(bb->get_name()));
	}
	if (bb->get_index() >= 0)
	    curr->set_index(bb->get_index());

	// Cleanest way to get the worldline index:

	curr->set_worldline_index(wb->find_index(bb));

	curr->set_mass(bb->get_mass());
    }

    ww->set_tree_node(curr);

    // Set standard values for some quantities:

    ww->set_t_curr(-VERY_LARGE_NUMBER);
}
