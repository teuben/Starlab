local inline void attach_new_node(worldbundle *wb, worldline *ww,
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
}
