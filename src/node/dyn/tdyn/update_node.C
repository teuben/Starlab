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
