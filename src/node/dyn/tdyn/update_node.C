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

    if (debug && curr->get_index() == -42) {	    // from attach_new_node
	cerr << "added "; PRL(b->format_label());
    }

    //-----------------------------------------------------------------
    // The following "static" quantities should be constant within a
    // segment, and thus need only be set when node curr is created.
    // Should not be necessary to copy these data at each update.
    //
    //		name/index
    //		mass

    // Very helpful to attach the worldline ID as an index to the pdyn,
    // but then we lose the connection between the index seen by the
    // display program and the index in the original data.  Save both!

    // Index is the original index; worldline ID is worldline_index.

#if 0

    // Moved to attach_new_node.C (Steve, 5/30/01):

    if (b->get_name()) {
	curr->set_name(b->get_name());
	curr->set_index(atoi(b->get_name()));
    }
    if (b->get_index() >= 0)
	curr->set_index(b->get_index());

    // Cleanest way to get the worldline index:

    curr->set_worldline_index(wb->find_index(b));

    curr->set_mass(b->get_mass());

#endif

    // The kepler flag is attached once the motion is unperturbed,
    // but unperturbed motion doesn't start a new segment.

    // Assume that the unperturbed approximation is OK even if only
    // the start of the current section of the worldline (node b) is
    // flagged as unperturbed:
    //
    //		    b		   b->next
    //
    //		(unpert) -----x--- (unpert)	OK -- normal case
    //
    //		(unpert) -----x--- (pert)	OK -- motion became perturbed
    //						at the end of the step
    //
    //		  (pert) -----x--- (unpert)	not OK -- motion perturbed
    //						throughout the step

    if (b->get_kepler()) {

	if (!curr->get_kepler()) {

	    // Default placeholder:

	    curr->set_kepler((kepler*)1);

#if 1
	    // Locate b's sister node.  Tree info is attached to the
	    // base nodes.

	    tdyn *bbsis = bb->get_binary_sister();
	    tdyn *sis = find_event(bbsis, t);

	    // Check that sis has the same time and also has kep set.

	    if (sis->get_time() != b->get_time())
		cerr << "warning: unsynchronized binary nodes" << endl;

	    if (!sis->get_kepler())
		cerr << "warning: binary sister has kep = NULL" << endl;

	    // Create a real kepler structure.  For now, maintain
	    // separate (redundant) keplers for each component.

	    kepler * k = new kepler;

	    k->set_time(b->get_time());
	    k->set_total_mass(b->get_mass() + sis->get_mass());
	    k->set_rel_pos(sis->get_pos() - b->get_pos());
	    k->set_rel_vel(sis->get_vel() - b->get_vel());

	    k->initialize_from_pos_and_vel();

	    curr->set_kepler(k);

	    if (debug)
		cerr << "created kepler for " << curr->format_label()
		     << " at time " << b->get_time() << endl;
#endif
	}

    } else {

	// Delete any old kepler (if real).

	if (curr->get_kepler() && curr->get_kepler() != (kepler*)1) {
	    delete curr->get_kepler();

	    if (debug)
		cerr << "deleted kepler for " << curr->format_label()
		     << " at time " << t << " (" << b->get_time() << ")"
		     << endl;
	}

	curr->set_kepler(NULL);
    }

    //-----------------------------------------------------------------

    // Note that we currently treat binary coponents independently,
    // doubling the work done at each time step...

    kepler *k = curr->get_kepler();
    real mass_fac;

    if (!k)

	curr->set_pos(interpolate_pos(b, t, bb));

    else {
	if (k == (kepler*)1)

	    curr->set_pos(b->get_pos());

	else {

	    // Unperturbed motion:

	    k->transform_to_time(t);

	    // By construction (for now), rel_pos runs from curr to
	    // the other component.

	    mass_fac = 1 - curr->get_mass()/k->get_total_mass();

	    curr->set_pos(-mass_fac * k->get_rel_pos());
	}
    }

    if (bb != top || vel) {

	// Need velocity information on low-level nodes.
	// Very inefficient!

	if (!k)

	    curr->set_vel(interpolate_vel(b, t, bb));

	else {

	    if (k == (kepler*)1)
		curr->set_vel(b->get_vel());
	    else
		curr->set_vel(-mass_fac * k->get_rel_vel());
	}

#if 0
	// Diagnostic output, invoked for the elder sister of a
	// perturbed binary.

	if (!b->get_elder_sister()
	    && !b->get_kepler())
	    print_binary_diagnostics(t, b, bb, curr);
#endif

    }

    // **** NEW for pdyn data (Steve, 5/01). ****

    curr->set_stellar_type(b->get_stellar_type());
    curr->set_temperature(b->get_temperature());
    curr->set_luminosity(b->get_luminosity());

    ww->set_t_curr(t);

    // cerr << "updated " << curr->format_label()
    //      << endl;
}
