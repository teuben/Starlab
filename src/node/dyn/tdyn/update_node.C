local inline void update_node(worldbundle *wb,
			      worldline *ww, real t,
			      tdyn *bb, tdyn *top, bool vel,
			      bool debug)
{
    // Copy or interpolate all relevant quantities from the base node
    // bb to the node curr in the interpolated tree.  For convenience,
    // top is the top-level node of bb.

    tdyn *b = find_event(bb, t);
    pdyn *curr = ww->get_tree_node();

    if (debug)
	cerr << "current node " << b << " "
	     << b->get_time() << " "
	     << b->format_label() << endl;


    if (NEW == 0) {

	//-----------------------------------------------------------------
	// The following "static" quantities should be constant within a
	// segment, and thus need only be set when node curr is created.
	// Should not be necessary to copy these data at each update.
	//
	//		name/index
	//		worldline index
	//-----------------------------------------------------------------

	// Moved to attach_new_node.C for NEW = 1 (Steve, 5/30/01):

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
    }

    // Allow possibility that mass may change along the worldline:

    curr->set_mass(b->get_mass());

    // See if we need to deal with unperturbed motion.

    // The kepler flag (1) is attached once the motion is unperturbed,
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

    if (b->get_kepler() == (kepler*)1) {	// 2 now reserved for lightly
						// perturbed binary motion
						// -- not yet implemented
	if (!curr->get_kepler()) {

	    // Default placeholder:

	    curr->set_kepler((kepler*)1);

	    if (NEW == 1) {

		// Locate b's sister node.  Tree info is attached to the
		// base nodes.  Compute a kepler only for the elder sister
		// of a pair of nodes.  Younger sister just gets a "1".

		if (!bb->get_elder_sister()) {

		    tdyn *bbsis = bb->get_younger_sister();
		    tdyn *sis = find_event(bbsis, t);

		    // Check that sis has the same time and also has kep set.

		    if (sis->get_time() != b->get_time())
			cerr << "warning: unsynchronized binary nodes" << endl;

		    if (!sis->get_kepler())
			cerr << "warning: binary sister has kep = NULL" << endl;

		    // Create a real kepler structure.

//		    cerr << "creating new kepler for " << bb->format_label()
//			 << " at time " << b->get_time() << endl;

		    kepler * k = new kepler;

		    k->set_circular_binary_limit(1.e-6);	// ~arbitrary
		    set_kepler_tolerance(2);			// repetitious

		    k->set_time(b->get_time());
		    k->set_total_mass(b->get_mass() + sis->get_mass());
		    k->set_rel_pos(sis->get_pos() - b->get_pos());
		    k->set_rel_vel(sis->get_vel() - b->get_vel());

		    k->initialize_from_pos_and_vel();

		    curr->set_kepler(k);

		    if (debug)
			cerr << "created kepler for " << curr->format_label()
			     << " at time " << b->get_time() << endl;
		}
	    }
	}

    } else {

	// Delete any old kepler (if real).

	if (curr->get_kepler() && curr->get_kepler() != (kepler*)1) {
	    delete curr->get_kepler();

	    if (debug && 
		!bb->get_elder_sister())
		cerr << "deleted kepler for " << curr->format_label()
		     << " at time " << t << " (" << b->get_time() << ")"
		     << endl;
	}

	curr->set_kepler(NULL);
    }

    // Now actually do the interpolation.  Take advantage of the binary
    // tree structure, and do nothing for the younger sister of a pair.

    kepler *k = curr->get_kepler();

    if (!k) {

	// This is a top-level node or a perturbed binary component.
	// New code: Take no action if this is the younger component
	// of a binary.

	if (NEW == 0 || bb == top || !bb->get_elder_sister())
	    curr->set_pos(interpolate_pos(b, t, bb));

	// Need velocity information on low-level nodes, or if forced:

	if (vel || (bb != top && (NEW == 0 || !bb->get_elder_sister())))
	    curr->set_vel(interpolate_vel(b, t, bb));

#if 0
	if (!bb->get_elder_sister())
	    print_binary_diagnostics(t, b, bb, curr);
#endif

    } else if (k != (kepler*)1) {

	// This is the elder component of an unperturbed binary, and
	// has a real kepler attached (should happen only for NEW = 1)

	k->transform_to_time(t);

	// By construction, rel_pos runs from curr to the
	// other (younger) binary component.

	real mass_fac = 1 - curr->get_mass()/k->get_total_mass();
	curr->set_pos(-mass_fac * k->get_rel_pos());
	curr->set_vel(-mass_fac * k->get_rel_vel());

    } else if (NEW == 0) {

	// Default case (old code only):

	curr->set_pos(b->get_pos());
	curr->set_vel(b->get_vel());
    }

    if (NEW == 1) {

	// Update the younger component of a binary (perturbed or not).

	if (bb != top && !bb->get_elder_sister()) {

	    // Find the sister node in the interpolated tree.
	    // Check each pointer; none should be NULL.

	    tdyn *bbsis = bb->get_younger_sister();	     // sister base node
	    if (bbsis) {
		worldline *wsis = wb->find_worldline(bbsis); // sister worldline
		if (wsis) {
		    pdyn *sis = wsis->get_tree_node();
		    if (sis) {
			if (sis->get_mass() > 0) {
			    real mass_fac = -curr->get_mass()/sis->get_mass();
			    sis->set_pos(mass_fac*curr->get_pos());
			    sis->set_vel(mass_fac*curr->get_vel());
			} else
			    cerr << "update_node: "
				 << "error: sister mass <= 0" << endl;
		    } else
			cerr << "update_node: "
			     << "error: sister tree node NULL for "
			     << bbsis->format_label() << " at " << t << endl;
		} else
		    cerr << "update_node: "
			 << "error: sister worldline NULL" << endl;
	    } else
		cerr << "update_node: "
		     << "error: base node sister NULL" << endl;
	}
    }

    curr->set_stellar_type(b->get_stellar_type());
    curr->set_temperature(b->get_temperature());
    curr->set_luminosity(b->get_luminosity());

    ww->set_t_curr(t);

    if (debug)
	cerr << "updated " << curr->format_label() << endl;
}
