local INLINE void set_kepler_from_tdyn_pair(kepler *k, tdyn *b, tdyn *sis)
{
    k->set_time(b->get_time());
    k->set_total_mass(b->get_mass() + sis->get_mass());
    k->set_rel_pos(sis->get_pos() - b->get_pos());
    k->set_rel_vel(sis->get_vel() - b->get_vel());
    k->initialize_from_pos_and_vel();
}

local INLINE void update_node(worldbundle *wb,
			      worldline *ww, real t,
			      tdyn *bb, tdyn *top, bool vel,
			      bool debug)
{
    // Copy or interpolate all relevant quantities from the base node
    // bb on worldline ww to the node curr in the interpolated tree.
    // For convenience, top is the top-level node of bb.

    tdyn *b = find_event(ww, bb, t);			// current event

    pdyn *curr = ww->get_tree_node();			// current node in the
							// interpolated tree

    if (debug) {
	cerr << "current node " << b << " "
	     << b->get_time() << " "
	     << b->format_label() << endl;
	cerr << "base node " << bb << " "
	     << bb->get_time() << " "
	     << bb->format_label() << endl;
    }

    //======================================================================

    // Preliminaries:

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

    //======================================================================

    // See if we need to deal with unperturbed  or lightly perturbed
    // motion (most of the work of this function).

    // Variables which may be set and used later...

    tdyn *bbsis = NULL;
    worldline *wwsis = NULL;
    tdyn *sis = NULL;

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
    //		  (pert) -----x--- (unpert)	not OK -- perturbed motion
    //						throughout the step
    //		  (pert) -----x--- (pert)	not OK -- perturbed motion

    if (b->get_kepler() == (kepler*)1) {	// 2 is reserved for lightly
						// perturbed binary motion

	// If the kepler flag for b is set to 1, we have unperturbed
	// motion during this interval.  If curr has no kepler, create
	// one.  If a kepler does exist, there is still the possibility
	// that it refers to a previous kepler segment and that we have
	// skipped a perturbed interval.  Must check this and update the
	// kepler data if necessary.			(Steve, 9/01)

	kepler *k = curr->get_kepler();

	if (!k) {

	    // Default placeholder:

	    curr->set_kepler((kepler*)1);

	    if (NEW == 1) {

		// Locate b's sister node.  Tree info is attached to the
		// base nodes.  Compute a kepler only for the elder sister
		// of a pair of nodes.  Younger sister just gets a "1".

		if (!bb->get_elder_sister()) {

		    bbsis = bb->get_younger_sister();
		    wwsis = wb->find_worldline(bbsis);
		    sis = find_event(wwsis, bbsis, t);

		    // Check that sis has the same time and also has kep set.

		    if (sis->get_time() != b->get_time())
			cerr << "warning: unsynchronized binary nodes" << endl;

		    if (!sis->get_kepler())
			cerr << "warning: binary sister has kep = NULL" << endl;

		    // Create a real kepler structure.

		    // cerr << "creating new kepler for " << bb->format_label()
		    //      << " at time " << b->get_time() << endl;

		    k = new kepler;

		    k->set_circular_binary_limit(1.e-6);	// ~arbitrary
		    set_kepler_tolerance(2);			// repetitious

		    set_kepler_from_tdyn_pair(k, b, sis);
		    curr->set_kepler(k);

		    if (debug)
			cerr << "created kepler for " << curr->format_label()
			     << " at time " << b->get_time() << endl;
		}
	    }

	} else {

	    if (NEW == 1) {

		// A kepler structure already exists.  Check to see if it
		// is consistent with the current pos and vel.  May be
		// expensive to do this every step -- need a better way of
		// verifying the consistency of the data.

		if (!bb->get_elder_sister()) {

		    // Locate b's sister node.  Tree info is attached to the
		    // base nodes.  Code follows that above (kepler == NULL).
		    // Note that this isn't so expensive, as all these data
		    // should be used again below.

		    bbsis = bb->get_younger_sister();
		    wwsis = wb->find_worldline(bbsis);
		    sis = find_event(wwsis, bbsis, t);

		    // Check that sis has the same time and also has kep set.

		    if (sis->get_time() != b->get_time())
			cerr << "warning: unsynchronized binary nodes" << endl;

		    if (!sis->get_kepler())
			cerr << "warning: binary sister has kep = NULL" << endl;

		    // For now, use angular momentum as a check -- better than
		    // energy, which is likely to be nearly conserved during an
		    // intervening period of perturbed motion.

		    // m1r1 + m2r2 = 0
		    // r2 = -m1r1/m2
		    // r1 - r2 = r1 (1 + m1/m2)
		    // v1 - v2 = v1 (1 + m1/m2)
		    // (r1-r2)^(v1-v2) = (1+m1/m2)^2 r1^v1

		    // Hard to see how to avoid this cost...

		    vec ang_mom = square(1+b->get_mass()/sis->get_mass())
					* (b->get_pos() ^ b->get_vel());

#define TWFAC 1.e-6
		    if (!twiddles(square(ang_mom),
				  square(curr->get_kepler()
					  ->get_angular_momentum()), TWFAC)) {

			// Old and new orbits are inconsistent.  Redefine
			// the orbital parameters.  Factor TWFAC comes from
			// trial and error, but seems to work.

			set_kepler_from_tdyn_pair(k, b, sis);

			if (debug)
			    cerr << "recreated kepler for "
				 << curr->format_label()
				 << " at time " << b->get_time() << endl;
		    }
		}
	    }
	}

#if 0		// turn this on to enable new kep2 stuff, but the new code
		// is incompatible with the old buggy kira output format...

    // (thinking aloud...)
    //
    // What about lightly perturbed (lpert) binaries (kep = 2)?
    // Outcome is pretty simple:
    //
    //		    b		   b->next
    //
    //		(unpert) -----x--- (unpert)	unperturbed	    (1 kepler)
    //		(unpert) -----x--- (lpert)	unperturbed
    //		(unpert) -----x--- (pert)	unperturbed
    //
    //		 (lpert) -----x--- (unpert)	lightly perturbed   (2 keplers)
    //		 (lpert) -----x--- (pert)	lightly perturbed
    //		 (lpert) -----x--- (lpert)	lightly perturbed
    //
    //		  (pert) -----x--- (unpert)	perturbed	    (0 keplers)
    //		  (pert) -----x--- (lpert)	perturbed
    //		  (pert) -----x--- (pert)	perturbed
    //
    // i.e. only b matters.
    //
    // Now curr must keep track of 2 kepler pointers (b and next) in
    // some coherent and efficient way.  Interpolation will entail
    // interpolation between keplers at two different times.
    //
    //		kepler #1 is b,  kepler #2 is b->next
    //
    // As with unperturbed case, keep track of existing keplers, for
    // efficiency.  Need new data structure to hold 2 kepler pointers.
    // As of 9/01, the pdyn class contains two kepler pointers: kep and
    // kep2.						(Steve, 9/01)
    //
    // Curr has 0 keplers:	create 2
    //		1 kepler:	was unperturbed: check 1, create 1
    //		2 keplers:	already lightly perturbed; check 2
    //
    // We expect that both keplers may change from one step to the next,
    // and the "lightly perturbed" period will cover several steps, so
    // for now (maybe forever) just do the following:
    //
    //		- if 0 or 1 kepler exists, create two new ones
    //		  (don't even bother checking to see if the unpert
    //		   kepler is usable?)			kep2 == NULL
    //
    //		- if 2 keplers exist, see if any can be used (same
    //		  events, indicated by saving time or event address),
    //		  and create new ones as necessary:	kep2 != NULL

    } else if (b->get_kepler() == (kepler*)2) {

	// Lightly perturbed motion.

	// Take action only if this is the elder sister of a binary.
	// Younger sister pointers are also set in the process.

	if (!bb->get_elder_sister()) {

	    // Default action is to create two kepler structures.

	    bool createkep = true, createkep2 = true;
	    tdyn *next = b->get_next();

	    if (curr->get_kepler2()) {

		// Two keplers already exist (assume that kepler2 implies the
		// existence of kepler).  See if either is still valid.

		tdyn *event = curr->get_kepevent();
		tdyn *event2 = curr->get_kepevent2();

		if (event == b)

		    createkep = false;

		else if (event == next) {

		    // Old event will become the new event2.
		    // Remove the old kep2 (event2 shouldn't be b!).

		    delete curr->get_kepler2();

		    // New kep2 is old kep.

		    curr->set_kepler2(curr->get_kepler());
		    pdyn *sister = curr->get_younger_sister();
		    sister->set_kepler2(curr->get_kepler());

		    curr->set_kepler(NULL);
		    sister->set_kepler(NULL);

		    createkep2 = false;
		}

		if (event2 == next)

		    createkep2 = false;

		else if (event2 == b) {

		    // Old event2 will become the new event.
		    // Remove the old kep (event shouldn't be next!).

		    delete curr->get_kepler();

		    // New kep is old kep2.

		    curr->set_kepler(curr->get_kepler2());
		    pdyn *sister = curr->get_binary_sister();
		    sister->set_kepler(curr->get_kepler2());

		    curr->set_kepler2(NULL);
		    sister->set_kepler2(NULL);

		    createkep = false;
		}
	    }

	    if (createkep || createkep2) {

		// Must construct at least one kepler structure.

		bbsis = bb->get_younger_sister();
		wwsis = wb->find_worldline(bbsis);
		sis = find_event(wwsis, bbsis, t);

		// Check that sis has the same time and also has kep set.
		// (Same check and output as above.)

		if (sis->get_time() != b->get_time())
		    cerr << "warning: unsynchronized binary nodes" << endl;

		if (!sis->get_kepler())
		    cerr << "warning: binary sister has kep = NULL" << endl;

		if (createkep) {
		    curr->rmkepler();
		    kepler *k = new kepler;
		    set_kepler_from_tdyn_pair(k, b, sis);
		    curr->set_kepler(k);
		}
			
		if (createkep2) {
		    curr->rmkepler2();
		    kepler *k = new kepler;
		    set_kepler_from_tdyn_pair(k, next, sis->get_next());

		    // sis->next should be OK...

		    curr->set_kepler2(k);
		}
	    }

	    curr->set_kepevent(b);
	    curr->set_kepevent2(next);	// don't bother setting sister pointers
	}

#endif

    } else {						// unperturbed motion

	// Delete any old keplers (if real).

	if (curr->get_kepler() == (kepler*)2)		// shouldn't occur
	    PRL(curr->get_kepler());

	if (curr->get_kepler()
	     && curr->get_kepler() != (kepler*)1
	     && curr->get_kepler() != (kepler*)2) {	// "2" shouldn't be
							// needed, but...
	    curr->rmkepler();

	    if (debug && 
		!bb->get_elder_sister())
		cerr << "deleted kepler for " << curr->format_label()
		     << " at time " << t << " (" << b->get_time() << ")"
		     << endl;
	}

	if (curr->get_kepler2()) {

	    curr->rmkepler2();

	    if (debug && 
		!bb->get_elder_sister())
		cerr << "deleted kepler2 for " << curr->format_label()
		     << " at time " << t << " (" << b->get_time() << ")"
		     << endl;
	}
    }

    //======================================================================

    // Now actually do the interpolation.  Take advantage of the binary
    // tree structure, and do nothing for the younger sister of a pair.

    kepler *k = curr->get_kepler();

    if (!k || k == (kepler*)2) {		// *** still to be completed ***

	// This is a top-level node or a perturbed binary component.
	// New code: Take no action if this is the younger component
	// of a binary.

	if (NEW == 0 || bb == top || !bb->get_elder_sister())
#ifndef NEW_INTERP
	    curr->set_pos(interpolate_pos(b, t, bb));
#else
	    set_interpolated_pos(b, t, curr, bb);
#endif

	// Need velocity information on low-level nodes, or if forced:

	if (vel || (bb != top && (NEW == 0 || !bb->get_elder_sister())))
	    curr->set_vel(interpolate_vel(b, t, bb));

#if 0
	if (!bb->get_elder_sister())
	    print_binary_diagnostics(t, b, bb, curr);
#endif

    } else if (k != (kepler*)1 && !bb->get_elder_sister()) {

	// This is the elder component of an unperturbed binary, and
	// has a real kepler attached (should happen only for NEW = 1)

	k->transform_to_time(t);

	// By construction, rel_pos runs from curr to the
	// other (younger) binary component.

	real mass_fac = 1 - curr->get_mass()/k->get_total_mass();
	curr->set_pos(-mass_fac * k->get_rel_pos());
	curr->set_vel(-mass_fac * k->get_rel_vel());

    } else

	if (NEW == 0) {

	    // Default case (old code only):

	    curr->set_pos(b->get_pos());
	    curr->set_vel(b->get_vel());
	}

    //======================================================================

    pdyn *csis = NULL;

    if (NEW == 1) {

	// Update the younger component of a binary (perturbed or not).

	if (bb != top && !bb->get_elder_sister()) {

	    // Find the sister node in the interpolated tree.
	    // Check each pointer; none should be NULL.

	    if (!bbsis) {
		bbsis = bb->get_younger_sister();	     // sister base node
		wwsis = wb->find_worldline(bbsis);	     // sister worldline
	    }

	    if (bbsis && wwsis) {

		pdyn *csis = wwsis->get_tree_node();

		if (csis) {
		    if (csis->get_mass() > 0) {
			real mass_fac = -curr->get_mass()/csis->get_mass();
			csis->set_pos(mass_fac*curr->get_pos());
			csis->set_vel(mass_fac*curr->get_vel());
		    } else
			cerr << "update_node: "
			     << "error: sister mass <= 0" << endl;
		} else
		    cerr << "update_node: "
			 << "error: sister tree node NULL for "
			 << bbsis->format_label() << " at " << t << endl;
	    } else
		cerr << "update_node: "
		    << "error: sister base node sister NULL" << endl;
	}
    }

    //======================================================================

    // Clean up.

    curr->set_stellar_type(b->get_stellar_type());
    curr->set_temperature(b->get_temperature());
    curr->set_luminosity(b->get_luminosity());

    ww->set_t_curr(t);
    ww->set_current_event(b);				// for find_event()

    if (bbsis) {
	if (!sis) sis = find_event(wwsis, bbsis, t);
	wwsis->set_t_curr(t);
	wwsis->set_current_event(sis);			// for find_event()
	if (csis) {
	    csis->set_stellar_type(sis->get_stellar_type());
	    csis->set_temperature(sis->get_temperature());
	    csis->set_luminosity(sis->get_luminosity());
	}
    }

    if (debug)
	cerr << "updated " << curr->format_label() << endl;
}
