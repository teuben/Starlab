
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Initiate, extend, or terminate slow binary motion.
//
// Externally visible functions:
//
//	void hdyn::startup_slow_motion
//	void hdyn::extend_or_end_slow_motion
//	bool has_binary_perturbers
//	void check_slow_consistency
//	void clear_perturbers_slow_perturbed
//
// Note from Steve, 8/99:
//
// Slow binary motion is handled by two data classes:
//
// 1. Slow motion itself is signified simply by attaching a slow_binary
//    pointer to the binary components (the data structure is shared
//    by the two components).  If the slow pointer is set, the binary
//    time scale is reduced by a slowdown factor kappa (i.e. dtau =
//    dt/kappa), and the perturbation is increased by the same factor.
//    The motion *must* continue for an integral number of orbits (see
//    Mikkola & Aarseth); it is reassessed for extension or termination
//    at the same phase (just after apocenter) of every orbit.  The
//    transition to unperturbed motion is handled by attaching a flag
//    indicating that the motion must terminate at the proper phase;
//    subsequently, unperturbed motion can begin normally.  This portion
//    of the motion is straightforward to program and maintain.
//
// 2. However, once a binary is undergoing slow motion, its center of
//    mass motion and its effect on its perturbers require special
//    attention.  The reason is that the internal motion causes a
//    "ripple" in the CM and external acc and jerk, and this must be
//    handled self-consistently.  When the internal motion is slowed, the
//    computed acc and jerk felt by the CM and its neighbors are not
//    consistent with the actual time development of the internal motion,
//    leading to large k (= j') and l (= j'') derivatives and excessively
//    short timesteps.  In fact, the opposite should be the case: the
//    slowed motion should *reduce* the ripple, allowing the CM and
//    neighbor steps to increase.  The corrector now treats interactions
//    involving slow binaries separately, reducing the relevant a, j, k,
//    and l by factors of kappa, kappa^2, etc., and does in fact allow
//    longer timesteps for the CM and perturber motion.
//
//    This portion of kira is significantly messier to code and maintain
//    than the slow motion itself.  The CM motion and perturber motion
//    are handled differently (and both are presently implemented only
//    for top-level nodes):
//
// 	* For the CM, instead of correcting the center of mass force
// 	  for the perturbation due to perturbers, we instead store the
// 	  perturbative part as acc_p and jerk_p in the slow_binary
// 	  structure attached to the *daughter* nodes.  Function
// 	  store_old_force also maintains old_acc_p and old_jerk_p in
// 	  this structure.  The corrector computes the higher
// 	  derivatives separately for the pertubative term, taking the
// 	  slowdown factors into account, and then combines it with the
// 	  usual correction terms.
//
// 	* For perturbers of a slow binary, things are more complex, as
// 	  the correction in principle is different for each slow
// 	  binary a particle happens to perturb.  Presently, I (Steve)
// 	  see no alternative to storing acc_p and jerk_p separately
// 	  (and maintaining the old_ versions) for every slow binary a
// 	  node perturbs.  Then, in the corrector phase, the high
// 	  derivatives are computed separately for each slow binary
// 	  interaction and combined into the total.  The acc_p and
// 	  jerk_p terms for each binary, together with other data
// 	  needed to manage them, are stored in the structure
// 	  slow_perturbed, arranged as a linked list attached to the
//	  *top-level* node of the perturber (since the correction is
//	  applied only to the motion of the top-level node, this seems
//	  the most natural place to store the information).  The
//	  structure also contains a pointer to the *center of mass*
//	  node of the slow binary in question.  Relevant list elements
// 	  are deleted each time the slow binary's perturber list is
//	  recomputed, and when slow motion ends.
//
// The procedures just described are currently implemented only for
// top-level nodes, principally because the perturbative component of
// the interparticle force is easy to determine in that case (since all
// forces on top-level nodes are initially computed in the center of
// mass approximation and the correction is easily absorbed into
// correct_acc_and_jerk).  It probably is possible to apply similar
// corrections for low-level binaries (e.g. in a triple), but these have
// not yet been implemented.  They probably should be implemented at
// least for binary sisters (to allow for, e.g. a triple containing a
// weakly perturbed binary).  For low-level slow binaries, it may be
// sufficient to use the kepler timestep criterion and avoid any further
// correction...  (Note that slow motion can always be applied; only the
// application of interparticle force corrections are affected.)
//
// Current conventions:
//
//	- slow motion may be applied to binaries at any level
//
//	- slow corrections are applied only to top-level nodes, so
//	  only top-level nodes should have slow_perturbed lists, and
//	  low-level binaries should not appear on those lists
//
//	- corrections to triple/multiple components are *not* yet
//	  implemented -- see hdyn_ev.C (correct_and_update)
//
//	  *** probably will use the slow binary structure attached ***
//	  *** to the (elder) binary component to implement this... ***
//
// The slow and slow_perturbed classes are defined in _dyn_.h and managed
// in _dyn_slow.C.  The only use is in kira, accessed via this file.

#include "hdyn.h"
#include "hdyn_inline.C"

bool has_binary_perturbers(hdyn* b)
{
    // Assume b is a low-level node (or a binary node) on entry.

    hdyn *pnode = b->find_perturber_node();

    if (pnode == NULL || !pnode->get_valid_perturbers()) {

	// Really should check that another binary exists somewhere in
	// the system, but this test is used to permit slow motion,
	// and pnode must exist in any case for that to happen.

	return true;
    }

    // See if any perturber is a binary CM or binary component.
    // Probably OK to allow unperturbed binaries to be perturbers, as
    // they will be handled in the CM approximation, but exclude them
    // for now.

    for (int k = 0; k < pnode->get_n_perturbers(); k++) {
	hdyn * p = pnode->get_perturber_list()[k];
	if (p->is_parent() || p->is_low_level_node()) return true;
    }

    return false;
}

// Store {2^k}^{1/3} in a fixed-length array (limits the slowdown factor
// to a maximum of 10^6!).  Array length is defined in kira_counters,
// but we insert an extra entry for kappa = 1 at the start.

static real k3[SLOW_BINS+1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

local void set_k3()
{
    int k2 = 1;
    for (int k = 0; k < SLOW_BINS; k++) {
	k3[k] = pow(k2, 1.0/3);
	k2 *= 2;
    }
}

local int get_slowdown_factor(hdyn* bi, real P = 0)
{
    // Initialize the local array.

    if (*k3 == 0) set_k3();

    // Require a valid list of perturbers.  Note that this will cause
    // existing slow motion to be terminated if no valid list is found.

    hdyn* pnode = bi->find_perturber_node();

    if (pnode == NULL || !pnode->get_valid_perturbers()) return -1;

    // if (has_binary_perturbers(bi)) return 0;

    // Estimate the number of orbital periods before the external
    // perturbation becomes unacceptably high, and set the binary
    // slowdown factor accordingly.

    // Don't use apert/jpert, since this is sensitive to the very
    // suborbital "wiggles" we are trying to average away!

    if (P <= 0)
	P = get_period(bi, bi->get_younger_sister());

    if (P <= 0 || P == VERY_LARGE_NUMBER) return -1;	// redundant...

    // Procedure: we can set the slowdown level to kappa if no member of
    // the perturber list will become a "kappa-perturber," with
    //
    //		kappa * perturbation  >  max_slow_perturbation
    //
    // during the next kappa orbits (i.e. in time kappa * P).  For each
    // perturber we compute the maximum permissible kappa.  The returned
    // value is the minimum of these kappa values.

    int kappa = bi->get_max_slow_factor();

    hdyn *par = bi->get_parent();
    real scale = binary_scale(par);	// "size" of the binary

    // Absolute pos, vel, and acc of the binary's top-level node:

    hdyn *top = bi->get_top_level_node();
    vector tpos = top->get_pred_pos();
    vector tvel = top->get_pred_vel();
    vector tacc = top->get_acc();

    // Loop over perturbers.

    for (int i = 0; i < pnode->get_n_perturbers(); i++) {

	hdyn *bj = pnode->get_perturber_list()[i];

	// Compute the critical separation corresponding to the maximum
	// allowed "slow" perturbation.  Note that, for some perturber
	// options, this quantity is independent of the perturber.
	// Live with that inefficiency (for now, at least).

	real rpert3 = crit_separation_cubed(par,
					    bj->get_mass(),
					    scale,
					    bi->get_max_slow_perturbation());

	real rpert = pow(rpert3, 1.0/3);	// live with this too...

	// (Note that rpert3 scales as 1 / gamma; hence the k3 = k^{1/3}
	// factor below.)

	// Compute relative pos, vel, and acc for use in the perturber
	// timescale estimates below.  Use top-level nodes if they are
	// different, in order to remove all internal motion from the
	// calculation.  If bi and bj are in the same clump, use par
	// and its binary sister in the calculation.

	vector dpos = 0, dvel = 0, dacc = 0;

	hdyn *jtop = bj->get_top_level_node();

	if (jtop != top) {
	    dpos = tpos - jtop->get_pred_pos();
	    dvel = tvel - jtop->get_pred_vel();
	    dacc = tacc - jtop->get_acc();
	} else {
	    hdyn *sis = par->get_binary_sister();
	    dpos = par->get_pos() - sis->get_pos();
	    dvel = par->get_vel() - sis->get_vel();
	    dacc = par->get_acc() - sis->get_acc();
	}

	real dr = abs(dpos);
	real vr = dvel * dpos / dr;
	real ar = dacc * dpos / dr;

	// Loop over k2 (= 2^k) and
	//
	//	1. determine the distance at which bj would be a
	//	   k2-perturber of pnode,
	//	2. estimate the time needed for bj to reach that
	//	   radius,
	//	3. stop when that time is less than k2*P.

	int k = 0, k2 = 1;
	while (k2 <= kappa) {

	    if (time_to_radius(dr - rpert * k3[k], vr, ar) < k2 * P)
		break;

	    k++;
	    k2 *= 2;
	}

	kappa = min(kappa, k2/2);
	if (kappa < 2) break;	
    }

    // *** Should we limit the slowdown increment to a factor of 2? ***

    return kappa;
}

// We could probably incorporate startup_slow_motion() into
// extend_or_end_slow_motion(), but leave the two functions
// separate for now...

void hdyn::startup_slow_motion()
{
    int k = get_slowdown_factor(this);

    if (k > 1) {

	create_slow(k);

	// Possibly should expand the perturber list to allow for
	// the slowdown factor, but that would increase the cost per
	// step and undo the benefit of the slowdown.  Assume for
	// now that we can ignore the extra weak perturbers so long
	// as we carry out this procedure one orbit at a time.

	if (diag->slow) {
	    cerr << endl
		 << "starting slow motion for "
		 << get_parent()->format_label()
		 << " at time " << get_time()
		 << "; kappa = " << get_kappa()
	         << endl;
	    cerr << "    perturbation = " << sqrt(perturbation_squared)
		 << endl
		 << "    timestep = " << timestep
		 << "  parent timestep = " << get_parent()->timestep
		 << endl;
	    if (get_parent()->get_nn())
		cerr << "    nn of parent = "
		     << get_parent()->get_nn()->format_label() << endl;
	}

	// Update slow_perturbed lists, if this is a top-level binary.

	hdyn *par = get_parent();
	if (par->is_top_level_node()) {

	    // Set/reset the slow_perturber lists of all perturbers,
	    // to pass consistency checks elsewhere in the code.

	    hdyn *pnode = find_perturber_node();

	    if (!pnode || !pnode->get_valid_perturbers()) {

		// No perturber node?  Shouldn't happen, since this is a
		// prerequisite for slow motion to begin...

		cerr << "startup_slow_motion: "
		     << "no valid perturber node for new slow binary "
		     << get_parent()->format_label() << endl;

	    } else {

		if (diag->slow)
		    cerr << "    perturber node is " << pnode->format_label()
			 << "  n_perturbers = " << pnode->get_n_perturbers()
			 << endl;

		for (int j = 0; j < pnode->get_n_perturbers(); j++) {

		    hdyn *pert_top = pnode->get_perturber_list()[j]
					  ->get_top_level_node();

		    // Only consider perturbers outside the present clump.
		    // (This test isredundant if we limit the slow_perturbed
		    // treatment to top-level binaries only.)

		    if (pert_top != par) {

			slow_perturbed *s = pert_top->find_slow_perturbed(par);
			if (!s)
			    s = pert_top->add_slow_perturbed(par,
						         diag->slow_perturbed);

			if (diag->slow && diag->slow_level > 0)
			    cerr << "    updated slow_perturbed list"
				 << " of top-level node "
				 << pert_top->format_label() << endl;
		    }
		}
	    }
	}
    }
}

static char* term_reason[4] = {"forced",	   // set_stop
			       "invalid plist",	   // get_slowdown_factor = -1
			       "binary perturber", // get_slowdown_factor =  0
			       "perturbed"};	   // get_slowdown_factor =  1

void hdyn::extend_or_end_slow_motion(real P)	// convenient to pass P (!)
						// as an optional argument
{
    int k = -2, old_k = get_kappa();

    if (!slow->get_stop())
	k = get_slowdown_factor(this, P);

    if (k > 1) {

	// Extend the slow motion.

	if (k != old_k) {

	    // Modify the value of kappa.

	    extend_slow(k);

	    // Notify all perturbers that kappa has changed.

	    hdyn *par = get_parent();
	    if (par->is_top_level_node()) {

		hdyn *pnode = find_perturber_node();
		if (pnode && pnode->get_valid_perturbers()) {

		    for (int j = 0; j < pnode->get_n_perturbers(); j++) {

			hdyn *pert_top = pnode->get_perturber_list()[j]
					      ->get_top_level_node();

			if (pert_top != par) {

			    slow_perturbed *s = pert_top
						  ->find_slow_perturbed(par);
			    if (!s)
				s = pert_top->add_slow_perturbed(par,
						       diag->slow_perturbed);

			    s->set_kappa(k);

			    if (diag->slow && diag->slow_level > 0)
				cerr << "    updated slow_perturbed list of"
				     << " top-level node "
				     << pert_top->format_label() << endl;
			}
		    }
		}
	    }

	} else

	    // Accept the current time as the new apocenter time
	    // (already done by extend_slow in other case).

	    slow->set_t_apo(get_time());

	if (diag->slow && (diag->slow_level > 0 || k != old_k))
	    cerr << endl
		 << "extending slow motion for "
		 << get_parent()->format_label()
		 << " at time " << get_time()
		 << "; new kappa = " << get_kappa()
		 << endl;

    } else {

	// End the slow motion: first correct the slow_perturbed lists of
	// all perturbers, if necessary, then delete the slow structure.

	hdyn *par = get_parent();

	if (diag->slow)
	    cerr << endl
		 << "ending slow motion for "
		 << par->format_label()
		 << " at time " << get_time()
		 << " (" << term_reason[k+2] << ")"
		 << endl;

	if (par->is_top_level_node()) {

	    hdyn *pnode = find_perturber_node();

	    if (pnode && pnode->get_valid_perturbers())
		for (int j = 0; j < pnode->n_perturbers; j++) {
		    hdyn *pert_top = pnode->perturber_list[j]
					  ->get_top_level_node();

		    if (pert_top != par) {
#if 0
			cerr << "correcting slow_perturbed list of "
			     << pert_top->format_label() << endl;
#endif
			pert_top->remove_slow_perturbed(par,
							diag->slow_perturbed);
		    }
		}
	}

	delete_slow();
    }
}

void clear_perturbers_slow_perturbed(hdyn * b)
{
    // Clear slow binary b from the slow_perturbed lists of all
    // of its perturbers.

    if (b->get_valid_perturbers())
	for (int j = 0; j < b->get_n_perturbers(); j++) {
	    hdyn * pert_top = b->get_perturber_list()[j]
			       ->get_top_level_node();
#if 0
	    cerr << "removing " << b->format_label();
	    cerr << " from slow_perturbed list of "
		 << pert_top->format_label()
		 << endl;
#endif
	    pert_top->remove_slow_perturbed(b, b->get_kira_diag()
						->slow_perturbed);
	}
}

void list_slow_perturbed(hdyn *b)
{
    // Print slow_perturbed data on all nodes below (and including) b.

    bool found = false;

    for_all_nodes(hdyn, b, bb) {
	slow_perturbed *s = bb->get_sp();
	if (s) {
	    found = true;
	    bb->dump_slow_perturbed("    ");
	}
    }

    if (!found)
	cerr << "    no slow_perturbed data found for "
	     << b->format_label() << endl;
}

void check_slow_consistency(hdyn *b)
{
    // Perform basic checks of slow and slow_perturber data structures.
    //
    //	    1. Check that all slow_perturber lists contain valid entries:
    //
    //		- only top-level nodes should have slow_perturbed lists
    //		- no duplicate entries
    //		- node pointer points to a valid node
    //		- node pointer points to a slow binary node
    //		- node pointer points to a slow binary perturbed
    //		  by 'this' node
    //		- 'this' kappa agrees with node kappa
    //
    //	    2. Check that each perturber of a top-level slow binary
    //	       contains the slow binary on its slow_perturber list.
    //
    // Do not perform corrections if problems are found -- report only.

    // First check slow_perturber lists.

    hdyn *root = b->get_root();

    for_all_nodes(hdyn, b, bi) {
	if (bi != root) {

	    slow_perturbed *s = bi->get_sp();

	    // Check for illegal pointers.

	    if (s && !bi->is_top_level_node())
		cerr << "check_slow_consistency: low-level node "
		     << bi->format_label() << " has a slow_perturbed pointer"
		     << endl;

	    // Check for duplicates.

	    while (s) {
		slow_perturbed *s1 = s->get_next();
		while (s1) {
		    if (s1->get_node() == s->get_node()) {
			cerr << "check_slow_consistency: node "
			     << s->get_node()->format_label()
			     << " (" << s->get_node()
			     << ")" << endl;
			cerr << "                        "
			     << "duplicated on slow_perturbed list of "
			     << bi->format_label() << endl;
		    }
		    s1 = s1->get_next();
		}
		s = s->get_next();
	    }

	    // Then apply basic consistency checks to top-level nodes only.

	    if (bi->is_top_level_node()) {

		s = bi->get_sp();

		while (s) {

		    hdyn *pert_node = (hdyn*)s->get_node();

		    if (!is_valid_slow(pert_node)	// checks valid and slow
			|| !pert_node->is_top_level_node()) {

			cerr << "check_slow_consistency: node "
			     << pert_node->format_label()
			     << " (" << pert_node
			     << ")" << endl;
			cerr << "                        "
			     << "on slow_perturbed list of "
			     << bi->format_label() << endl
			     << "                        "
			     << "is ";
			if (!pert_node->is_valid())
			    cerr << "invalid";
			else if (pert_node->is_top_level_node())
			    cerr << "not a slow binary CM";
			else
			    cerr << "not a top-level binary";
			cerr << " at time " << b->get_system_time() << endl;

		    } else {

			// Pert_node is a valid slow top-level binary.
			// Check that bi (or a component) perturbs it.

			// Note that find_perturber_node starts checking at
			// the parent level...

			hdyn *pnode = pert_node->get_oldest_daughter()
					       ->find_perturber_node();
			if (!pnode) {

			    cerr << "check_slow_consistency: no perturber node"
				 << " for " << pert_node->format_label()
				 << " (" << pert_node << ")"
				 << endl;
			    cerr << "                        "
				 << "on slow_perturbed list of "
				 << bi->format_label()
				 << endl;
			    cerr << "                        "
				 << "at time " << b->get_system_time()
				 << endl;

			} else {

			    // Convention is that binaries are resolved into
			    // components (unless unperturbed) on perturber
			    // lists, but the slow_perturbed list is attached
			    // to the top-level node of the perturber and points
			    // to the center of mass node of the slow binary.

			    bool found = false;
			    for (int k = 0;
				 k < pnode->get_n_perturbers(); k++) {
				hdyn* pk_top = pnode->get_perturber_list()[k]
				    		    ->get_top_level_node();
				if (pk_top == bi) {
				    found = true;
				    break;
				}
			    }

			    if (!found) {

				cerr << "check_slow_consistency: node "
				     << pert_node->format_label()
				     << " (" << pert_node
				     << ") on slow_perturbed list"
				     << endl;
				cerr << "                        "
				     << "of " << bi->format_label()
				     << " is not perturbed by it"
				     << endl;

				pert_node->print_perturber_list();

			    } else {

				// Check kappa.

				real pk = pert_node->get_oldest_daughter()
				    		   ->get_kappa();

				if (pk != s->get_kappa()) {
				    cerr << "check_slow_consistency: node "
					 << pert_node->format_label()
					 << " on slow_perturbed list of ";
				    cerr << bi->format_label()
					 << " has kappa = " << pk
					 << " != " << s->get_kappa()
					 << endl;
				}
			    }
			}
		    }

		    s = s->get_next();

		}
	    }
	}
    }

    // Check consistency in the reverse direction.

    for_all_leaves(hdyn, b, bi) {
	if (bi->is_low_level_node()
	    && bi->get_elder_sister() == NULL
	    && bi->get_slow()) {

	    hdyn *pnode = bi->find_perturber_node();
	    hdyn *par = bi->get_parent();

	    if (!pnode || !pnode->get_valid_perturbers()) {

		cerr << "check_slow_consistency: "
		     << "no valid perturber node for slow binary "
		     << par->format_label() << endl;

	    } else {

		if (par->is_top_level_node()) {

		    for (int k = 0; k < pnode->get_n_perturbers(); k++) {

			hdyn *pert_top = pnode->get_perturber_list()[k]
			    		      ->get_top_level_node();

			if (pert_top != par) {	// redundant, now...

			    // Perturbers not in this clump should have
			    // slow_perturbed top-level data.

			    bool found = false;
			    slow_perturbed *s = pert_top->get_sp();

			    while (s) {
				if (s->get_node() == par) found = true;
				s = s->get_next();
			    }
			    if (!found) {

				cerr << "check_slow_consistency: "
				     << "no reference to slow binary CM "
				     << par->format_label() << endl;
				cerr << "                        "
				     << "on slow_perturbed list of perturber "
				     << pert_top->format_label()
				     << endl;
				cerr << "                        "
				     << "at time " << b->get_system_time()
				     << endl;

				bi->print_pert();		// component MF
				par->print_perturber_list();	// CM MF!

				pert_top->print_slow_perturbed();
			    }
			}
		    }
		}
	    }
	}
    }
}
