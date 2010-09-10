
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //           _\|/_
//=======================================================//             /|\ ~

//
//  hdyn_tree.C: functions related to tree evolution in kira.
//.............................................................................
//    version 1:  Nov 1995   Jun Makino, Steve McMillan
//.............................................................................
//
//  Externally visible functions:
//
//	real   hdyn::distance_to_sister_squared
//	hdyn*  hdyn::new_sister_node
//	int    hdyn::adjust_tree_structure
//
//	void   split_top_level_node
//	void   combine_top_level_nodes
//
//  The functions that actually make changes to the tree structure are
//
//	void split_top_level_node
//	void combine_top_level_nodes
//	void combine_low_level_nodes	(local)
//
//.............................................................................
//
//  NOTE (SLWM, 7/97): when a binary node is created, update_binary_sister
//  is not called, so the components' timesteps, etc. may not be up to date
//  until the next time step is taken.  This is not generally a problem,
//  but may cause difficulties in set_unperturbed_timestep if a lower-level
//  node takes a step before the higher-level binary is advanced.
//
//.............................................................................

#include "hdyn.h"
#include <star/dstar_to_kira.h>
#include <star/single_star.h>


extern bool temp_debug_flag;			// defined in kira.C


// Note:  The *only* places in kira where "close" is defined, for
//	  purposes of making binary trees, are too_close() and
//	  new_sister_node().
//
// Close_criterion is now a variable, part of kira_options:
//
//		0: "close" criterion based on distance
//		1:			      potential
//		2:			      force

// Both criteria 1 and 2 expand the binary distance criterion as masses
// increase.  Both appear to work (Steve, 7/29/98).

real hdyn::distance_to_sister_squared()
{
    vec d_sep = pos - get_binary_sister()->get_pos();
    return d_sep * d_sep;
}

local void init_predictor_tree(hdyn * b)
{
    if (! b->is_root()) {
	b->init_pred();
	b->store_old_force();
    }
    for_all_daughters(hdyn,b,bb) init_predictor_tree(bb);
}

// synchronize_branch:  Synchronize parent and sister nodes up to, but
//                      not including, the specified ancestor (so the
//			ancestor node is synchronized....).

local void synchronize_branch(hdyn * bi, hdyn * ancestor)
{
    if (!is_descendent_of(bi, ancestor, 0))
	err_exit("synchronize_branch: bi not descended from ancestor");

    for (hdyn * bb = bi; bb != ancestor; bb = bb->get_parent()) {

	// Synchronize sister node:

	bb->get_binary_sister()->synchronize_node();

	// Synchronize parent:

	bb->get_parent()->synchronize_node();
    }
}

local inline real max_mass(hdyn* b)
{
    real m_max = 0;
    if (b->is_parent()) {
	for_all_daughters(hdyn, b, bb)
	    if (bb->get_mass() > m_max)
		m_max = Starlab::max(m_max, max_mass(bb));
    } else
	m_max = b->get_mass();

    return m_max;
}

local inline real n_mass(hdyn* b)
{
    // Choose which mass we use in the binary criterion.

    // Must be careful to avoid "polymerization" of clumps containing
    // very massive stars.  Try limiting the mass of a node to the
    // maximum mass of any of its leaves...

    // return b->get_mass();

    return max_mass(b);
}

local bool number_of_daughter_nodes_exceeds(node * b, int n)
{
    int nn = 0;

    for_all_daughters(node, b, bb)
	if (++nn > n)
	    return true;

    return false;
}

local void halve_timestep(hdyn * b)
{
    b->set_timestep(0.5 * b->get_timestep());
//     cerr << "halved timestep for " << b->format_label()
// 	 << " to " << b->get_timestep() << endl;
}

local void print_relative_energy_and_period(hdyn* bi, hdyn* bj)
{
    real E, P;
    get_total_energy_and_period(bi, bj, E, P);

    cerr << "relative E/mu = " << E << flush;

    if (E < 0)
    	cerr << "  P = " << P;				// no endl, note
}


local bool too_close(hdyn * bi, hdyn * bj, real limit_sq,
		     bool check_big = false,
		     bool print = true)
{
    if (bi == bj)
	return false;

    // Note: pred_pos necessary here to take care of runaways where a
    // neighbor time step runs away to zero.  (Steve, 4/03).

    vec ipos = bi->get_pred_pos();
    vec jpos = bj->get_pred_pos();
    vec sep = ipos - jpos;
    real gap_sq = square(sep);

    // Binary criterion modified 7/98 by SLWM and SPZ to depend on
    // acceleration; modified again by SLWM and JM to include potential.
    // The criteria below are chosen to reduce to the simple distance
    // expressions when all masses are equal.  Incorporation of eps2
    // effectively suppresses binaries when eps > d_min.

    bool close;
    real mass_factor = 1;		// default (DISTANCE criterion:
					//		rij  <  d_min)

    if (bi->get_kira_options()->close_criterion > 0) {

	mass_factor = 0.5 * (n_mass(bi) + n_mass(bj)) / bi->get_mbar();

	if (bi->get_kira_options()->close_criterion == 1) {

	    // POTENTIAL criterion:  (mi + mj) / rij  >  2 <m> / d_min
	    //
	    //		==>  rij  <  0.5 * (mi+mj) * d_min / <m>

	    mass_factor *= mass_factor;

	} else if (bi->get_kira_options()->close_criterion == 2) {

	    // FORCE criterion:  (mi + mj) / rij^2  >  2 <m> / d_min^2
	    //
	    //		==>  rij^2  <  0.5 * (mi + mj) * d_min^2 / <m>

	}

    }

    close = (gap_sq + bi->get_eps2() < mass_factor * limit_sq);

#if 1

    // EXPERIMENTAL: Probably shouldn't keep/accept strongly perturbed
    // binaries.  Note that, if we apply this test, we MUST make
    // acceptance (too close) stronger than rejection (too big).

    real pert2 = 0;
    if (check_big) {

	// Call came from too_big() and bi and bj are binary sisters.
	// Use the perturbation as is.  Note that valid_perturbers()
	// is too strong a condition, as a perturbation will be
	// computed even if the list is not valid.  Better simply to
	// look for a non-negative perturbation.  (Steve, 8/07)  Note
	// also a bug in the commented code: bi should be
	// bi->get_parent().

//	if (bi->get_valid_perturbers() && bi->get_n_perturbers() > 0)
	if (bi->get_perturbation_squared() >= 0)
	    pert2 = bi->get_perturbation_squared();

    } else if (close) {

	// bi and bj are close and not binary sisters.  Compute a
	// perturbation from their accelerations.  Write ai = aij +
	// aip, aj = aji + ajp, solve for aip and ajp, then compute
	// the perturbation as in hdyn_ev.C.  Note that, for a binary
	// component, ai = aij + mj*(aip-ajp)/(mi+mj), mi*ai+mj*aj =
	// 0, and perturbation = |mj*(aip-ajp)/(mi+mj)|/|aij|

	real mi = bi->get_mass(), mj = bj->get_mass();
	vec ai = hdyn_something_relative_to_root(bi, &hdyn::get_acc);
	vec aj = hdyn_something_relative_to_root(bj, &hdyn::get_acc);
	real r2eps = gap_sq + bi->get_eps2();
	vec aij = -mj*sep/(r2eps*sqrt(r2eps));
	vec aji = -mi*aij/mj;
	vec aip = ai - aij;
	vec ajp = aj - aji;

	pert2 = square(mj*(aip-ajp)/(mi+mj))/square(aij);
	pert2 *= 1.2;			// harder to combine than to split
    }

    close &= (pert2 < 0.1);		// choice of 0.1 here is ~arbitrary

#endif

    if (close && print && bi->get_kira_diag()->tree
		       && bi->get_kira_diag()->tree_level > 0) {
	cerr << endl << "*** "; bi->print_label(cerr);
	cerr << " too close to "; bj->print_label(cerr);
	cerr << " (criterion " << bi->get_kira_options()->close_criterion
	     << ")" << endl;
    }

    return close;
}

local bool too_big(hdyn * bi, real limit_sq)
{
    hdyn *od = bi->get_oldest_daughter();
    if (!od) return false;

    // Don't split an unperturbed binary even if it is large.

    if (od->get_kepler()) return false;

    // Expand this to extend binaries with small perturbation which are
    // not actually unperturbed.  May need other checks on nn distance,
    // timestep, etc...  (Steve, 12/12/01)

    if (bi->get_valid_perturbers()) {	// was od->... BUG! SMcM 3/08
	real pert2 = od->get_perturbation_squared();
	// if (pert2 >= 0 && pert2 <= od->get_gamma2()) return false;
	if (pert2 >= 0 && pert2 <= 1.e-10) return false;
    }

    bool big = !too_close(od, od->get_younger_sister(), limit_sq,
			  true, false);

    // Don't allow a slow binary to be split, but schedule the slow
    // motion to be stopped at next apocenter.

    if (big && od->get_slow()) {
	if (!od->get_slow()->get_stop()) {
	    cerr << "too_big: scheduling end of slow motion for "
		 << bi->format_label() << " at time " << bi->get_system_time()
		 << endl;
	    od->get_slow()->set_stop();
	}
	return false;
    }

    if (big && bi->get_kira_diag()->tree
	    && bi->get_kira_diag()->tree_level > 0)
	cerr << endl << "*** " << bi->format_label() << " too big" << endl;

    return big;
}

// new_sister_node:  Check to see if a given node is too near any other.
//                   Return a pointer to the node that should become the
//                   new sister, or NULL if the tree is OK as is.
//
//		     **** Used only by adjust_low_level_node.		****
//		     **** Called twice, with 'this' = either component.	****

static int pp3count = 0;

hdyn* hdyn::new_sister_node(bool & top_level_combine)
{
    top_level_combine = false;	// if true, combine top level nodes on return;
    				// if false, combine lower level nodes if a
				//    non-NULL value is returned

    // Note on multiples from Steve, 10/21/98; extended discussion
    // from Steve for Steve added 6/6-13/00.  Some care is needed in
    // dealing with multiples because of a "feature" of function
    // calculate_acc_and_jerk_on_low_level_node (hdyn_ev.C):
    //
    // If (A, B) is the top-level binary of a multiple, and if A is
    // a leaf with nn under B, then both A and B's nn pointers will
    // be correctly set (A to the leaf under B, B to A).  However, if
    // A is a composite node containing B's nn, then A's nn will be
    // B (or a leaf under B), but B's nn will be A (efficiency issue
    // in hdyn_ev).  Note that 'this' in this function may be either
    // A or B, leading (for a triple) to four distinct possibilities.
    // Suppose (a,b) is a binary where b is about to become the sister
    // of a perturber c. Then, depending on the ordering of the
    // particles and which one is 'this', we can have:
    //
    //
    //       O		A = (a,b), B = c, 'this' = A
    //      / \		nn of A is c, nn of c is A
    //     A   c	nn of a is b, nn of b is c
    //    / \							~
    //   a   b		exchange will not be detected until
    // 			component a is advanced (return at *)
    //
    //
    //       O		A = (a,b), B = c, 'this' = c
    //      / \		nn of A is c, nn of c is A
    //     A   c	nn of a is b, nn of b is c
    //    / \							~
    //   a   b		exchange will not be detected until
    // 			component a is advanced (return at *)
    //
    //
    //       O		A = c, B = (a,b), 'this' = c
    //      / \		nn of c is b, nn of B is c
    //     c   B	nn of a is b, nn of b is c
    //        / \						~
    //       a   b	exchange will be detected when c
    // 			is advanced
    //
    //
    //       O		A = c, B = (a,b), 'this' = B
    //      / \		nn of c is b, nn of B is c
    //     c   B	nn of a is b, nn of b is c
    //        / \						~
    //       a   b	exchange will be detected when c
    // 			is advanced
    //
    // Thus, a close encounter between (e.g.) triple members will not be
    // picked up on integration of the outer binary if A is composite,
    // but it will be if A is a leaf.  An exchange in the former case must
    // be detected when A's internal motion is integrated (leaves under A
    // have the correct nn pointers).  This is OK, as there is currently
    // no provision for such changes to be triggered by the parent node.
    // However, proper account must then be taken of the fact that B's nn
    // may not be a leaf, but A.  Failure to do this can lead to runaway
    // time steps (--> 0) when a close encounter goes undetected.
    //
    // Additional note (Steve, 6/14/00).  It is possible that the
    // internal motion never gets the chance to trigger the exchange,
    // as the parent node time step may run away first.  Specific
    // example: A is a receding unbound pair and B happens to pass
    // close to A's center of mass.  We flag this by checking for (a)
    // a strongly perturbed binary with nn lying between its
    // components (i.e. the cause), or (b) disparate center of mass
    // and component time steps in a strongly perturbed binary (the
    // effect).  (See the use of "check_binary_params" below.)  We
    // then synchronize the daughter nodes and reduce their time steps
    // in order to guarantee that they will be on the integration list
    // next time around, and identify the configuration at the
    // component level by including an explicit check (**) on the
    // perturbation within this function.

    if (nn->get_parent() == parent)	// already sisters, so no change needed
	return NULL;			// (* - see note above)

    real mass_factor = 1;		// default (DISTANCE criterion)

    if (options->close_criterion > 0) {

	real m = n_mass(this);
	mass_factor = (m + n_mass(get_binary_sister())) / (m + n_mass(nn));

	if (options->close_criterion == 1) {

	    // POTENTIAL criterion:

	    mass_factor *= mass_factor;

	} else if (options->close_criterion == 2) {

	    // FORCE criterion:

	}
    }

    bool nn_too_close = (distance_to_sister_squared()
				> d_nn_sq * lag_factor * mass_factor);

    if (!nn_too_close) {	// (** Steve 6/14/00)  Relax the distance
				// criterion for a strongly perturbed binary,
				// as described above.  Perhaps better to
				// check explicitly the position of nn
				// relative to the parent node?

	nn_too_close = (perturbation_squared > 1
			&& distance_to_sister_squared() > d_nn_sq);

	if (nn_too_close && diag->tree && diag->tree_level > 1) {
	    cerr << "in new_sister_node for " << format_label()
		 << " at time " << system_time << ":" << endl;
	    cerr << "    candidate new sister for strongly perturbed binary is "
		 << nn->format_label() << endl;
	}
    }

    if (nn_too_close) {
	
	// Nearest neighbor is closer than binary sister.
	//
	// Should check whether or not this and nn are mutual nearest
	// neighbors.  If they are, then we are done.  If not, then we
	// must be a little more careful.  We should further check if
	// one of nn's ancestors can become the new sister.
	//
	//    Criterion:  If the nearest neighbor of any ancestor
	//    of nn is a component of this node, connect these two.

	hdyn * new_sister = nn;
	hdyn * root = get_root();

	// The nn of new_sister may be a node (see note above).  Take care
	// of this possibility before proceeding.  Expect that the nn is
	// the binary sister (should be the only way in which this can
	// occur -- see hdyn_ev), but do the search generally, anyway.

	hdyn *snn = new_sister->get_nn();

	if (snn && snn->is_valid() && !snn->is_leaf()) {

	    // The true nn lies below snn.  Locate it, and modify
	    // new_sister->nn accordingly.  Its distance should
	    // be new_sister->d_nn_sq, but let's not tempt fate...

	    real d_sq_min = VERY_LARGE_NUMBER;
	    hdyn * local_nn = NULL;
	    vec s_pos = hdyn_something_relative_to_root(new_sister,
							   &hdyn::get_pos);

	    for_all_leaves(hdyn, snn, bb) {
		vec b_pos = hdyn_something_relative_to_root(bb,
							       &hdyn::get_pos);
		real d_sq = square(b_pos - s_pos);
		if (d_sq < d_sq_min) {
		    local_nn = bb;
		    d_sq_min = d_sq;
		}
	    }

	    if (local_nn) {

		snn = local_nn;
		new_sister->set_nn(snn);

		if (diag->tree && diag->tree_level > 0) {
		    cerr << "new_sister_node:  computed local nn for "
			 << new_sister->format_label() << endl;
		    PR(local_nn);
		    cerr << "  (" << snn->format_label() << ")" << endl;
		    PRC(new_sister->get_d_nn_sq()); PRL(d_sq_min);
		}
	    }
	}

	// Now snn is the leaf closest to new_sister and is new_sister->nn.
	// Attach to the highest possible ancestor of new_sister.  If snn
	// is this (or part of this), then we are done.  Otherwise, move
	// up the tree until new_sister->nn is part of or contains this.

	if (snn && snn->is_valid() && common_ancestor(snn, this) != this) {

	    hdyn* ns = new_sister;

	    do {
		new_sister = new_sister->get_parent();

	    } while(new_sister != root
		    && new_sister->get_nn() != NULL
		    && new_sister->get_nn()->is_valid()
		    && new_sister->parent != parent
		    && common_ancestor(new_sister->get_nn(), this) != this
		    && common_ancestor(new_sister->get_nn(), this)
				!= new_sister->get_nn());

	    // (Last line above added by Steve 9/04.  The logic here
	    //  is getting out of hand...)

	    if (new_sister == root
		|| new_sister->get_nn() == NULL
		|| !new_sister->get_nn()->is_valid()
		|| new_sister->parent == parent
		|| common_ancestor(new_sister, this)== this
		|| common_ancestor(new_sister, this)== new_sister
		) {

		return NULL;
	    }

	    if (diag->tree && diag->tree_level > 0) {
		cerr << "new_sister_node:  ascended tree from "
		     << ns->format_label();
		cerr << " to " << new_sister->format_label() << endl;
	    }
	}

	hdyn * ancestor = common_ancestor(this, new_sister);

	if (ancestor->is_root()) {

	    top_level_combine = TRUE;
	    return new_sister;

	} else {

	    // The new sister node is in the same subtree.
	    // Attach this node to the location as high up as possible.

	    if (parent != ancestor) {
		while (new_sister->get_parent() != ancestor) {
		    new_sister = new_sister->get_parent();
		}
	    }

	    if (new_sister->kep == NULL) {
		if (diag->tree && diag->tree_level > 0) {
		    cerr << "new sister node for "; pretty_print_node(cerr);
		    cerr << " (nn = "; nn->pretty_print_node(cerr);
		    cerr << ") is " ; new_sister->pretty_print_node(cerr);
		    PRI(2); PRL(top_level_combine);
		}
		return new_sister;
	    }
	}
    }
    return NULL;
}


local void check_merge_esc_flags(hdyn *bi, hdyn *bj)
{
    // Update any escaper flags in two nodes being combined.
    // The CM node is flagged as escaping iff both components are.

    bool iesc = find_qmatch(bi->get_dyn_story(), "t_esc");
    bool jesc = find_qmatch(bj->get_dyn_story(), "t_esc");

    if (iesc && jesc) {

	// Both components are flagged as escapers.

	real t_esc = Starlab::max(getrq(bi->get_dyn_story(), "t_esc"),
			          getrq(bj->get_dyn_story(), "t_esc"));
	putrq(bi->get_parent()->get_dyn_story(), "esc", 1);
	putrq(bi->get_parent()->get_dyn_story(), "t_esc", t_esc);

    } else {

	// One component wasn't flagged as escaping.

	putrq(bi->get_parent()->get_dyn_story(), "esc", 0);

	if (iesc) {
	    putrq(bi->get_dyn_story(), "esc", 0);
	    rmq(bi->get_dyn_story(), "t_esc");
	}

	if (jesc) {
	    putrq(bj->get_dyn_story(), "esc", 0);
	    rmq(bj->get_dyn_story(), "t_esc");
	}
    }
}

void combine_top_level_nodes(hdyn * bj, hdyn * bi,
			     int full_dump)		// default = 0
{
    // Combine two top-level nodes into a binary.  Original node names
    // and addresses are retained, so neighboring perturber lists are
    // unaffected.  bi is the node whose step caused the change; bj is
    // its nearest neighbor.  The new CM node will be (bj, bi).

    if (bj->get_kira_diag()->tree) {

        cerr << endl << "combine_top_level_nodes: attaching ",
	bj->pretty_print_node(cerr);
	cerr << " to ", bi->pretty_print_node(cerr);
	int p = cerr.precision(INT_PRECISION);
	cerr << " at time " << bj->get_system_time();
	cerr.precision(p);

	if (bj->get_kira_diag()->tree_level > 0) {

	    cerr << endl;
	    print_binary_from_dyn_pair(bj, bi);
	    cerr << endl;
	    if (bj->is_parent()) {
		print_binary_from_dyn_pair(bj->get_oldest_daughter(),
					   bj->get_oldest_daughter()
					     ->get_younger_sister());
		cerr << endl;
	    }
	    if (bi->is_parent()) {
		print_binary_from_dyn_pair(bi->get_oldest_daughter(),
					   bi->get_oldest_daughter()
					     ->get_younger_sister());
		cerr << endl;
	    }

	    if (bj->get_kira_diag()->tree_level > 1) {

		// Not really necessary if pp3 (new node) occurs below:
		//
		// cerr << endl;
		// pp3(bi, cerr);
		// pp3(bj, cerr);

		if (bj->get_kira_diag()->tree_level > 2) {
		    cerr << endl;
		    put_node(bj, cerr, bj->get_kira_options()->print_xreal);
		    put_node(bi, cerr, bi->get_kira_options()->print_xreal);
		}

		cerr << endl;
		// plot_stars(bj);	// *** may cause problems in code
					// *** as of 8/03 (maybe bug in
					// *** perturber lists...)
	    }

	} else {

	    cerr << endl << "                         ";
	    print_relative_energy_and_period(bj, bi);
	    cerr << endl;
	}
    }

    // Make sure bi and bj are up to date.

    predict_loworder_all(bi->get_root(), bi->get_system_time());

    // print_recalculated_energies(bi->get_root());
    // pp3((bi->get_root()), cerr);
    // synchronize_tree(bi->get_root(), time);
    // print_recalculated_energies(bi->get_root());
    // pp3((bi->get_root()), cerr);

    bi->synchronize_node();
    bj->synchronize_node();

    if (full_dump) {

	// Components are synchronized, but lower levels are not.
	// The 4tree software needs trees to be synchronous, but
	// synchronizing everything here will likely lead to runaway
	// timesteps.  Function put_node in this case must use
	// *predicted* quantities and system time as time.
	// All necessary predictions have already been done.)

	// Dump out the "before" system (bi and bj), for use in 4tree
	// applications, with "defunct" flags attached (final "3").

	// cerr << "combine_top_level_nodes: time " << bi->get_system_time();
	// cerr << "  put_node for " << endl << "    " << bi->format_label();
	// cerr << " and " << bj->format_label() << endl;

	put_node(bi, cout, false, 3);
	put_node(bj, cout, false, 3);
    }

    // Remove any slow_perturbed references to the components from
    // other top-level nodes.  Only possible if a component is a slow
    // binary...

    for (hdyn *bb = bi; bb != bj; bb = bj)
	if (bb->get_kappa() > 1) {
	    hdyn *pnode = bb->find_perturber_node();
	    if (pnode && pnode->get_valid_perturbers())
		for (int k = 0; k < pnode->get_n_perturbers(); k++) {
		    hdyn *p = pnode->get_perturber_list()[k]
				   ->get_top_level_node();
		    if (p->get_sp())
			p->remove_slow_perturbed(bb);		
		}
	}

    create_binary_from_toplevel_nodes(bi, bj);		// --> (bj, bi)

    hdyn *par = bj->get_parent();
//    pp3(par);

    // Copy any slow_perturbed lists to the new top-level node, and
    // delete the low-level lists.

    bi->copy_slow_perturbed(par);
    bi->clear_slow_perturbed();
    bj->copy_slow_perturbed(par);
    bj->clear_slow_perturbed();

    // Delete any slow_perturbed element in the new CM node that refers
    // to either component.

    if (par->get_sp()) {
	par->remove_slow_perturbed(bi);
	par->remove_slow_perturbed(bj);
    }

    // Reset all perturber lists below the newly-formed top-level
    // node -- probably not strictly necessary for bi and bj if
    // low-level perturber lists are allowed, but cleaner and much
    // more convenient to clear and recompute the lists.  Note that 
    // this will result in perturbations being computed using the
    // entire system until the top-level list is computed.

    for_all_nodes(hdyn, par, bb) {
	bb->set_valid_perturbers(false);
	bb->set_perturbation_squared(-1);
    }

    // Parent step was set equal to the minimum daughter step in
    // function create_binary_from_toplevel_nodes(). Old code reduced
    // it here to make it less than the daughter steps, but consistent
    // with the current time, so that the parent perturber list will
    // be computed before the internal motion is advanced.  May not be
    // necessary, since kira integrates top-level nodes before
    // components.
    //
    // Note that this almost certainly starts the CM off with a step
    // far shorter than necessary.  Assume that it is OK to let the
    // timestep algorithm double the step back to the proper value
    // over the next 10 or so time steps.
    //
    // An alternative would be to compute the perturber lists
    // etc. explicitly, but we can't do this here, as the GRAPE
    // doesn't yet know about the CM.  Could do it from the top-level
    // in kira.  (Steve, 3/03)
    //
    // A bad side effect of halving here is that if the pair is
    // immediately split up again, the effect is to reduce all time
    // steps by a factor of two.  This can cause zero-length steps,
    // but should only occur if comething is wrong with the
    // combine/split logic...


    // halve_timestep(par);	// *** commented out by Steve 9/8/10 ***

    real dt_par = par->get_timestep();
    // PRL(dt_par);

    // Note from Steve (8/03): for multiples, really want to have the
    // CM step precede any low-level step, not just the daughters.
    // May not be convenient if a daughter step is very short, or just
    // happens to place the daughter in a bad timestep block.

    // SPZ: changed real min_time to xreal min_time
    xreal min_time = VERY_LARGE_NUMBER;
    for_all_nodes(hdyn, bj, bb)
        if (bb != bj && !bb->get_elder_sister()) {
	    real bb_time = bb->get_next_time();
	    if (bb_time < min_time) min_time = bb_time;
	}
    for_all_nodes(hdyn, bi, bb)
        if (bb != bi && !bb->get_elder_sister()) {
	    real bb_time = bb->get_next_time();
	    if (bb_time < min_time) min_time = bb_time;
	}

    real dt_min = min_time - bj->get_system_time();
    // PRL(min_time);
    if (dt_min < dt_par) {
      if (dt_min > dt_par/16) {			// 16 is arbitrary
	  dt_par = dt_min;
	  par->set_timestep(dt_par);
	  cerr << "reduced timestep for " << par->format_label()
	       << " to " << dt_par << endl;

      }
      // else cerr << "combine_top_level_nodes: "
      //           << "Can't reduce parent timestep sufficiently..." << endl;
    }

//    PRC(dt_par); PRC(dt_min); PRL(par->get_timestep());
//    pp3(par);

    bi->init_pred();
    bj->init_pred();
    bj->get_parent()->init_pred();

    init_predictor_tree(bj->get_parent());	// for diagnostic output...

    check_merge_esc_flags(bj, bi);

    if (full_dump) {

	// Dump out the "after" system (parent), for use in 4tree
	// applications.  Must predict this portion of the tree again.

	predict_loworder_all(bj->get_parent(), bj->get_system_time());

	// cerr << "combine_top_level_nodes: time " << bi->get_system_time();
	// cerr << "  put_node for " << bj->get_parent()->format_label()
	//      << endl;

	put_node(bj->get_parent(), cout, false, 2);
    }

    bi->get_kira_counters()->top_level_combine++;
}

void split_top_level_node(hdyn * bi,
			  int full_dump)	// default = 0
{
    // Split a top-level node into two components.  The original node is
    // deleted, so any references to it (in nn/coll pointers or perturber
    // lists) must be corrected.

    if (!bi->is_top_level_node())
	err_exit("split_top_level_node: not at top level");

    if (bi->get_oldest_daughter()->get_slow()) {

	// Shouldn't happen -- flag and defer the split.

	cerr << "split_top_level_node: warning: trying to split slow binary "
	     << bi->format_label() << endl
	     << "    at time " << bi->get_system_time() << endl;

	// bi->get_oldest_daughter()->delete_slow();

	bi->get_oldest_daughter()->get_slow()->set_stop();
	return;
    }

    if (bi->get_kira_diag()->tree) {

        cerr << endl << "split_top_level_node: splitting ";
	bi->pretty_print_node(cerr);
	int p = cerr.precision(INT_PRECISION);
	cerr << " at time " << bi->get_system_time();
	cerr.precision(p);

	if (bi->get_kira_diag()->tree_level >= 1) {
	    cerr << endl;
	    pp2(bi, cerr);
	    print_binary_from_dyn_pair(bi->get_oldest_daughter(),
				       bi->get_oldest_daughter()
				         ->get_younger_sister());
	    cerr << endl << flush;

	    if (bi->get_kira_diag()->tree_level > 1) {
		cerr << endl;
		pp3(bi, cerr);

		// Excessive output?

		if (bi->get_kira_diag()->tree_level > 2)
		    put_node(bi, cerr, bi->get_kira_options()->print_xreal);
	    }

	} else {

	    cerr << endl << "                      ";
	    print_relative_energy_and_period(bi->get_oldest_daughter(),
					     bi->get_oldest_daughter()
					       ->get_younger_sister());
	    cerr << endl << flush;
	}

	// Info on perturbers:

	cerr << "                     ";
	hdyn *pnode = bi->get_oldest_daughter()->find_perturber_node();
	if (pnode && pnode->is_valid()) {
	    cerr << "pnode = " << pnode->format_label();
	    if (pnode->get_valid_perturbers()) {

		// Test here for pert2 > 0.  It shouldn't be negative, and
		// shouldn't be a problem even if it is, but this has been
		// found to put at least one compiler (gcc 3.4.2) on at
		// least one OS (RedHat Fedora core 3) into an infinite
		// loop!

		real pert2 = bi->get_oldest_daughter()
			       ->get_perturbation_squared();
		if (pert2 > 0) pert2 = sqrt(pert2);

	        cerr << ", " << pnode->get_n_perturbers()
		     << " perturber(s) (" << pert2 << ")";
	    } else
	        cerr << ", perturbers unknown";
	} else
	    cerr << "pnode invalid";

	    cerr << endl << flush;
    }

    hdyn *od = bi->get_oldest_daughter();
    hdyn *yd = od->get_younger_sister();

    // Synchronize the nodes involved.

    predict_loworder_all(bi->get_root(), bi->get_system_time());

    // for debugging ...
    // synchronize_tree(bi->get_root());
    // synchronize_tree(bi);

    bi->synchronize_node();
    od->synchronize_node();
    update_binary_sister(od);

    if (full_dump) {

	// Components are synchronized, but lower levels are not.

	// Dump out the "before" system, for use in 4tree applications,
	// using predicted quantities.

	// cerr << "split_top_level_node: time " << bi->get_system_time();
	// cerr << "  put_node for " << bi->format_label() << endl;

	put_node(bi, cout, false, 3);
    }

#if 0

    // Compute the CM and component accelerations directly, using the CM
    // approximation for other stars; then repeat component calculation
    // using perturbers.

    vec xod = bi->get_pos() + od->get_nopred_pos();
    vec xyd = bi->get_pos() + yd->get_nopred_pos();
    vec acm = 0, aod = 0, ayd = 0;
    for_all_daughters(hdyn, od->get_root(), bb) {
      if (bb != bi) {
	vec x = bb->get_nopred_pos() - bi->get_pos();
	real r2 = x*x + bi->get_eps2();
	acm += bb->get_mass()*x/(r2*sqrt(r2));
	x = bb->get_nopred_pos() - xod;
	r2 = x*x + od->get_eps2();
	aod += bb->get_mass()*x/(r2*sqrt(r2));
	x = bb->get_nopred_pos() - xyd;
	r2 = x*x + yd->get_eps2();
	ayd += bb->get_mass()*x/(r2*sqrt(r2));
      }
    }
    vec atop = (od->get_mass()*aod + yd->get_mass()*ayd)/bi->get_mass();
    vec x12 = yd->get_nopred_pos() - od->get_nopred_pos();
    real r2 = x12*x12 + bi->get_eps2();
    vec a12 = yd->get_mass()*x12/(r2*sqrt(r2));
    aod += a12;
    ayd -= a12*od->get_mass()/yd->get_mass();

    vec aodp = 0;
    int np = bi->get_n_perturbers();
    hdyn **plist = bi->get_perturber_list();
    cerr << "perturbers: ";
    for (int k = 0; k < np; k++) {
        hdyn *p = plist[k];
	cerr << p->format_label() << " ";
	vec xp = hdyn_something_relative_to_root(p, &hdyn::get_nopred_pos);
	vec x = xp - xod;
	r2 = x*x + od->get_eps2();
	aodp += p->get_mass()*x/(r2*sqrt(r2));
	x = xp - xyd;
	r2 = x*x + yd->get_eps2();
	aodp -= p->get_mass()*x/(r2*sqrt(r2));
    }
    cerr << endl;
    aodp *= yd->get_mass()/bi->get_mass();
    aodp += a12;

    PRL(bi->get_acc());
    PRL(acm);
    PRL(atop);
    PRL(atop-bi->get_acc());

    PRL(od->get_acc());
    PRL(aod-atop);
    PRL(aod-atop-od->get_acc());

    PRL(yd->get_acc());
    PRL(ayd-atop);
    PRL(ayd-atop-yd->get_acc());

    PRL(a12);
    PRL(np);
    PRL(aodp);
    PRL(aodp-od->get_acc());

#endif

    // Express od and yd quantities relative to root.

    od->set_pos(hdyn_something_relative_to_root(od, &hdyn::get_pos));
    od->set_vel(hdyn_something_relative_to_root(od, &hdyn::get_vel));
    yd->set_pos(hdyn_something_relative_to_root(yd, &hdyn::get_pos));
    yd->set_vel(hdyn_something_relative_to_root(yd, &hdyn::get_vel));

    od->set_acc(hdyn_something_relative_to_root(od, &hdyn::get_acc));
    od->set_jerk(hdyn_something_relative_to_root(od, &hdyn::get_jerk));
    yd->set_acc(hdyn_something_relative_to_root(yd, &hdyn::get_acc));
    yd->set_jerk(hdyn_something_relative_to_root(yd, &hdyn::get_jerk));

    // Reconstructing acc and jerk this way explicitly neglects the
    // differential contribution of non-perturbers, introducing errors
    // at the level of the perturbation threshold.  Ordinarily this is
    // negligible, and is recovered after a few steps.  However, in
    // exceptional circumstances (or due to bugs elsewhere), it can
    // happen that we repeatedly form and destroy a binary, in which
    // case repeatedly introcuce this (constant) error and may force
    // the time step to zero.  Safest to recompute acc and jerk on od
    // and yd once they are restored below.  (SLWM, 5/07)

    od->store_old_force();
    yd->store_old_force();

    od->init_pred();
    yd->init_pred();

    // May not be necessary to do this if low-level perturbers are
    // allowed, but cleaner to recompute the lists.  Clear all
    // perturber lists below the new top-level nodes.

    for_all_nodes(hdyn, od, bb) {
	bb->set_valid_perturbers(false);
	bb->set_perturbation_squared(-1);
    }
    for_all_nodes(hdyn, yd, bb) {
	bb->set_valid_perturbers(false);
	bb->set_perturbation_squared(-1);
    }

    // Dissolve the parent node and attach od and yd at top level.

    hdyn *root = bi->get_root();

    detach_node_from_general_tree(bi);

    // Imperfect correction of neighbor pointers:

    if (bi->get_nn()->get_nn() == bi)
	bi->get_nn()->set_nn(NULL);

    // Should also check for the unlikely possibility that:
    //
    //	(1) binary bi is unperturbed, and
    //	(2) binary bi is on the perturber list of some other binary.
    //
    // If so, bi must be removed from the perturber list before being
    // deleted.

    if (!RESOLVE_UNPERTURBED_PERTURBERS)
	expand_perturber_lists(bi->get_root(), bi,
			       bi->get_kira_diag()->tree
			           && bi->get_kira_diag()->tree_level > 0);

    // Copy any slow_perturbed list to the components.

    bi->copy_slow_perturbed(od, true);		// "true" ==> overwrite
    bi->copy_slow_perturbed(yd, true);

    delete bi;

    add_node(od, root);
    add_node(yd, root);

#if 1

    //    if (od->get_system_time() > 361.15673100948) {
    if (od->get_system_time() > 0.) {

	// Recompute acc and jerk on od and yd (see note above).  Should
	// be possible to use GRAPE to do this, but defer for now to avoid
	// possible internal bookkeeping problems in the kira_ev/GRAPE
	// code -- to be checked.  (Steve, 5/07)

        // Note that the front-end calculation must deviate from that
        // done on the GRAPE at some level of precision.  In that
        // case, there will always be a discrepancy in the
        // accelerations at some level, which will eventually
        // translate into the time step going to zero.  Only solution
        // may be to avoid this circumstance by being more careful
        // about combine/split commands...

	od->calculate_acc_and_jerk(true);
	yd->calculate_acc_and_jerk(true);

	od->store_old_force();
	yd->store_old_force();

	od->init_pred();
	yd->init_pred();
    }

#endif

#if 0
      temp_debug_flag = true;
      cerr << "setting temp_debug_flag" << endl;
#endif

    // Bad idea to routinely reduce time steps.

    // halve_timestep(od);
    // halve_timestep(yd);

    if (full_dump) {

	// Dump out the "after" system (od and yd), for use in 4tree
	// applications.

	predict_loworder_all(od, od->get_system_time());
	predict_loworder_all(yd, yd->get_system_time());

	// cerr << "split_top_level_node: time " << od->get_system_time();
	// cerr << "  put_node for " << endl << "    " << od->format_label();
	// cerr << " and " << yd->format_label() << endl;

	put_node(od, cout, false, 2);
	put_node(yd, cout, false, 2);
    }

    bi->get_kira_counters()->top_level_split++;
}

local void combine_low_level_nodes(hdyn * bi, hdyn * bj,
				   int full_dump = 0)
{
    // Rearrange the tree structure within a multiple system to try
    // to make bi and bj sister nodes.

    if (bi->get_slow()) {

	// Not sure if this can happen, but flag it just in case...
	// Hmmm, looks like it can happen (Steve, 9/00).

	// Old code:
	//
	// cerr << "combine_low_level_nodes: warning: splitting slow binary "
	//      << bi->get_parent()->format_label() << " at time "
	//      << bi->get_system_time()
	//      << endl;
	// bi->delete_slow();	// will likely cause a crash later!

	// New:

	// As in split_top_level_node, defer the split
	// and schedule the slow motion for termination...

	bi->get_slow()->set_stop();
	return;
    }

    bool pp3_at_end = false;
//    pp3_at_end = (bi->get_time() > 82 && bi->get_time() < 82.5);

    if (bi->get_kira_diag()->tree && bi->get_kira_diag()->tree_level > 0) {
	cerr << endl << "combine_low_level_nodes:  combining ";
	bi->pretty_print_node(cerr);
	cerr << " and ";
	bj->pretty_print_node(cerr);
	int p = cerr.precision(HIGH_PRECISION);
	cerr << " at time " << bi->get_system_time() << "\n";
	cerr.precision(p);
//	pp3_at_end = true;
    }

    hdyn *ancestor;
    ancestor = common_ancestor(bi, bj);

#if 0
    if (bi->get_time() > 2.6 && bi->get_time() < 2.7) {
	cerr << endl << "------------------------------" << endl;
	cerr << "combine_low_level_nodes:" << endl;
	cerr << "bi = " << bi->format_label() << endl;
	cerr << "bj = " << bj->format_label() << endl;
	cerr << "ancestor = " << ancestor->format_label() << endl;
	pp3(ancestor);
	pp3_at_end = true;
    }
#endif

    // Make sure bi and bj are up to date.

    predict_loworder_all(bi->get_root(), bi->get_system_time());

    bi->synchronize_node();
    bj->synchronize_node();

    synchronize_branch(bi, ancestor);
    synchronize_branch(bj, ancestor);

    hdyn* old_top_level_node = bi->get_top_level_node();

//    if (full_dump == 1) {
    if (full_dump) {

	// Dump out the "before" system (top-level), for use in 4tree
	// applications.

	// cerr << "combine_low_level_nodes: time " << bi->get_system_time();
	// cerr << "  put_node for " << old_top_level_node->format_label()
	//      << endl;

	put_node(old_top_level_node, cout, false, 3);
    }

    // Note from Steve, 3/99:
    //
    // Function move_node will remove the parent node of bi and
    // reattach node bi as a sister of bj.  If the binary sister
    // of bi is also a binary and contains bj, then it will become
    // an ancestor node of the new system.  If, as is likely, bi is
    // a perturber of its binary sister, the perturber list of the
    // new system will then be corrupted, as the ancestor node will
    // still contain bi as a perturber.  Avoid this unpleasantness
    // by forcing the perturber lists of the entire subtree to be
    // recomputed.  However, retain the top-level perturber list,
    // which should be unchanged, so we at least don't have to
    // compute perturbations using the entire system.
    // Messy, messy...
    //
    // Messier still, if bi's parent node happens to be the top-level
    // node, then the top-level node will be deleted and its perturber
    // information lost.  Save the perturber info and restore it to
    // the eventual top-level node.

    bool vp = old_top_level_node->get_valid_perturbers();
    if (!vp) vp = old_top_level_node->get_valid_perturbers_low();
    int np = 0;
    real p_sq = -1;
    real p_fac = 0;
    hdyn** pert_list = NULL;

    if (vp) {
	if (old_top_level_node->get_valid_perturbers())
	    np = old_top_level_node->get_n_perturbers();
	else
	    np = old_top_level_node->get_n_perturbers_low();
	p_sq = old_top_level_node->get_perturbation_squared();
	p_fac = old_top_level_node->get_perturbation_radius_factor();
	pert_list = new hdyn *[MAX_PERTURBERS];
	for (int i = 0; i < np; i++)
	    pert_list[i] = old_top_level_node->get_perturber_list()[i];
	// cerr << "saved top-level perturber information" << endl;
    }

    move_node(bi, bj);
    hdyn *top_level_node = bi->get_top_level_node();

    if (top_level_node != old_top_level_node) {

	if (vp) {

	    // Restore original top-level perturbation information.
	    // (Might have been easier simply to copy the node...)

	    top_level_node->remove_perturber_list();
	    top_level_node->new_perturber_list();

	    top_level_node->set_valid_perturbers(old_top_level_node
						 ->get_valid_perturbers());
	    top_level_node->set_n_perturbers(old_top_level_node
					     ->get_n_perturbers());
	    top_level_node
		->set_valid_perturbers_low(old_top_level_node
					   ->get_valid_perturbers_low());
	    top_level_node
		->set_n_perturbers_low(old_top_level_node
				       ->get_n_perturbers_low());
	    top_level_node->set_perturbation_squared(p_sq);
	    top_level_node->set_perturbation_radius_factor(p_fac);

	    for (int i = 0; i < np; i++)
		top_level_node->get_perturber_list()[i] = pert_list[i];

	    delete pert_list;
	    cerr << "restored top-level perturber information, np = "
		 << np << endl;

	} else {

	    top_level_node->set_valid_perturbers(vp);
	    top_level_node->set_perturbation_squared(-1);

	}
    }	

    // halve_timestep(bi);
    // halve_timestep(bj->get_parent());

    predict_loworder_all(bi->get_root(), bi->get_system_time());

    // Could force all low-level perturber lists in the present subtree 
    // to be rebuilt by setting valid_perturbers false, but better to
    // rebuild now...

    // "False" here means rebuild all low-level lists, even if non-null.

    if (ALLOW_LOW_LEVEL_PERTURBERS)
	top_level_node->create_low_level_perturber_lists(false);
    
    bi->get_kira_counters()->low_level_combine++;

    // Make sure that all CM node names are up to date.

    hdyn *bn = bi->get_parent();
    while (bn != bi->get_root()) {
	label_binary_node(bn);
	bn = bn->get_parent();
    }

    check_merge_esc_flags(bi, bi->get_binary_sister());

//    if (full_dump == 1) {
    if (full_dump) {

	// Dump out the "after" system (new top-level), for use in 4tree
	// applications.

	predict_loworder_all(top_level_node, bi->get_system_time());

	// cerr << "combine_low_level_nodes: time " << bi->get_system_time();
	// cerr << "  put_node for " << top_level_node->format_label() << endl;

	put_node(top_level_node, cout, false, 2);
    }

    if (pp3_at_end) {
	cerr << endl << "after combine_low_level_nodes:" << endl;
	pp3(bi->get_top_level_node());
	cerr << endl << "------------------------------" << endl;
    }
}

local int adjust_low_level_node(hdyn * bi, int full_dump = 0)
{
    int status = 1;

    bool top_level_combine;	// TRUE iff top-level nodes are to be combined.
				// Takes precedence over the returned value of
				// new_sister_node.  If new_sister_node is
				// non-NULL, it is the node which should become
				// the new sister of node bi within the same
				// binary subtree.

    hdyn *sister = bi->new_sister_node(top_level_combine);

#if 0
    if (sister != NULL) {
	cout << "adjust low_level, "; bi->pretty_print_node(cout);
	cout << "-- "; sister->pretty_print_node(cout);
	cout << "-- "; sister->get_nn()->pretty_print_node(cout); cout<<endl;
    }
#endif

    if (sister != NULL)
	predict_loworder_all(bi->get_root(), bi->get_system_time());

    if (top_level_combine) {

	if (number_of_daughter_nodes_exceeds(bi->get_root(), 2)) {

	    if (bi->get_kira_diag()->tree
		    && bi->get_kira_diag()->tree_level > 1) {
		cerr << "\nadjust_low_level_node: "
		     << "normal top level combine at time "
		     << bi->get_time() << endl;
		// print_recalculated_energies(bi->get_root());
	    }

	    // cerr<< "Call combine_top_level_nodes from adjust_low_level_node"
	    //     << endl;

	    combine_top_level_nodes(sister->get_top_level_node(),
				    bi->get_top_level_node(), full_dump);
	    status = 2;

	} else {

	    if (bi->get_kira_diag()->tree
		    && bi->get_kira_diag()->tree_level > 1) {
		cerr << "adjust_low_level_node:"
		     << " top level combine with two top level nodes at time "
		     << bi->get_time() << endl;
		// print_recalculated_energies(bi->get_root());
	    }

	    hdyn* t = bi->get_top_level_node();
	    hdyn* od = t->get_oldest_daughter();

	    if (!od)
		status = 0;		// shouldn't happen
	    else {

		if (od->get_slow()) {

		    if (!od->get_slow()->get_stop()) {
			cerr << "adjust_low_level_node: scheduling end "
			     << "of slow motion for " << bi->format_label()
			     << " at time " << bi->get_system_time()
			     << endl;
			od->get_slow()->set_stop();
		    }
		    status = 0;

		} else {

		    // cerr << "call split_top_level_node 1" << endl;

		    split_top_level_node(t, full_dump);
		    status = 3;

		}
	    }
	}

    } else if (sister != NULL) {

	combine_low_level_nodes(bi, sister, full_dump);

    } else

	status = 0;

    return status;
}


local inline bool check_binary_params(hdyn *b)
{
    if (b->is_low_level_node()) {

	hdyn *od = b->get_oldest_daughter();

	if (od						// b is a low-level,
	    && od->get_perturbation_squared() > 1) {	// strongly perturbed,
	    						// binary node

	    bool sync = false;
	    int reason = 0;

	    if (od->get_timestep()			// very disparate
		    > 128*b->get_timestep())		// component timesteps
							// (128 is arbitrary)
		sync = true;

	    else if (b->distance_to_sister_squared()	// components are too
	    	     < od->distance_to_sister_squared()	// far apart, with some
		    && od->get_timestep()		// timestep disparity
			> 4*b->get_timestep()) {	// (4 is arbitrary)

		sync = true;
		reason = 1;

	    }

	    // Synchronize the binary components if they are not on the
	    // current integration list, and reset the binary timestep
	    // to match the parent node.

	    if (od->get_time() == b->get_time()) sync = false;

	    if (sync) {

		od->synchronize_node();
		od->set_timestep(b->get_timestep());
		cerr << "set component timestep for " << od->format_label()
		     << " to " << b->get_timestep() << endl;


		if (b->get_kira_diag()->tree
		    && b->get_kira_diag()->tree_level > reason)
		    cerr << "check_binary_params: synchronizing ("
			 << reason << ") "
			 << b->format_label() << " components at time "
			 << b->get_system_time() << endl;

		return true;
	    }
	}
    }

    return false;
}

local void print_tree_structure(hdyn *bb)		// for debugging
{
    hdyn *b = bb->get_top_level_node();
    hdyn *od = b->get_oldest_daughter();

    if (od) {
	hdyn *yd = od->get_younger_sister();

	cerr << endl << "tree structure for " << b->format_label()
	     << " at time " << b->get_system_time() << endl;
	cerr << "output triggered by " << bb->format_label() << endl;

	// Daughters:

	cerr << "daughters are " << od->format_label() << " and ";
	cerr << yd->format_label();
	cerr << ", separation = " << abs(od->get_pos()-yd->get_pos()) << endl;
	cerr << "daughter neighbors are " << od->get_nn()->format_label()
	     << " and ";
	cerr << yd->get_nn()->format_label() << endl;

	// Granddaughters:

	for_all_daughters(hdyn, b, d) {			// d is od or yd
	    hdyn *od1 = d->get_oldest_daughter();
	    if (od1) {
		cerr << "    daughters of " << d->format_label();
		hdyn *od2 = od1->get_younger_sister();
		cerr << " (separation = "
		     << abs(od1->get_pos()-od2->get_pos()) << "):" << endl;
		hdyn *s = d->get_binary_sister();
		for_all_daughters(hdyn, d, dd) {	// dd is od1 or od2
		    if (s->is_parent()) {
			for_all_daughters(hdyn, s, ss) {
			    cerr << "        distance from "
				 << dd->format_label();
			    cerr << " to " << ss->format_label() << " is ";
			    cerr << abs(d->get_pos()+dd->get_pos()
					-s->get_pos()-ss->get_pos())
				 << endl;
			}
		    } else {
			cerr << "        distance from " << dd->format_label();
			cerr << " to " << s->format_label() << " is ";
			cerr << abs(d->get_pos()+dd->get_pos()-s->get_pos())
			     << endl;
		    }
		}
	    }
	}
    }	
}

int hdyn::adjust_tree_structure(int full_dump)		// default = 0
{
    int status = 0;

    // Value of status (return value) indicates adjustment that occurred:
    //
    //		0	no adjustment
    //		1	low-level combine
    //		2	top-level combine (direct)
    //		3	top-level split (direct)
    //		4	top-level combine (indirect)
    //		5	top-level split (indirect)
    //		6	low-level timestep change (no tree change)
    //
    // (Similar coding also applies to adjust_low_level_node.)

    // cerr << "\n*** in adjust_tree_structure for " << format_label() << endl;

    hdyn *br = get_root();
    if (is_top_level_node()) {		// this is a top-level node

	bool cm = false;
	if (oldest_daughter) cm = true;

	if (nn == NULL) {
            // pretty_print_node(cerr); cerr << " nn is NULL " << endl;
	    return 0;
        }

	hdyn *nn_top = nn->get_top_level_node();

	if (too_close(this, nn_top, d_min_sq)) {

	    if (number_of_daughter_nodes_exceeds(parent, 2)) {

		predict_loworder_all(get_root(), system_time);

		// cerr << "\ntime = " << system_time << endl;
		// print_recalculated_energies(get_root());
		// pp3(get_root(), cerr);

		// cerr << "Call combine_top_level_nodes from "
		//      << "adjust_tree_structure"<<endl;

		combine_top_level_nodes(nn_top, this, full_dump);

		 // cerr << "Time = " << system_time<< endl;
		 // print_recalculated_energies(get_root());
		 // pp3(get_root(), cerr);

		status = 2;

	    } else {

		status = 0;
	    }

	} else if (too_big(this, d_min_sq * lag_factor)) {

	    // print_recalculated_energies(get_root());

	    predict_loworder_all(get_root(), system_time);

	    // Don't check for slow binary, as too_big should return false
	    // if this is slow.

	    // cerr << "call split_top_level_node 2" << endl;

	    split_top_level_node(this, full_dump);

	    // Better exit next -- 'this' has already been deleted...

	    status = 3;

	} else {

	    status = 0;
	}

    } else {			// this is a low-level node
				// *** check both this and our sister ***
	if (kep == NULL) {

	    hdyn *s = get_binary_sister();
	    hdyn *old_top_level_node = get_top_level_node();

	    status = adjust_low_level_node(this, full_dump);
	    if (!status) status = adjust_low_level_node(s, full_dump);

	    if (status) {

		if (status > 1) status += 2;
		if (old_top_level_node != get_top_level_node())
		    get_top_level_node()->zero_perturber_list();

	    } else {

		// Check sizes of low-level binaries.  In some circumstances
		// it is possible for the center of mass timestep of a
		// binary to become very small as its components recede,
		// with the result that the binary is never split up.
		// Explicitly attempt to avoid that here (Steve, 6/00).

		if (check_binary_params(this) || check_binary_params(s))
		    status = 6;
	    }

#if 0
	    if (system_time > 2295.905
		&& node_contains(get_top_level_node(), "9530"))
		print_tree_structure(this);
#endif

	} else {

	    status = 0;
	}
    }

    return status;
}
