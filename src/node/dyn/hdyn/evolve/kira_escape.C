
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions associated with escapers.
//
// Externally visible functions:
//
//	void check_and_remove_escapers

#include "hdyn.h"
#include <star/dstar_to_kira.h>

local int tree_level(hdyn *b)
{
    int level = 0;		// root = 0
    while (b->get_parent()) {
	level++;
	b = b->get_parent();
    }
    return level;
}

local void indent(int l)
{
    for (int i = 0; i < l; i++) cerr << " ";
}

local bool remove_escapers(hdyn* b,		// root node
			   real rmax,		// maximum allowed distance
			   vec center_pos,	// relative to root
			   vec center_vel)
{
    bool correct_dynamics = false;
    real rmax2 = rmax*rmax;
    int n_esc = 0;

    for_all_daughters(hdyn, b, bi)
	if (square(bi->get_pos()- center_pos) > rmax2)
	    n_esc++;

    if (n_esc <= 0) {
	cerr << endl << "  No escapers" << endl;
	return false;
    }

    cerr << endl << n_esc << " escaper(s):\n";

    hdyn** esc_list = new hdynptr[n_esc];

    real epot0, ekin0, etot0;
    calculate_energies_with_external(b, epot0, ekin0, etot0);

    n_esc = 0;
    for_all_daughters(hdyn, b, bj) {

	// Escape criterion (note that we do NOT check E > 0):
	
	if (square(bj->get_pos()- center_pos) > rmax2) {

	    // Print out hierarchical information on the escaper.

	    for_all_nodes(hdyn, bj, bb) {

		int level = tree_level(bb);

		indent(2*level+2);
		cerr << bb->format_label() << " " << bb->get_mass() << endl;

		if (level <= 1) {
		    indent(2*level+2);
		    cerr << "pos: " << bb->get_pos() - center_pos
			 << "   |pos| = " << abs(bb->get_pos() - center_pos)
			 << endl;
		    indent(2*level+2);
		    cerr << "vel: " << bb->get_vel() - center_vel
			 << "   |vel| = " << abs(bb->get_vel() - center_vel)
			 << endl;
		} else {
		    indent(2*level+2); cerr << "pos: " << bb->get_pos() << endl;
		    indent(2*level+2); cerr << "vel: " << bb->get_vel() << endl;
		}
	    }

	    if (b->get_use_sstar())
		if (has_sstar(bj)) {
		    cerr << "    ";
		    put_state(make_star_state(bj));
		    cerr << endl;
		}
	
//	    if (bj->get_starbase() != NULL) {
//		bj->get_starbase()->print_star_story(cerr);
//	    }

	    bool newl = true;
	    for_all_daughters(hdyn, b, bb) {
		if (bb->get_nn() == bj) {
		    if (newl) {
			cerr << endl;
			newl = false;
		    }
		    cerr << "    reset NB for " << bb->format_label()<<endl;
		    bb->set_nn(NULL);
		}
	    }

	    // Correct the list of perturbed top_level binaries.

	    if (bj->is_perturbed_cm())
		bj->remove_from_perturbed_list();

	    esc_list[n_esc++] = bj;
	    detach_node_from_general_tree(bj);

	    b->inc_mass(-bj->get_mass());	// Not done by detach_node...
	}
    }

    // Note from Steve (7/01):  Used to reset the CM because of the
    // chance we might compute escapers relative to the origin of
    // coordinates (rather than, say, the density center).  Without
    // resetting, recoil would eventually carry the entire system
    // beyond the tidal radius!  Now we should *never* use the
    // coordinate origin in determining escapers, so there is no need
    // to correct here.

    // b->to_com();

    // However, the root node is now always at the center of mass of the
    // particles in the system, so correct the particles and root here.

    b->reset_com();

    // Correct all perturber lists.

    correct_perturber_lists(b, esc_list, n_esc);

    // Should possibly recompute accs and jerks on top-level nodes here too.
    // Not necessary if we are going to reinitialize next.

#if 0
    calculate_acc_and_jerk_on_all_top_level_nodes(b);
#endif

    if (b->n_leaves() > 1) {

	real epot = 0, ekin = 0, etot = 0;
	calculate_energies_with_external(b, epot, ekin, etot);

	cerr << endl;
	PRI(2); PRC(epot0); PRC(ekin0); PRL(etot0);
	PRI(2); PRC(epot); PRC(ekin); PRL(etot);

	// pp3(cm, cerr);

	real de_total = etot - etot0;

	PRI(2); PRL(de_total);

	// Update statistics on energy change components:

	b->get_kira_counters()->de_total += de_total;

	predict_loworder_all(b, b->get_system_time());	// Unnecessary?

	cerr << "initialize_system_phase2 called from remove_escapers\n";
	initialize_system_phase2(b, 4);		// default set_dt
    }

    return true;
}



void check_and_remove_escapers(hdyn* b,
			       xreal& ttmp, hdyn** next_nodes,
			       int& n_next, bool& tree_changed)
{
    if (b->get_scaled_stripping_radius() <= 0) return;

    // Note: Definition of scaled_stripping_radius
    // 	 incorporates initial mass factor.

    real stripping_radius = b->get_scaled_stripping_radius()
				* pow(total_mass(b), 1.0/3.0);

    // Note also that this scaling ensures that omega
    // remains constant as the system evolves if the Jacobi
    // radius is identified with the stripping radius.

    real mass0 = total_mass(b);

    vec center_pos, center_vel;
    get_std_center(b, center_pos, center_vel);

    center_pos -= b->get_pos();			// std_center quantities
    center_vel -= b->get_vel();			// include the root node

    // Note from Steve (7/01):  If we have a tidal field, we should
    // probably determine the center self-consistently as the center
    // of mass of particles within the Jacobi surface.
    //
    // To do...  (See also refine_cluster_mass().)

    if (remove_escapers(b, stripping_radius, center_pos, center_vel)) {

	PRI(2); PRL(center_pos);
	PRI(2); PRC(b->get_scaled_stripping_radius());
	PRL(stripping_radius);

	tree_changed = true;
	fast_get_nodes_to_move(b, next_nodes, n_next,
			       ttmp, tree_changed);
	tree_changed = true;

	cerr << "\n  New N = " << b->n_leaves()
	     <<  "  mass = " << total_mass(b) << endl;
	PRI(2); print_recalculated_energies(b);

	// Update statistics.

	b->get_kira_counters()->dm_escapers += mass0 - total_mass(b);
    }
}
