
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// hdyn_tt:  Starlab hdyn-specific node-handling functions.

// Externally visible functions:
//
//	hdyn* hdyn::next_node
//	hdyn* hdyn::find_perturber_node
//	bool hdyn::nn_stats
//	real hdyn::print_pert
//	void hdyn::setup_binary_node
//
//	void print_nn
//	void print_coll

#include "hdyn.h"

#ifndef TOOLBOX

// Note from Steve 8/20/98.  This is based on the node version, but is
// modified to include the possiblilty of treating unperturbed centers
// of mass as leaves.

hdyn* hdyn::next_node(hdyn* base,
		      bool resolve_unperturbed) // default = true
{
    // If resolve_unperturbed is false, we only descend to the next level
    // in the tree if the oldest daughter is perturbed (Steve, 8/20/98).

    if (oldest_daughter
	&& (resolve_unperturbed
	    || get_oldest_daughter()->kep == NULL
	    || !get_oldest_daughter()->fully_unperturbed)) {

	return get_oldest_daughter();

    } else if (this == base) {		// in case base is a leaf...

	return NULL;

    } else if (younger_sister) {

	return get_younger_sister();

    } else {

	// Can't go down or right.  See if we can go up.

	if (this == base)		// in case base is a leaf...
	    return NULL;

	if (parent == NULL)		// in case b is atomic...
	    return NULL;

	hdyn* tmp = get_parent();
	while (tmp->get_younger_sister() == NULL) {
	    if (tmp == base)
		return NULL;
	    else
		tmp = tmp->get_parent();
	}

	// Now tmp is the lowest-level ancestor with a younger sister.

	if (tmp == base)
	    return NULL;
	else
	    return tmp->get_younger_sister();
    }

    return NULL;	// To keep some compilers happy... 
}

// find_perturber_node:  return a pointer to the first ancestor node
//			 found that contains a valid perturber list,
//			 or NULL otherwise.

hdyn* hdyn::find_perturber_node()
{
    if (is_root() || is_top_level_node()) return NULL;

    hdyn* pnode;
    hdyn* top_level = get_top_level_node();

    if (!ALLOW_LOW_LEVEL_PERTURBERS) {

	// Same result, so avoid the extra work in this case.

	pnode = top_level;

    } else {

	pnode = get_parent();
	while (pnode != top_level && !pnode->valid_perturbers)
	    pnode = pnode->get_parent();

    }

    if (!pnode->valid_perturbers) pnode = NULL;
    return pnode;
}

// print_nn: print out the ID of b's nearest neighbor.

void print_nn(hdyn* b, int level, ostream& s)
{
    if (level > 0)
	s << "nn of " << b->format_label() << " is ";

    if (b->get_nn()) {
	s << "(" << b->get_nn() << ") ";
	if (b->get_nn()->is_valid()) {
	    s << b->get_nn()->format_label();
	    if (level > 1)
		s << " (d = " << sqrt(b->get_d_nn_sq()) << ")";
	} else
	    s << "invalid";
    } else
	s << "(NULL)";

    if (level > 0)
        s << endl;
}

// print_coll: print out the ID of b's coll particle.

void print_coll(hdyn* b, int level, ostream& s)
{
    if (level > 0)
	s << "coll of " << b->format_label() << " is ";

    if (b->get_coll()) {
	s << "(" << b->get_coll() << ") ";
	if (b->get_coll()->is_valid()) {
	    s << b->get_coll()->format_label();
	    if (level > 1)
		s << " (d = " << sqrt(b->get_d_coll_sq())
		  << ", R = " << b->get_coll()->get_radius()
		  << ")";
	} else
	    s << "invalid";
    } else
	s << "(NULL)";

    if (level > 0)
        s << endl;
}

bool hdyn::nn_stats(real energy_cutoff, real kT,
		    vector center, bool verbose,
		    bool long_binary_output,		// default = true
		    int which)				// default = 0
{
    if (which == 0 && verbose)
	cerr << "\n  Bound nn pairs:\n";

    if (!nn) return false;

    // Avoid printing mutual nearest neighbors twice, and print
    // the data from the primary (where possible):

    if (nn->nn == this) {
	if (mass < nn->mass)			    // fails for equal masses
	    return false;
	else if (mass == nn->mass) {

	    // The following *can* occur in the GRAPE-4 version:

#if 0
	    if (nn == this)
		cerr << "    warning: " << format_label()
		     << " is its own nearest neighbor..." << endl;
#endif

	    if (this >= nn)			    // always legal?
	        return false;			    // compare labels if not...
	}
    }

    real M = mass + nn->mass;
    if (M <= 0) return false;

    vector dx, dv;

    if (parent == nn->parent) {

        dx = pos - nn->pos;
	dv = vel - nn->vel;

    } else {

	// Won't happen if nn search is limited to top-level nodes
	// (and probably unnecessary otherwise, as binary should
	// already have been picked up)...

        dx = hdyn_something_relative_to_root(this, &hdyn::get_pos)
	   - hdyn_something_relative_to_root(nn, &hdyn::get_pos);
        dv = hdyn_something_relative_to_root(this, &hdyn::get_vel)
	   - hdyn_something_relative_to_root(nn, &hdyn::get_vel);

    }

    real mu = mass * nn->mass / M;
    real E = mu*(0.5*dv*dv - M / abs(dx));

    // Convention: energy_cutoff is in terms of kT if set,
    //		   in terms of E/mu if kT = 0.

    bool found = false;

    if ((kT > 0 && E < -energy_cutoff*kT)
	|| (kT == 0 && E <  -energy_cutoff*mu)) {

	// cerr << " nn of " << format_label();
	// cerr << " is " << nn->format_label() << endl << flush;

	print_binary_from_dyn_pair(this, nn, kT, center,
				   verbose, long_binary_output);
	cerr << endl;

	found = true;
    }

    return found;
}

real hdyn::print_pert(bool long_binary_output,		// default = true
		      int indent)			// default = BIN_INDENT
{
    // Print out information on perturbers and nn of this low-level node.
    // Formatting is to match sys_stats output.
    // Return zero or the energy of this binary, if unperturned.

    // NOTE: 'this' is a binary component.

    if (long_binary_output) {

	hdyn* parent = get_parent();
	hdyn* nn = parent->get_nn();

	// cerr << endl;		// omitted by print_binary_params...

	if (timestep > 0 && parent->timestep > 0) {
	    PRI(indent);
	    cerr << "cpt timestep = " << timestep
		 << "  CM timestep = " << parent->timestep << endl;
	}

	PRI(indent);
	real pert_sq = 0;
	if(get_perturbation_squared()>=0)
	  pert_sq = sqrt(get_perturbation_squared());
	cerr << "pert = " << pert_sq;

	if (slow) cerr << " [" << get_kappa() << "]";

	cerr << " (";
	hdyn* pnode = find_perturber_node();
	if (pnode && pnode->get_valid_perturbers())
	    cerr << pnode->get_n_perturbers()
		 << " perturbers)";
	else
	    cerr << "no perturber list)";

	cerr << "  nn is ";
	if (nn == NULL || nn == parent)
	    cerr << "unknown" << endl;
	else {
	    cerr << nn->format_label() << endl;
	    real rnn = abs(parent->get_pos() - nn->get_pos());
	    if (rnn > 0) {
		real energy =
		    -(parent->get_mass() + nn->get_mass()) / rnn
			+ 0.5 * square(parent->get_vel()
				       - nn->get_vel());

		// Some of this information may already have been
		// printed out for bound nn pairs, but repeat it
		// here for uniformity.

		// if (energy >= 0) {
		PRI(indent);
	        cerr << "nn mass = " << nn->get_mass()
		     << "  dist = " << rnn
		     << "  E/mu = " << energy << endl;
		// }
	    }
	}

    } else {

	if (slow) cerr << " [" << get_kappa() << "]";
	cerr << endl;

    }

    real e_unp = 0;

    if (kep)
	e_unp -= mass * get_binary_sister()->mass
	    	      / (2 * kep->get_semi_major_axis());

    return e_unp;
}

#define TOLERANCE (1e-12)

inline int within_tolerance(real x, real scale)
{
    if (scale != 0.0) {
	return abs(x) / abs(scale) < TOLERANCE;
    } else {
	return abs(x) < TOLERANCE;
    }
}

void check_consistency_of_node(hdyn * node,
			       hdyn_VMF_ptr get_something,
			       char *id)
{
    hdyn *first_child = node->get_oldest_daughter();
    hdyn *second_child = first_child->get_younger_sister();

    if (second_child->get_younger_sister() != NULL) {
	node->pretty_print_node(cerr);
	cerr << " is not a binary node \n";
	return;
    }
    vector d = (first_child->*get_something) ()
	        - (second_child->*get_something) ();
    vector cm = (first_child->*get_something) () * first_child->get_mass()
	        + (second_child->*get_something) () * second_child->get_mass();
    real mass_of_node = first_child->get_mass() + second_child->get_mass();

    if (!within_tolerance(mass_of_node - node->get_mass(), node->get_mass())) {

	node->pretty_print_node(cerr);
	cerr << " has incorrect total mass "
	     << mass_of_node << " " << node->get_mass();

    } else {

	cm /= mass_of_node;
	if (!within_tolerance(abs(cm), abs(d))) {
	    node->pretty_print_node(cerr);
	    cerr << " has incorrect center of mass " << id
		 << ": " << cm << endl;
	}
    }
}

void check_consistency_of_nodes(hdyn * node)
{
    if (node->get_parent() != NULL) {

	// node is not root

	if (node->get_oldest_daughter() != NULL) {
	  check_consistency_of_node(node, &hdyn::get_pos, "pos");
	  check_consistency_of_node(node, &hdyn::get_vel, "vel");
	  check_consistency_of_node(node, &hdyn::get_pred_pos, "ppos");
	  check_consistency_of_node(node, &hdyn::get_pred_vel, "pvel");
	  check_consistency_of_node(node, &hdyn::get_acc, "acc");
	  check_consistency_of_node(node, &hdyn::get_jerk, "jerk");
	  check_consistency_of_node(node, &hdyn::get_old_acc, "old_acc");
	  check_consistency_of_node(node, &hdyn::get_old_jerk, "old_jerk");
	}
    }

    if (node->get_oldest_daughter() != NULL) {
	for (hdyn * bb = node->get_oldest_daughter(); bb != NULL;
	     bb = bb->get_younger_sister()) {
	    check_consistency_of_nodes(bb);
	}
    }
}

void hdyn::setup_binary_node()
{
    hdyn *younger_daughter = (hdyn *) oldest_daughter->get_younger_sister();
    hdyn *older_daughter = (hdyn *) oldest_daughter;

    // Check if the times of daughters are exactly the same.

    if (older_daughter->get_time() != younger_daughter->get_time()) {
	cerr << "setup_binary_node, times of daughters are different:"
	     << older_daughter->get_time() << " "
	     << younger_daughter->get_time() << "\n";
	exit(1);
    }

    // Set up physical quantities of the parent node.
    // Note that pot, old_acc, old_jerk are not set -- must be recalculated.

    time = older_daughter->time;
    timestep = min(older_daughter->timestep, younger_daughter->timestep);

    mass = older_daughter->mass + younger_daughter->mass;

    real f1 = older_daughter->mass / mass;
    real f2 = younger_daughter->mass / mass;
    pos = f1 * older_daughter->pos + f2 * younger_daughter->pos;
    vel = f1 * older_daughter->vel + f2 * younger_daughter->vel;
    acc = f1 * older_daughter->acc + f2 * younger_daughter->acc;
    jerk = f1 * older_daughter->jerk + f2 * younger_daughter->jerk;

    store_old_force();
    pred_pos = f1 * older_daughter->pred_pos + f2 * younger_daughter->pred_pos;
    pred_vel = f1 * older_daughter->pred_vel + f2 * younger_daughter->pred_vel;

    // Set up the physical coordinates of the daughters.

    older_daughter->pos -= pos;
    older_daughter->vel -= vel;
    older_daughter->acc -= acc;
    older_daughter->jerk -= jerk;
    older_daughter->store_old_force();
    older_daughter->pred_pos -= pred_pos;
    older_daughter->pred_vel -= pred_vel;

    younger_daughter->pos -= pos;
    younger_daughter->vel -= vel;
    younger_daughter->acc -= acc;
    younger_daughter->jerk -= jerk;
    younger_daughter->store_old_force();
    younger_daughter->pred_pos -= pred_pos;
    younger_daughter->pred_vel -= pred_vel;
}

void create_binary_from_toplevel_nodes(hdyn * bi, hdyn * bj)
{
    if (!bi->is_top_level_node())
	err_exit("create_binary_from_toplevel_nodes: bi not top level node");
    if (!bj->is_top_level_node())
	err_exit("create_binary_from_toplevel_nodes: bj not top level node");

    detach_node_from_general_tree(*bj);
    bj->set_younger_sister(NULL);
    bj->set_elder_sister(NULL);

    hdyn *new_n = new hdyn();
    insert_node_into_binary_tree(*bj, *bi, *new_n);
    bj->get_parent()->setup_binary_node();

    label_binary_node(bj->get_parent());
}

local vector change_of_absolute_something_of_parent(hdyn * node,
						    vector & dx,
						    real dm,
						    hdyn_VMF_ptr something)
{
    return ((node->get_mass() + dm) * dx + (node->*something) () * dm)
	 / (node->get_mass() + node->get_binary_sister()->get_mass() + dm);
}

local void adjust_relative_something_of_sister(vector
					         absolute_change_of_parent,
					       real dm,
					       hdyn * node,
					       hdyn_MF_ptr inc_something)
{
    hdyn *sister;
    sister = (hdyn *) (node->get_binary_sister());

    (sister->*inc_something) (-absolute_change_of_parent);
    real dummy = dm;		// to keep the compiler happy

}

local void adjust_parent_and_sister(hdyn * node,
				    hdyn * ancestor,
				    vector d_something,
				    real dm,
				    hdyn_VMF_ptr get_something,
				    hdyn_MF_ptr inc_something,
				    int modify_node)
{
//    cerr << "in adjust_parent_and_sister...\n" << flush;

    vector d_parent = 0;

    if (node != ancestor) {
	if(node->is_low_level_node()){
	    d_parent = change_of_absolute_something_of_parent(node,
							      d_something,
							      dm,
							      get_something);
	    adjust_relative_something_of_sister(d_parent, dm, node,
						inc_something);
	    if (node->get_parent() != ancestor)
		adjust_parent_and_sister(node->get_parent(), ancestor,
					 d_parent, dm,
					 get_something, inc_something, 1);
	}
    } else {
	d_parent = 0;
    }

    if (modify_node)
	(node->*inc_something) (d_something - d_parent);
}

local void init_pred_of_parent_and_sister(hdyn * node,
					  hdyn * ancestor)
{
    if (node != ancestor) {
	if(node->is_low_level_node()){
	    node->get_parent()->init_pred();
	    node->get_binary_sister()->init_pred();
	    if (node->get_parent() != ancestor)
	    init_pred_of_parent_and_sister(node->get_parent(), ancestor);
	}
    }
    node->init_pred();
}

local void adjust_parent_mass(hdyn * node,
			      hdyn * ancestor,
			      real dm)
{
    node->set_mass(node->get_mass() + dm);

    if (node->get_parent() != ancestor)
	adjust_parent_mass(node->get_parent(), ancestor, dm);
}

void remove_node_and_correct_upto_ancestor(hdyn * ancestor, hdyn * node)
{
    // cerr << "in remove_node_and_correct_upto_ancestor" << endl;

    if (node->is_top_level_node()) {

	detach_node_from_general_tree(*node);

    } else {

	real dm = -node->get_mass();
	vector dz;
	dz = 0;

	hdyn *parent = node->get_parent();
	hdyn *sister = (hdyn *) (node->get_binary_sister());

	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_pos,
				 &hdyn::inc_pos, 0);
	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_vel,
				 &hdyn::inc_vel, 0);
	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_acc,
				 &hdyn::inc_acc, 0);
	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_old_acc,
				 &hdyn::inc_old_acc, 0);
	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_jerk,
				 &hdyn::inc_jerk, 0);
	adjust_parent_and_sister(node, ancestor, dz, dm,
				 &hdyn::get_old_jerk,
				 &hdyn::inc_old_jerk, 0);

	adjust_parent_mass(node, ancestor, dm);
	init_pred_of_parent_and_sister(node, ancestor);

	sister->set_pos(parent->get_pos());
	sister->set_vel(parent->get_vel());
	sister->set_acc(parent->get_acc());
	sister->set_jerk(parent->get_jerk());
	sister->store_old_force();

	// Sister pred_xxx may be needed before the next step.
	// Not currently correct -- fix that here.

	sister->init_pred();

	detach_node_from_binary_tree(*node);

	// parent must be deleted.... (17-Aug-1996)

	delete parent;
    }
}

vector something_relative_to_ancestor(hdyn * bj,
				      hdyn * bi,
				      hdyn_VMF_ptr get_something)
{
#ifdef DEBUG
    cerr << "something_relative_to_ancestor\n";
#endif
    vector d_something = 0;
    for (hdyn * b = bi; b != bj; b = b->get_parent()) {
	if (b == NULL) {
	    cerr << "something_relative_to_ancestor: Error, bj is not the "
		 << "ancestor of bi\n";
	    exit(1);
	}
	d_something += (b->*get_something)();
    }
    return d_something;
}

vector hdyn_something_relative_to_root(hdyn * bi,
				       hdyn_VMF_ptr get_something)
{
#ifdef DEBUG
    cerr << "hdyn_something_relative_to_root\n";
#endif

    // Special treatment of top-level nodes:

    if (bi->is_top_level_node())

	return (bi->*get_something) ();

    else {

	vector d_something = 0;
	for (hdyn * b = bi; b->get_parent() != NULL; b = b->get_parent()) {
	    if (b == NULL)
		err_exit(
		"hdyn_something_relative_to_root: bj not the ancestor of bi");
	    d_something += (b->*get_something)();
	}
	return d_something;
    }
}

local vector relative_something(hdyn * ancestor,
				hdyn * bj,
				hdyn * bi,
				hdyn_VMF_ptr get_something)
{
#ifdef DEBUG
    cerr << "relative_something\n";
#endif
    return (something_relative_to_ancestor(ancestor, bj, get_something)
	    - something_relative_to_ancestor(ancestor, bi, get_something));
}

local void calculate_new_physical_quantities(hdyn * ancestor,
					     hdyn * node,
					     hdyn * new_sister)
{
    node->set_pos(relative_something(ancestor, node, new_sister,
				     &hdyn::get_pos));
    node->set_vel(relative_something(ancestor, node, new_sister,
				     &hdyn::get_vel));
    node->set_acc(relative_something(ancestor, node, new_sister,
				     &hdyn::get_acc));
    node->set_old_acc(relative_something(ancestor, node, new_sister,
					 &hdyn::get_old_acc));
    node->set_jerk(relative_something(ancestor, node, new_sister,
				      &hdyn::get_jerk));
    node->set_old_jerk(relative_something(ancestor, node, new_sister,
					  &hdyn::get_old_jerk));
}

local void insert_node_and_correct_upto_ancestor(hdyn * ancestor,
						 hdyn * node,
						 hdyn * new_sister)
{
    if (new_sister->get_parent() == NULL) {
	add_node(*node, *new_sister);	// Convention for new top-level node

    } else {
	hdyn *new_n = new hdyn();
	insert_node_into_binary_tree(*node, *new_sister, *new_n);

	hdyn *parent = node->get_parent();
	parent->set_mass(new_sister->get_mass());
	parent->set_pos(new_sister->get_pos());
	parent->set_vel(new_sister->get_vel());
	parent->set_acc(new_sister->get_acc());
	parent->set_jerk(new_sister->get_jerk());
	parent->store_old_force();
	new_sister->set_pos(0);
	new_sister->set_vel(0);
	new_sister->set_acc(0);
	new_sister->set_jerk(0);
	new_sister->set_old_acc(0);
	new_sister->set_old_jerk(0);

	// Attach the node with zero mass (so no corrections needed),
	// then correct the mass to the proper value.

	real dm = node->get_mass();
	node->set_mass(0);
	static vector dx, dv, da, dj;
	dx = dv = da = dj = 0;

	adjust_parent_and_sister(node, ancestor, dx, dm,
				 &hdyn::get_pos,
				 &hdyn::inc_pos, 1);
	adjust_parent_and_sister(node, ancestor, dv, dm,
				 &hdyn::get_vel,
				 &hdyn::inc_vel, 1);
	adjust_parent_and_sister(node, ancestor, da, dm,
				 &hdyn::get_acc,
				 &hdyn::inc_acc, 1);
	adjust_parent_and_sister(node, ancestor, da, dm,
				 &hdyn::get_old_acc,
				 &hdyn::inc_old_acc, 1);
	adjust_parent_and_sister(node, ancestor, dj, dm,
				 &hdyn::get_jerk,
				 &hdyn::inc_jerk, 1);
	adjust_parent_and_sister(node, ancestor, dj, dm,
				 &hdyn::get_old_jerk,
				 &hdyn::inc_old_jerk, 1);
	adjust_parent_mass(node, ancestor, dm);
	init_pred_of_parent_and_sister(node, ancestor);

	label_binary_node(ancestor);
    }
}

void correct_leaf_for_change_of_mass(hdyn * node, real dm)  
{
    hdyn *parent = node->get_parent();
    hdyn *ancestor = node->get_root();
    static vector dx, dv, da, dj;
    dx = dv = da = dj = 0;
    if(node->is_low_level_node()){
	adjust_parent_and_sister(node, ancestor, dx, dm,
				 &hdyn::get_pos,
				 &hdyn::inc_pos, 1);
	adjust_parent_and_sister(node, ancestor, dv, dm,
				 &hdyn::get_vel,
				 &hdyn::inc_vel, 1);
	adjust_parent_and_sister(node, ancestor, da, dm,
				 &hdyn::get_acc,
				 &hdyn::inc_acc, 1);
	adjust_parent_and_sister(node, ancestor, da, dm,
				 &hdyn::get_old_acc,
				 &hdyn::inc_old_acc, 1);
	adjust_parent_and_sister(node, ancestor, dj, dm,
				 &hdyn::get_jerk,
				 &hdyn::inc_jerk, 1);
	adjust_parent_and_sister(node, ancestor, dj, dm,
				 &hdyn::get_old_jerk,
				 &hdyn::inc_old_jerk, 1);
	adjust_parent_mass(node, ancestor, dm);
	init_pred_of_parent_and_sister(node, ancestor);
    }else{
	node->set_mass(node->get_mass() + dm);
    }
}

void correct_leaf_for_change_of_vector(hdyn * node,
				       vector d_something,
				       hdyn_VMF_ptr get_something,
				       hdyn_MF_ptr inc_something)
{
    hdyn *ancestor = node->get_root();
    adjust_parent_and_sister(node, ancestor, d_something, 0,
			     get_something, inc_something, 1);
    init_pred_of_parent_and_sister(node, ancestor);
}


// move_node: reattach a node as a sister of another node.

void move_node(hdyn * node_to_move,
	       hdyn * place_to_insert)
{
    // cerr << "in move_node" << endl << flush;

    // Make a copy of node_to_move prior to any tree operations.

    static hdyn tmp;		// NECESSARY to avoid an occasional unexplained
			        // Bus Error on return with g++ version 2.2.2!!

    tmp = *node_to_move;	// NECESSARY to ensure proper initialization!

    hdyn *ancestor = (hdyn *) common_ancestor(node_to_move, place_to_insert);

    calculate_new_physical_quantities(ancestor, &tmp, place_to_insert);

    // The parent of node_to_move will be deleted by function
    // remove_node_and_correct_upto_ancestor.  Correct any nn
    // references, at least within the local clump, before the
    // tree structure is modified.

    for_all_nodes(hdyn, node_to_move->get_top_level_node(), bb) {
	if (bb->get_nn() == node_to_move->get_parent()) {
	    bb->set_nn(NULL);
	    // cerr << "set nn of " << bb->format_label() << " NULL" << endl;
	}
    }

    hdyn* temp_ancestor = NULL;
    if (node_to_move->get_parent() == ancestor) {

	// In this case, the ancestor will be removed, so we must reset it.

	temp_ancestor = (hdyn *) node_to_move->get_binary_sister();
    }
    remove_node_and_correct_upto_ancestor(ancestor, node_to_move);

    real temp_mass;
    if (temp_ancestor) {
	temp_mass = ancestor->get_mass();
	ancestor = temp_ancestor;
    }

    // Restore the original node_to_move.

    *node_to_move = tmp;
    insert_node_and_correct_upto_ancestor(ancestor, node_to_move,
					  place_to_insert);

    if (temp_ancestor)
	ancestor->set_mass(temp_mass);

    // Bus error sometimes occurs AFTER last executable line here...

    // Update times and labels:

    node_to_move->get_parent()->set_time(place_to_insert->get_time());
    node_to_move->get_parent()->set_timestep(place_to_insert->get_timestep());

    // cerr << "leaving move_node" << endl << flush;
}


// move_node_by_index: reattach node #i as a sister of node #j

void move_node_by_index(int i, int j, hdyn * root)
{
    move_node((hdyn *) node_with_index(i, root),
	      (hdyn *) node_with_index(j, root));
}

//-------------------------------------------------------------------------

// These are hdyn versions of functions in dyn/util/dyn_di.C -- they
// perform the same operations and have the same calling sequences,
// but they also set the hdyn::pot member data.

local void accumulate_energies(hdyn * root, hdyn * b, real eps2,
			       real & epot, real & ekin, real & etot,
			       bool cm = false)
{
    if (b->get_oldest_daughter()) {

	// Start/continue the recursion.  In the cm = true case,
	// expand daughters only if root really is root.

	if (!cm || b->is_root())
	    for_all_daughters(hdyn, b, bb)
		accumulate_energies(root, bb, eps2, epot, ekin, etot, cm);

	etot = epot + ekin;
    }

    // Do the work.

    if (!b->is_root()) {
	real m = b->get_mass();

	// pot_on_general_node is a dyn function returning the total
	// potential per unit mass of a top-level node, or the potential
	// per unit mass of a low-level node due to its binary sister

	real pot = m * pot_on_general_node(root, b, eps2, cm);

	epot += 0.5 * m * pot_on_general_node(root, b, eps2, cm);
	ekin += 0.5 * m * square(b->get_vel());

	b->set_pot(pot);	// will only be used for top-level nodes
    }
}

void calculate_energies(hdyn * root, real eps2,
			real & epot, real & ekin, real & etot,
			bool cm)	// default = false
{
    epot = ekin = etot = 0;
    accumulate_energies(root, root, eps2, epot, ekin, etot, cm);
}

#endif
