
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// hdyn_sync:  Print out various diagnostic diagnostics on the time
////             and time steps in the input system.
////
////             The user should modify this program to correct any
////             problems revealed by the diagnostics (too hard to
////             automatize this process).

#include "hdyn.h"
#include "util_io.h"

#ifdef TOOLBOX

#define PRX(x) cerr << #x << " = "; xprint(x);

local bool is_power_of_two(real dt)
{
    real p = 1;
    int count = 0;

    while (p != dt && ++count < 64) p /= 2;

    if (p == dt)
	return true;
    else
	return false;
}

main(int argc, char **argv)
{
    check_help();

    hdyn *b;

    if (b = get_hdyn(cin)) {

	// Determine min and max top-level time.

	xreal system_time = b->get_system_time();
	xreal tmin, tmax;

	PRX(system_time);

	// Top-level nodes:

	int n_top = 0;
	tmin = 1.e10, tmax = 0;

	for_all_daughters(hdyn, b, bb) {
	    n_top++;
	    xreal tt = bb->get_time();
	    if (tt < tmin) tmin = tt;
	    if (tt > tmax) tmax = tt;
	}

	cerr << n_top << " top-level nodes" << endl;
	PRI(4); PRX(tmin);
	PRI(4); PRX(tmax);

	// Low-level perturbed nodes:

	int n_low = 0;
	tmin = 1.e10, tmax = 0;

	for_all_leaves(hdyn, b, bb)
	    if (bb->is_low_level_node()
		&& bb->get_younger_sister()
		&& !bb->get_kepler()) {
		n_low++;
		xreal tt = bb->get_time();
		if (tt < tmin) tmin = tt;
		if (tt > tmax) tmax = tt;
	    }

	cerr << n_low << " low-level perturbed nodes" << endl;
	if (tmax > 0) {
	    PRI(4); PRX(tmin);
	    PRI(4); PRX(tmax);
	} else
	    cerr << "    (none)" << endl;

	// Low-level unperturbed nodes:

	int n_pert = 0;
	tmin = 1.e10, tmax = 0;

	for_all_leaves(hdyn, b, bb)
	    if (bb->is_low_level_node()
		&& bb->get_younger_sister()
		&& bb->get_kepler()) {
		n_pert++;
		xreal tt = bb->get_time();
		if (tt < tmin) tmin = tt;
		if (tt > tmax) tmax = tt;
	    }

	cerr << n_pert << " low-level unperturbed nodes" << endl;
	PRI(4); PRX(tmin);
	PRI(4); PRX(tmax);

	// Timesteps (top-level/perturbed nodes):

	real dtmin = 1.e10, dtmax = 0;
	bool powers_of_two = true;
	int nbad2 = 0;

	int noncom = 0, com = 0;

	// Assume steps are powers of 2 and less than 1 in
	// determining commensurability.

	for_all_nodes(hdyn, b, bb) 
	    if (bb->is_top_level_node()
		|| (bb->get_younger_sister()
		    && !bb->get_kepler())) {
		real dt = bb->get_timestep();
		dtmin = Starlab::min(dtmin, dt);
		dtmax = Starlab::max(dtmax, dt);
		if (!is_power_of_two(dt)) {
		    powers_of_two = false;
		    nbad2++;
		} else {
		    if (fmod2(system_time, dt) == 0)
			com++;
		    else
			noncom++;
		}
	    }

	bool commensurate = (noncom == 0);

	cerr << "perturbed steps:" << endl;
	PRI(4); PRC(dtmin); PRC(dtmax); PRL(nbad2);
	if (nbad2 == 0) {
	    PRI(4); PRC(commensurate); PRC(com); PRL(noncom);
	}

	// Now, what to do if there are problems?
	//
	// Probably need to modify this program as actual problems arise.
	// Expect (in normal circumstances):
	//
	//	all times <= system time		(necessary)
	//	all top-level and perturbed low-level
	//	    node times = system time		(expected)
	//	unperturbed low-level nodes time may be
	//	    less than system time, but time steps
	//	    should span system time		(necessary)
	//	time steps should be commensurate with
	//	    system time				(necessary)
	//
	// What to do if these conditions are not met?
	//
	//     span in top-level/pert low-level node times
	//     unperturbed time step doesn't reach system time
	//     time steps not commensurate
	//
	// Hard to define general actions.  User may simply need
	// to modify this program to make corrections.  Code below
	// deals with simple case of non-commensurate times.
	//
	// Templates:

	// May need to modify the unperturbed times by the same increment
	// as the other times...

	xreal dt = (xreal)41 - system_time;

	// Top-level nodes:

	b->set_system_time(41);
	b->set_time(41);		// bad that these can differ...

	for_all_daughters(hdyn, b, bb) {

	    // ...

	    bb->set_time(41);
	}

	// Low-level perturbed nodes:

	for_all_leaves(hdyn, b, bb)
	    if (bb->is_low_level_node()
		&& bb->get_younger_sister()
		&& !bb->get_kepler()) {

		// ...

		bb->set_time(41);
	    }

	// Low-level unperturbed nodes:

	for_all_leaves(hdyn, b, bb)
	    if (bb->is_low_level_node()
		&& bb->get_younger_sister()
		&& bb->get_kepler()) {

		// ...

		real t = bb->get_time()+dt;
		bb->set_time(t);
		bb->get_younger_sister()->set_time(t);
	    }

	put_node(cout, *b);
    }
}

#endif
