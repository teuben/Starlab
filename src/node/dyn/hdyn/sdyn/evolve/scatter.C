
//// scatter:  simple N-body scattering program.
////
//// Options:      -A    specify accuracy parameter [0.02]
////               -c    specify CPU time check interval (s) [3600]
////               -C    specify cube size for snapshot output [10]
////               -d    specify log output interval [none]
////               -D    specify snapshot output interval [none]
////               -i    specify mkscat-style string to initialize
////                        the integration ["-M 1 -v 2 -t -p"]
////               -p    read initial state from cin [false]
////               -t    specify time span of integration [1]
////               -v    verbose mode [false]

#ifdef TOOLBOX

#include "scatter.h"

// tree_is_unbound: return TRUE iff all the top-level nodes of the specified
//                  tree are unbound, widely separated, and receding.
//
// Specifically, require that all top-level pairs be unbound and receding,
// and that no node be a significant perturber of any other pair.
//
// Since deeper levels of structure are already taken into account in the
// tree structure, we only need to check that all nodes at the top level
// are unbound and separating.
//
// Note: If we want to apply analytic extension of nearly unbound top-level
//	 nodes, it should be done here.

local bool tree_is_unbound(sdyn* root, int debug)
{
    if (debug) {
	cerr << "Top level nodes: ";
	for_all_daughters(sdyn, root, bb) cerr << " " << id(bb);
	cerr << endl;
    }

    root->to_com();      // Move to the center-of-mass frame.

    real kin = 0;
    real pot = 0;

    for_all_daughters(sdyn, root, bi) {
	for (sdyn* bj = bi->get_younger_sister();
	     bj != NULL; bj = bj->get_younger_sister()) {

	    if (debug) cerr << "checking i = " << id(bi)
			    << "  j = " << id(bj) << endl;
	    
	    // Test radial velocity, separation, and relative energy of (i,j):

	    if ((bi->get_pos() - bj->get_pos())
		 * (bi->get_vel() - bj->get_vel()) < 0) return FALSE;

	    real rij = abs(bi->get_pos() - bj->get_pos());

	    if (debug) cerr << "    rij = " << rij << endl;

	    real mij = bi->get_mass() + bj->get_mass();
	    real mu_scale = TIDAL_TOL_FACTOR * bi->get_mass() / mij;

	    // (mij here for the same reason as in scatter3/triple_escape.)

	    real rlimit = max(LARGE_SEPARATION,
			      bi->get_radius() * pow(mu_scale, -1/3.0));

	    if (debug) cerr << "    rlimit = " << rlimit << endl;

	    if (rij < rlimit) return FALSE;

	    real scaled_safety_factor = ENERGY_SAFETY_FACTOR * rlimit / rij;

	    real kij = 0.5 * square(bi->get_vel() - bj->get_vel());
	    real pij = mij / rij;
	    real eij = kij - pij;

	    if (debug) cerr << "    kij = " << kij
		 << "  pij = " << pij
		 << "  eij = " << kij - pij
		 << endl;

	    if (eij < 0) return FALSE;
	    if (abs(kij/pij - 1) < scaled_safety_factor) return FALSE;

	    real aij = 0.5 * mij / eij;
	    vector cmij = (bi->get_mass() * bi->get_pos()
			     + bj->get_mass() * bj->get_pos()) / mij;

	    // Check the perturbations of all other particles on (i,j).

	    for_all_daughters(sdyn, root, bk)
		if (bk != bi && bk != bj) {

		    real rk = abs(bk->get_pos() - cmij);

		    if (debug) cerr << "    checking perturber " << id(bk)
			 << " at distance " << rk << "...";

		    if (rk < aij * pow(TIDAL_TOL_FACTOR * mij
				       / (bk->get_mass() + mij), -1/3.0)) {
			if (debug) cerr << "too close" << endl;
			return FALSE;
		    }
		    if (debug) cerr << endl;
		}

	    if (debug) cerr << "    done" << endl;

	    pot -= bi->get_mass() * bj->get_mass() / rij;
	}
	kin += 0.5 * bi->get_mass() * square(bi->get_vel());
    }

    // Finally, check total energy.

    // cerr << "total system energy = " << kin + pot << endl;

    if (kin + pot <= 0) return FALSE;
    return TRUE;
}

#define DYNAMICS	1
#define STABILITY	1
#define K_MAX		2
#define DT_CHECK	20

// scatter: Take the system with root node b and integrate it forward
// 	    in time up to time delta_t.  Check the state of the system
//	    every DT_CHECK time units.

void scatter(sdyn* b, real eta,
	     real delta_t, real dt_out, real cpu_time_check,
	     real dt_snap, real snap_cube_size,
	     int debug)
{
    int stop_at_cpu = 0;
    real cpu_save = cpu_time();

    for (real t = 0; t <= delta_t; t += DT_CHECK) {

	tree_evolve(b, 20, dt_out, dt_snap, snap_cube_size, eta,
		    cpu_time_check);

	if (debug) {
	    real kin, pot;
	    calculate_energy(b, kin, pot);
	    int p = cerr.precision(STD_PRECISION);
	    cerr << "\nStatus at time " << b->get_time();
	    cerr.precision(INT_PRECISION);
	    cerr << " (energy = " << kin + pot << "):\n";
	    cerr.precision(p);
	}

	// Check to see if the scattering is over.

	if (cpu_time_check < 0
	    && cpu_time() - cpu_save > abs(cpu_time_check)) return;

	make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);

	if (tree_is_unbound(b, debug)) return;
   }
}

local void pp(sdyn* b, ostream & s, int level = 0) {

    int p = s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    b->pretty_print_node(s);
    s << " \t"<< b->get_mass() << " \t"
      << b->get_pos() << "   "
      << b->get_vel() <<endl;

    for (sdyn * daughter = b->get_oldest_daughter();
	 daughter != NULL;
	 daughter = daughter->get_younger_sister())
	pp(daughter, s, level + 1);	

    cerr.precision(p);
}

main(int argc, char **argv)
{
    sdyn* b;             // pointer to the nbody system
    
    real  delta_t = VERY_LARGE_NUMBER;
                         // time span of the integration
    real  eta = 0.02;    // time step parameter (for fixed time step,
                         //   equal to the time step size; for variable
                         //   time step, a multiplication factor)

    real  dt_out = VERY_LARGE_NUMBER;
                         // output time interval

    bool  D_flag = FALSE;
    real  dt_snap = VERY_LARGE_NUMBER;
                         // snapshot time interval
    real snap_cube_size = 10;
    real cpu_time_check = 3600;

    char* default_init
		= "-M 1 -v 2 -t -p";    // Head-on collision between identical
    char* init_string = default_init;	// unit-mass binaries with v_inf = 2.

    int pipe = 0;	// Default is to make data with mkscat

    int debug = 0;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:i:pt:v";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'A': eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = atof(poptarg);
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      break;
	    case 'D': D_flag = TRUE;
		      dt_snap = atof(poptarg);
		      break;
	    case 'i': init_string = poptarg;
		      break;
	    case 'p': pipe = 1 - pipe;
		      break;
	    case 't': delta_t = atof(poptarg);
		      break;
	    case 'v': debug = 1 - debug;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
	}            

    if (!D_flag) dt_snap = delta_t; // Guarantee output at end

    cpu_init();

    // Construct the system according to the input string, or
    // read it from stdin:

    if (pipe)
	b = get_sdyn(cin);
    else
	b = mkscat(init_string);

    b->log_history(argc, argv);

    // Force zero radii (note that in the two-body case, the target will
    // have unit radius, according to mkscat...).

    for_all_leaves(sdyn, b, bb) bb->set_radius(0);

    // Flatten and remake the tree, to confirm that mktree did what we expect.

    b->flatten_node();
    real kin, pot;
    real etot_init = calculate_energy_from_scratch(b, kin, pot);
    make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);

    // Initial output:

    cerr << "*** Initial configuration (random seed = "
	 << get_initial_seed() << "):\n";
    b->set_name("root");
    pp(b, cerr);

    cerr << "Total energy = " << etot_init << endl;
    cerr << "Normal form:  ";
    print_normal_form(b, cerr);

    b->flatten_node();	// Needed because low-level routines don't
			// understand tree structure (yet).

    // Integrate the system to completion:

    scatter(b, eta,
	    delta_t, dt_out, cpu_time_check,
	    dt_snap, snap_cube_size,
	    debug);

    cerr.precision(STD_PRECISION);
    if (debug) cerr << endl;

    b->flatten_node();
    real etot_error = calculate_energy_from_scratch(b, kin, pot) - etot_init;
    make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);

    // Final output:

    cerr << "*** Final system configuration (time = "
	 << b->get_time() << "):\n";
    b->set_name("root");
    pp(b, cerr);

    cerr << "Energy error = " << etot_error << endl;

    cerr << "Normal form:  ";
    print_normal_form(b, cerr);
}

#endif
