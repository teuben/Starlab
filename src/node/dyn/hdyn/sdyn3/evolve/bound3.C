
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// bound3:   Perform three-body experiments starting from bound
////           configurations.
////
////           The program initializes a three-body system in virial
////           equilibrium within a specified region, and by
////           default integrates the system until a clear dynamical
////           outcome is recognized.
////
////           Default program output is intermediate and final states.
////
////           Units: G = 1, total mass = 1, virial radius = 1.
////
//// Options:   -A    specify accuracy parameter [0.05]
////            -b    print positions and velocities in report [no]
////            -c    specify CPU time check, in hours [1]
////            -C    specify snap cube size [10]
////            -d    specify log output interval [none]
////            -D    specify snap output interval [none]
////            -e    specify initial eccentricity [0]
////            -g    specify tidal tolerance [1.e-6]
////            -L    specify minimum initial outer separation [10]
////            -m    specify secondary mass (binary mass = 1) [0.5]
////            -M    specify incomer mass [0.5]
////            -n    specify number of experiments [1]
////            -N    specify random number count [0]
////            -o    specify outer orbit orientation [random]
////            -p    force planar prograde [not forced]
////            -P    force planar retrograde [not forced]
////            -q    minimal output [false]
////            -Q    intermediate amount of output [false]
////            -r    specify incomer impact parameter [0]
////            -R    specify separation at which to start and stop [none]
////            -s    specify random seed [taken from system clock]
////            -S    specify separation at which to stop [none]
////            -U    specify maximum initial outer separation [none]
////            -v    specify incomer velocity at infinity [0 ==> Etot = 0]
////            -x    specify primary radius [0]
////            -y    specify secondary radius [0]
////            -z    specify incomer radius [0]
////
//// As a convenient shorthand, any "dt" interval specified less than zero
//// is interpreted as a power of 2, i.e. "-d -3" sets dt_out = 0.125.

// Input to the function is an initial state, output is intermediate
// and final states.

// Note the function bound3 is completely deterministic.
// No randomization is performed at this level.

// Starlab library function.

#include "scatter3.h"

#ifdef TOOLBOX

// check_init:  Make sure an initial state is sensible.

local void check_init(initial_state3 & init)
{
    if (init.m2 > 1) err_exit("check_init: m2 > 1");
    if (init.m2 < 0) err_exit("check_init: m2 < 0");
    if (init.m3 > 1) err_exit("check_init: m3 > 0");
    if (init.m3 < 0) err_exit("check_init: m3 < 0");
    if (init.r1 < 0) err_exit("check_init: r1 < 0");
    if (init.r2 < 0) err_exit("check_init: r2 < 0");
    if (init.r3 < 0) err_exit("check_init: r3 < 0");
    if (init.eta <= 0) err_exit("check_init: eta <= 0");
}

// bound3:  Perform a three-body experiment, initializing from a specified
//          state and returning intermediate- and final-state structures.

void bound3(sdyn3 * b,
	    initial_state3 & init,
	    intermediate_state3 & inter,
	    final_state3 & final,
	    real cpu_time_check,
	    real dt_out,         // diagnostic output interval
	    real dt_snap,        // snapshot output interval
	    real snap_cube_size) // limit output to particles within cube
{
    inter.id = final.id = init.id;	    // For bookkeeping...

    mksphere(b, 3);			    // Note that all systems are
    initialize_bodies(inter.system);	    // initialized as "non-particles"
    initialize_bodies(final.system);	    // until explicitly set.

    // Save some initial data:

    sdyn3_to_system(b, init.system);

//    put_node(b, cerr);
//    print_bodies(cerr, init.system, STD_PRECISION);

    real e_init = energy(b);
    real cpu_init = cpu_time();
    real cpu = cpu_init;
    
    real kinetic, potential;
    get_top_level_energies(b, 0.0, potential, kinetic);

//    PRC(e_init); PRC(potential); PRL(kinetic);

    // Evolve the system until the interaction is over or a collision occurs.

    inter.n_kepler = 0;

    do {

	int n_stars = 0;
	for_all_daughters(sdyn3, b, bb) n_stars++;

	real check_int = CHECK_INTERVAL;
//	check_int = 0.5;

	tree3_evolve(b, check_int, dt_out, dt_snap, snap_cube_size,
		     init.eta, cpu_time_check); // , dt_print, p);

	// Check the CPU time.  Note that the printed CPU time is the
	// time since this routine was entered.

	if (cpu_time() - cpu > cpu_time_check) {
	    if (cpu == cpu_init) cerr << endl; // (May be waiting in mid-line!)

	    int p = cerr.precision(STD_PRECISION);
	    cerr << "bound3:  CPU time = " << cpu_time() - cpu_init
		 << "  time = " << b->get_time()
		 << "  n_steps = " << b->get_n_steps()
		 << endl;
	    cerr << "           n_osc = " << b->get_n_ssd_osc()
		 << "  n_kepler = " << inter.n_kepler
		 << endl << flush;
	    cpu = cpu_time();
	    cerr.precision(p);
	}

	inter.n_osc = b->get_n_ssd_osc();	// For convenience

	sdyn3_to_system(b, inter.system);
	merge_collisions(b);

	// See if the number of stars has decreased.

	for_all_daughters(sdyn3, b, bbb) n_stars--;
	if (n_stars > 0) {
	    if (dt_snap < VERY_LARGE_NUMBER
		&& system_in_cube(b, snap_cube_size)) {
		put_node(b);
		cout << flush;
	    }
	}

//	cout << b->get_time() << " " << energy(b)-e_init << endl;

    } while (// b->get_n_steps() < 1000000 &&
	     !extend_or_end_scatter(b, init, &inter, &final));

    // Set intermediate state:

    inter.r_min_min = VERY_LARGE_NUMBER;
    inter.n_stars = 0;

    for_all_daughters(sdyn3, b, bb) {
        inter.index[inter.n_stars] = bb->get_index();
        inter.r_min[inter.n_stars] = sqrt(bb->get_min_nn_dr2());
	inter.r_min_min = Starlab::min(inter.r_min_min,
				       inter.r_min[inter.n_stars]);
        inter.n_stars++;
    }

    if (b->get_n_ssd_osc() <= 1)
	inter.descriptor = non_resonance;
    else {
	int n_no_nn_change = 0;
	for_all_daughters(sdyn3, b, bbb)
	    if (bbb->get_nn_change_flag() == 0)
	        n_no_nn_change++;
	if (n_no_nn_change >= 2)
	    inter.descriptor = hierarchical_resonance;
	else
	    inter.descriptor = democratic_resonance;
    }

    // Set final state:

    final.time = b->get_time();
    final.n_steps = (int) b->get_n_steps();
    final.error = energy(b) - e_init; // *Absolute* error, note!

    sdyn3_to_system(b, final.system);

    //  Check for integration errors (relax tolerance in the case of mergers).

    if (abs(final.error) > ENERGY_TOLERANCE)
	if ((   final.descriptor != merger_binary_1
	     && final.descriptor != merger_binary_2
	     && final.descriptor != merger_binary_3
	     && final.descriptor != merger_escape_1
	     && final.descriptor != merger_escape_2
	     && final.descriptor != merger_escape_3
	     && final.descriptor != triple_merger)
	    || abs(final.error)
	         > MERGER_ENERGY_TOLERANCE * abs(potential_energy(b)))
	    final.descriptor = error;

    // In all cases, the final system should be a properly configured
    // 1-, 2-, or 3-body system.  Print it for the last time, to reflect
    // final state, if the "-D" command-line option was specified.

    if (dt_snap < VERY_LARGE_NUMBER && system_in_cube(b, snap_cube_size)) {
	put_node(b);
	cout << flush;
    }
}

// Set up a template initial state.  We don't actually use all the
// elements of the initial_state3 structure, but we will clean this
// up later if necessary...

void make_standard_bound_init(initial_state3 & init)
{
    init.m2 = 1.0/3;   	        // mass of second particle
    init.m3 = 1.0/3;   		// mass of third particle (m1 + m2 + m3 = 1)
    init.r1 = 0;   		// radius of star #1
    init.r2 = 0;   		// radius of star #2
    init.r3 = 0;   		// radius of star #3

    init.r_stop
          = VERY_LARGE_NUMBER;  // final separation
    init.tidal_tol_factor = DEFAULT_TIDAL_TOL_FACTOR;
    init.eta
          = DEFAULT_ETA;  	// accuracy parameter

    initialize_bodies(init.system);
    init.id = get_initial_seed() + get_n_rand();
}

main(int argc, char **argv)
{
    initial_state3 init;
    make_standard_bound_init(init);

    int  seed 	    = 0;    	// seed for random number generator
    int n_rand      = 0;        // number of times to invoke the generator
                                // before starting for real
    int  n_experiments = 1;     // default: only one run
    real dt_out     =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	// output time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int planar_flag = 0;
    real psi = 0;

    bool  b_flag = FALSE;
    bool  q_flag = FALSE;
    bool  Q_flag = FALSE;
    bool  i_flag = true;
    bool  s_flag = false;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:bc:C:d:D:e:g:L:m:M:n:N:o:pPqQr:R:s:S:U:v:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'b': b_flag = 1 - b_flag;
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      if (dt_out < 0) dt_out = pow(2.0, dt_out);
		      break;
	    case 'D': dt_snap = atof(poptarg);
	    	      if (dt_snap < 0) dt_snap = pow(2.0, dt_snap);
		      break;
	    case 'e': init.ecc = atof(poptarg);
		      break;
	    case 'g': init.tidal_tol_factor = atof(poptarg);
		      break;
	    case 'i': i_flag = !i_flag;
	    	      break;
	    case 'L': init.r_init_min = atof(poptarg);
		      break;
	    case 'm': init.m2 = atof(poptarg);
		      break;
	    case 'M': init.m3 = atof(poptarg);
		      break;
	    case 'n': n_experiments = atoi(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'p': planar_flag = 1;
		      break;
	    case 'P': planar_flag = -1;
		      break;
	    case 'q': q_flag = 1 - q_flag;
		      break;
	    case 'Q': Q_flag = 1 - Q_flag;
		      break;
	    case 'r': init.rho = atof(poptarg);
		      break;
	    case 'R': init.r_stop = atof(poptarg);
		      init.r_init_min = init.r_init_max = abs(init.r_stop);
		      break;
	    case 's': seed = atoi(poptarg);
	    	      s_flag = true;
		      break;
	    case 'S': init.r_stop = atof(poptarg);
		      break;
	    case 'U': init.r_init_max = atof(poptarg);
		      break;
	    case 'v': init.v_inf = atof(poptarg);
		      break;
	    case 'x': init.r1 = atof(poptarg);
		      break;
	    case 'y': init.r2 = atof(poptarg);
		      break;
	    case 'z': init.r3 = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
	}            

    if (Q_flag) q_flag = TRUE;

    if (init.m2 + init.m3 > 1) {
	cerr << "bound3: "; PRC(init.m2); PRL(init.m3);
	exit(1);
    }

    check_init(init);

    // Create the 3-body tree for all calculations.

    sdyn3 *b, *by, *bo;
    b = new sdyn3();
    bo = new sdyn3();
    if (i_flag) bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    // Default from Simo:

    bo->set_mass(1);
    bo->set_pos(vector(-0.97000436, 0.24308753, 0));
    bo->set_vel(-0.5*vector(0.93240737, 0.86473146, 0));

    for (int i = 1; i < 3; i++) {
        by = new sdyn3();
	if (i_flag) by->set_label(i+1);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
	by->set_parent(b);

	if (i == 1) {
	    by->set_mass(1);
	    by->set_pos(-bo->get_pos());
	    by->set_vel(bo->get_vel());
	} else {
	    by->set_mass(1);
	    by->set_pos(0);
	    by->set_vel(-2*bo->get_vel());
	}
        bo = by;
    }

    b->log_history(argc, argv);

    if (s_flag == FALSE) seed = 0;
    int random_seed = srandinter(seed, n_rand);

    cpu_init();

    for (int i = 0; i < n_experiments; i++) {

	if (n_experiments > 1) cerr << i+1 << ": ";

	cerr << "Random seed = " << get_initial_seed()
	     << "  n_rand = " << get_n_rand() << flush;
	init.id = get_initial_seed() + get_n_rand();

	intermediate_state3 inter;
	final_state3 final;

	real cpu = cpu_time();	
	bound3(b, init, inter, final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size);
	cpu = cpu_time() - cpu;

	cerr << ":  ";
//	print_bound3_outcome(inter, final, cerr);
	print_final(cerr, final);
    }
}

#endif
