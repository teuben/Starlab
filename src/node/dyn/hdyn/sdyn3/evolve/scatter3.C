
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// scatter3:  Perform three-body scattering experiments.
////
////            The program scatter3 randomizes all phase angles prior
////            to invoking the scatter3 function.  Other parameters are
////            fully determined, either by default or on the command line.
////
////            The program initializes a binary-single star scattering
////            at a separation determined by the tidal tolerance, and
////            by default integrates the system until a clear dynamical
////            outcome is recognized (escaping star or three unbound
////            stars).
////
////            Default program output is initial, intermediate and final
////            scattering states.
////
////            Units: G = 1, binary mass = 1, binary semi-major axis = 1.
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
////            -n    specify number of scattering experiments [1]
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

// No reference is made to scatter profiles in this file.

// Note the function scatter3 is completely deterministic.
// No randomization is performed at this level.

// Starlab library function.

#include "scatter3.h"

#ifndef TOOLBOX

// init_to_sdyn3:  Convert an initial state to an sdyn3 for integration.

local sdyn3 * init_to_sdyn3(initial_state3 & init, final_state3 & final)
{
    set_kepler_tolerance(2);

    kepler k1;		// inner binary.

    real peri = 1;	// default value (unimportant unless init.ecc = 1).
    if (init.ecc == 1) peri = 0;

    make_standard_kepler(k1, 0, 1, -0.5 / init.sma, init.ecc,
			 peri, init.phase.mean_anomaly, 1);  // "1" ==> align
							     // with x, y axes

#if 0
    cout << endl << "Inner binary:" << endl;
    k1.print_all(cout);
#endif

    // Multiply the incoming velocity by vcrit.

    real m2 = init.m2;
    real m3 = init.m3;
    real mtotal = 1 + m3;
    real v_inf = init.v_inf * sqrt( (1 - m2) * m2 * mtotal / m3 );
    
    real energy3 = .5 * v_inf * v_inf;
    real ang_mom3 = init.rho * v_inf;

    real ecc3 = (energy3 == 0 ? 1
		              : sqrt( 1 + 2 * energy3 
				            * pow(ang_mom3/mtotal, 2)));

    // Special case: if v_inf = 0, assume rho is periastron.

    real virial_ratio = init.rho;

    kepler k3;		// outer binary.

    // NOTE: passing a mean anomaly of 0 here will cause problems on
    // Linux systems if the orbit is linear (i.e. rho = 0), as this
    // implies r = 0, v = infinity.  Most systems will accept v = Inf
    // (which is actually OK, as we later set the separation to something
    // more reasonable).  However, Linux fails with a floating excetion.
    // Any nonzero value should be OK, as far as initialization is
    // concerned, but any given value of mean_anomaly will present
    // difficulties as v_inf --> 0 (because the corresponding separation
    // will be outside the desired value, leading to problems with the
    // function return_to_radius() below).  Use a negative value of
    // mean_anomaly here to keep us on the incoming branch, and take
    // special precautions later.

    real mean_anomaly = -0.01;
    make_standard_kepler(k3, 0, mtotal, energy3, ecc3, virial_ratio,
			 mean_anomaly);

    // Radius for "unperturbed" inner binary:

    real r_unp = (init.sma + init.r3)
	             * pow(init.tidal_tol_factor / mtotal, -1/3.0);

    if (r_unp <= k3.get_periastron()) {
	final.descriptor = preservation;
	final.sma = k1.get_semi_major_axis();
	final.ecc = k1.get_eccentricity();
	final.outer_separation = k3.get_periastron();
	final.escaper = 3;
	final.error = 0;
	final.time = 0;
	final.n_steps = 0;
	final.virial_ratio = k3.get_energy() * k3.get_periastron() 
	                                        / k3.get_total_mass() - 1;

	// Sufficiently large stellar radii (sum larger than the binary
	// periastron) will cause an immediate collision to occur, so we
	// must check for it explicitly here.

	if (k1.get_periastron() < init.r1 + init.r2) {
	    final.descriptor = merger_escape_3;
	    final.sma = final.ecc = -1;
	}

	return NULL;
    }

    init.r_init = max(init.r_init_min, min(init.r_init_max, r_unp));

    // NOTE: Can't assume that k3.separation is less than r_init.

    if (k3.get_separation() < init.r_init)
        k3.return_to_radius(init.r_init);
    else
        k3.advance_to_radius(init.r_init);

    k1.transform_to_time(k3.get_time());
    set_orientation(k3, init.phase);

#if 0
    cout << endl << "Outer binary:" << endl << flush;
    k3.print_all(cout);
#endif

    sdyn3 * b;
    b = set_up_dynamics(m2, m3, k1, k3);

//  Set up radii:

    sdyn3 * bb = b->get_oldest_daughter();
    bb->set_radius(init.r1);
    bb = bb->get_younger_sister();
    bb->set_radius(init.r2);
    bb = bb->get_younger_sister();
    bb->set_radius(init.r3);

    return  b;
}

// check_init:  Make sure an initial state is sensible.

local void check_init(initial_state3 & init)
{
    if (init.m2 > 1) err_exit("check_init: m1 < 0");
    if (init.m2 < 0) err_exit("check_init: m2 < 0");
    if (init.m3 < 0) err_exit("check_init: m3 < 0");
    if (init.sma <= 0) err_exit("check_init: semi-major axis <= 0");
    if (init.ecc < 0) err_exit("check_init: eccentricity < 0");
    if (init.ecc > 1) err_exit("check_init: eccentricity > 1");
    if (init.r1 < 0) err_exit("check_init: r1 < 0");
    if (init.r2 < 0) err_exit("check_init: r2 < 0");
    if (init.r3 < 0) err_exit("check_init: r3 < 0");
    if (init.v_inf < 0) err_exit("check_init: v_inf < 0");
    if (init.rho < 0) err_exit("check_init: rho < 0");
    if (init.eta <= 0) err_exit("check_init: eta <= 0");
}

// scatter3:  Perform a three-body scattering, initializing from a specified
//            state and returning intermediate- and final-state structures.

void scatter3(initial_state3 & init,
	      intermediate_state3 & inter,
	      final_state3 & final,
	      real cpu_time_check,
	      real dt_out,         // diagnostic output interval
	      real dt_snap,        // snapshot output interval
	      real snap_cube_size, // limit output to particles within cube
	      real dt_print,       // print output interval
	      sdyn3_print_fp p)    // function to print output
{
    check_init(init);
    inter.id = final.id = init.id;	    // For bookkeeping...

    sdyn3 * b = init_to_sdyn3(init, final); // Note that all systems are
    initialize_bodies(inter.system);	    // initialized as "non-particles"
    initialize_bodies(final.system);	    // until explicitly set.

    if (!b) {

	// Initial separation was less than pericenter, so just set
        // set up the intermediate state and quit.

        inter.n_osc = 1;
	inter.n_kepler = 0;
	inter.r_min[0] = init.sma * (1 - init.ecc);
	inter.r_min[1] = inter.r_min[0];
	inter.r_min[2] = final.sma * (final.ecc - 1); // (Approximately)
	inter.r_min_min = inter.r_min[0];
	inter.descriptor = non_resonance;

	return;
    }

    // Save some initial data:

    sdyn3_to_system(b, init.system);

    real e_init = energy(b);
    real cpu_init = cpu_time();
    real cpu = cpu_init;
    
    // Evolve the system until the interaction is over or a collision occurs.

    inter.n_kepler = 0;

    do {

	int n_stars = 0;
	for_all_daughters(sdyn3, b, bb) n_stars++;

	tree3_evolve(b, CHECK_INTERVAL, dt_out, dt_snap, snap_cube_size,
		     init.eta, cpu_time_check, dt_print, p);

	// Check the CPU time.  Note that the printed CPU time is the
	// time since this routine was entered.

	if (cpu_time() - cpu > cpu_time_check) {
	    if (cpu == cpu_init) cerr << endl; // (May be waiting in mid-line!)

	    int p = cerr.precision(STD_PRECISION);
	    cerr << "scatter3:  CPU time = " << cpu_time() - cpu_init
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
		put_node(cout, *b);
		cout << flush;
	    }
	}

    } while (!extend_or_end_scatter(b, init, &inter, &final));

    // Set intermediate state:

    inter.r_min_min = VERY_LARGE_NUMBER;
    inter.n_stars = 0;

    for_all_daughters(sdyn3, b, bb) {
        inter.index[inter.n_stars] = bb->get_index();
        inter.r_min[inter.n_stars] = sqrt(bb->get_min_nn_dr2());
	inter.r_min_min = min(inter.r_min_min, inter.r_min[inter.n_stars]);
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
	put_node(cout, *b);
	cout << flush;
    }

    // Delete the 3-body sytem.

    sdyn3 * bi = b->get_oldest_daughter();
    while (bi) {
	sdyn3 * tmp = bi->get_younger_sister();
	delete bi;
	bi = tmp;
    }
    delete b;
}

#else

main(int argc, char **argv)
{
    initial_state3 init;
    make_standard_init(init);

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
    int psi_flag = 0;
    real psi = 0;

    bool  b_flag = FALSE;
    bool  q_flag = FALSE;
    bool  Q_flag = FALSE;

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
	    case 'o': psi = atof(poptarg);
		      psi_flag = 1;
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

    if (init.m2 > 1) {
	cerr << "scatter3: init.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
    }

    cpu_init();
    int random_seed = srandinter(seed, n_rand);

    for (int i = 0; i < n_experiments; i++) {

	if (n_experiments > 1) cerr << i+1 << ": ";

	cerr << "Random seed = " << get_initial_seed()
	     << "  n_rand = " << get_n_rand() << flush;
	init.id = get_initial_seed() + get_n_rand();

	// Normally, we just want random angles and phases.
	// However, in the planar case, it may be desireable to
	// control the orientation of the outer orbit.  The
	// angle between the inner and outer orbital axes in
	// the planar case is init.phase.psi.  Specify psi
	// in *degrees*, with psi = 0 meaning that the outer
	// periastron lies along the positive x-axis.

	randomize_angles(init.phase);

	if (planar_flag == 1) {
	    init.phase.cos_theta = 1;	// Planar prograde
	    if (psi_flag) init.phase.psi = psi * M_PI / 180.0;
	} else if (planar_flag == -1) {
	    init.phase.cos_theta = -1;	// Planar retrograde
	    if (psi_flag) init.phase.psi = psi * M_PI / 180.0;
	}

	intermediate_state3 inter;
	final_state3 final;

	real cpu = cpu_time();	
	scatter3(init, inter, final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size);
	cpu = cpu_time() - cpu;

	cerr << ":  ";
	print_scatter3_outcome(inter, final, cerr);

	if (Q_flag) print_scatter3_summary(inter, final, cpu, cerr);

	if (!q_flag) print_scatter3_report(init, inter, final,
					   cpu, b_flag, cerr);
    }
}

#endif
