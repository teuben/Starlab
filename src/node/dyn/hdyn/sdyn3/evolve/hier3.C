
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// hier3:    Study the dynamics of hierarchical 3-body systems.
////
////            All phase angles set prior to invoking hier3.  Inner binary
////            starts at random phase; outer binary at apocenter.  In non-
////            planar case, outer orbit orientation is random.
////
////            Other parameters are fully determined, either by default
////            or on the command line.
////
////            The program initializes a hierarchical triple system and
////            follows its dynamics until a clear outcome is recognized
////            (escaping star).
////
////            Units: G = 1, binary mass = 1, binary semi-major axis = 1.
////
//// Options:   -a    specify outer semi-major axis [10]
////            -A    specify accuracy parameter [0.05]
////            -c    specify CPU time check, in hours [1]
////            -C    specify snap cube size [10]
////            -d    specify log output interval [none]
////            -D    specify snap output interval [none]
////            -e    specify initial inner eccentricity [0]
////            -E    specify initial outer eccentricity [0]
////            -g    specify tidal factor for extension/termination [1.e-6]
////            -m    specify secondary mass (binary mass = 1) [0.5]
////            -M    specify incomer mass [0.5]
////            -n    specify number of experiments [1]
////            -N    specify random number count [0]
////            -p    force planar prograde [not forced]
////            -P    force planar retrograde [not forced]
////            -q    specify outer pericenter rather than eccentricity [no]
////            -Q    suppress most diagnostic output [don't suppress]
////            -R    specify stopping radius [100]
////            -s    specify random seed [taken from system clock]
////            -S    specify separation at which to stop [none]
////            -t    specify output interval [none]
////            -T    specify log10 duration of run, in units of the outer
////                  period (Eggleton & Kiseleva's "n") [2]
////            -x    specify primary radius [0]
////            -X    specify Eggleton & Kiseleva period ratio X [none]
////            -y    specify secondary radius [0]
////            -X    specify Eggleton & Kiseleva periastron ratio Y [none]
////            -z    specify incomer radius [0]
////
//// As a convenient shorthand, any "dt" interval specified less than zero
//// is interpreted as a power of 2, i.e. "-d -3" sets dt_out = 0.125.

// Note the function hier3 is completely deterministic.
// No randomization is performed at that level.

#include "scatter3.h"

#ifndef TOOLBOX

// No library functions.

#else

// make_hier_init:  Set up a template initial state.

#define DEFAULT_A_OUT	10
#define DEFAULT_R_STOP 100

local void make_hier_init(initial_state3 & init)
{
    init.m2 = 0.5;   	        // mass of secondary in target binary
    init.m3 = 0.5;   		// mass of projectile (m1 + m2 = 1)
    init.r1 = 0;   		// radius of star #1
    init.r2 = 0;   		// radius of star #2
    init.r3 = 0;   		// radius of star #3

    init.sma = 1;		// inner binary semi-major axis
    init.ecc = 0;    		// inner binary eccentricity

    // Not used (for unbound systems):

    init.v_inf = 0;
    init.rho = 0;

    init.a_out = DEFAULT_A_OUT;	// outer binary semi-major axis
    init.e_out = 0;    		// outer binary eccentricity
    init.r_stop = DEFAULT_R_STOP;
				// final separation
    init.tidal_tol_factor = DEFAULT_TIDAL_TOL_FACTOR;
    init.eta
          = DEFAULT_ETA;  	// accuracy parameter

    initialize_bodies(init.system);
    init.id = get_initial_seed() + get_n_rand();
}

// init_to_sdyn3:  Convert an initial state to an sdyn3 for integration.

local sdyn3 * init_to_sdyn3(initial_state3 & init)
{
    set_kepler_tolerance(2);

    kepler k1;		// inner binary

    real mean_anomaly = 0;
    if (init.ecc == 1)
	mean_anomaly = -M_PI;
    else
	mean_anomaly = randinter(-M_PI, M_PI);

    make_standard_kepler(k1, 0, 1, -0.5 / init.sma, init.ecc,
			 1,	// value unimportant unless init.ecc = 1
			 mean_anomaly,
			 1);	// 1 ==> align with x, y axes

#if 0
    cerr << endl << "Inner binary:" << endl;
    k1.print_all(cerr);
#endif

    if (init.e_out == 1)
	err_exit("Linear outer binary not allowed!");

    kepler k3;		// outer binary.

    real m2 = init.m2;
    real m3 = init.m3;
    real mtotal = 1 + m3;

    make_standard_kepler(k3, 0, mtotal, -0.5 * mtotal / init.a_out, init.e_out,
			 1,	// value unimportant unless init.e_out = 1
			 M_PI);	// start at apastron
    set_orientation(k3, init.phase);

#if 0
    cerr << endl << "Outer binary:" << endl << flush;
    k3.print_all(cerr);
#endif

    sdyn3 * b;
    b = set_up_dynamics(m2, m3, k1, k3);

    // Set up radii:

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
    if (init.sma <= 0) err_exit("check_init: inner semi-major axis <= 0");
    if (init.ecc < 0) err_exit("check_init: inner eccentricity < 0");
    if (init.ecc > 1) err_exit("check_init: inner eccentricity > 1");
    if (init.a_out <= 0) err_exit("check_init: outer semi-major axis <= 0");
    if (init.e_out < 0) err_exit("check_init: outer eccentricity < 0");
    if (init.e_out > 1) err_exit("check_init: outer eccentricity > 1");
    if (init.r1 < 0) err_exit("check_init: r1 < 0");
    if (init.r2 < 0) err_exit("check_init: r2 < 0");
    if (init.r3 < 0) err_exit("check_init: r3 < 0");
    if (init.eta <= 0) err_exit("check_init: eta <= 0");
}

local void print_elements(sdyn3* b)
{
    // Compute and print out the essential orbital elements of the
    // assumed ((b1,b2),b3) hierarchical system.  Very inefficient!

    sdyn3* b1 = b->get_oldest_daughter();
    sdyn3* b2 = b1->get_younger_sister();
    sdyn3* b3 = b2->get_younger_sister();

    // Create a center of mass node.

    sdyn3 cm;
    cm.set_mass(b1->get_mass()+b2->get_mass());
    cm.set_pos((b1->get_mass()*b1->get_pos()+b2->get_mass()*b2->get_pos())
	        / cm.get_mass());
    cm.set_vel((b1->get_mass()*b1->get_vel()+b2->get_mass()*b2->get_vel())
	        / cm.get_mass());

    // Construct keplers describing the inner and outer orbits.

    kepler k1, k3;
    initialize_kepler_from_dyn_pair(k1, b1, b2, true);
    initialize_kepler_from_dyn_pair(k3, b3, &cm, true);

    // Print out vital information.

    cerr << b->get_system_time()      << "  "
	 << k1.get_semi_major_axis()  << "  "
	 << k3.get_semi_major_axis()  << "  "
	 << k1.get_eccentricity()     << "  "
	 << k3.get_eccentricity()     << "  "
	 << 0.1 * k3.get_periastron() << "  "
	 << 0.1 * k3.get_periastron() / k1.get_apastron() << "  "
	 << endl;
}

#define PERT_TOL 0.5

local real period_ratio(sdyn3* a, sdyn3* b, sdyn3* c, bool debug = false)
{
    // Return the ratio of the squares of the periods of the inner
    // and outer orbits, if both are bound.  Return 0 otherwise.

    if (debug) {
	cerr << "period_ratio:  a = " << a->format_label();
	cerr << "  b = " << b->format_label();
	cerr << "  c = " << c->format_label() << endl;
    }

    // First determine the relative energy of the a-b motion.

    real ma = a->get_mass();
    real mb = b->get_mass();
    real mc = c->get_mass();

    real mab = ma + mb;
    vector rab = b->get_pos() - a->get_pos();
    vector vab = b->get_vel() - a->get_vel();

    real dab = abs(rab);
    real eab = 0.5*square(vab) - mab/dab;
    if (debug) {
	PRI(4); PRC(dab); PRC(eab);
    }

    if (eab >= 0) {
	if (debug) cerr << endl;
	return 0;
    }

    // Now look at the "outer" pair ((a,b),c).

    // Determine the position and velocity of the CM of a and b.

    vector cmpos = (ma*a->get_pos() + mb*b->get_pos()) / mab;
    vector cmvel = (ma*a->get_vel() + mb*b->get_vel()) / mab;

    real mabc = mab + mc;
    vector rc = c->get_pos() - cmpos;
    vector vc = c->get_vel() - cmvel;

    // EXTRA CONDITION: require that c not be a dominant perturber of (a,b).

    real dc = abs(rc);

    // Multiplies here because of annoying pow bug in RH Linux...

    if (mc*dab*dab*dab > PERT_TOL*mab*dc*dc*dc) return 0;

    real ec = 0.5*square(vc) - mabc/dc;
    if (debug) {
	PRC(abs(rc)), PRL(ec);
    }

    if (ec >= 0) return 0;

    // Now we know that both the inner and the outer pairs are bound.
    // Compute the squared period ratio.

    real ee = eab/ec;	// defeat annoying pow bug...

    return (mabc/mab) * ee * ee;
}

local int outer_component(sdyn3* b, bool debug = false)
{
    // Determine the index of the outer component of the "most hierarchical"
    // system, according to the criterion of Eggleton & Kiseleva (1995).

    // Hmmm...  Seems that this can fail quite easily, giving a large period
    // to a spurious "outer" system and hence identifying the wrong star as
    // the outermost member (Steve 5/99).

    // E&K criterion is simply to maximize the ratio of outer to inner
    // periods.  Probably better to incorporate the perturbation due to
    // the third member as part of the criterion (Steve 5/99).

    sdyn3* b1 = b->get_oldest_daughter();
    sdyn3* b2 = b1->get_younger_sister();
    sdyn3* b3 = b2->get_younger_sister();

    // Test all possible combinations of inner and outer binaries.

    real ratio1 = period_ratio(b2, b3, b1, debug);
    real ratio2 = period_ratio(b3, b1, b2, debug);
    real ratio3 = period_ratio(b1, b2, b3, debug);

    if (debug) {
	cerr << "outer_component:  ";
	PRC(ratio1), PRC(ratio2), PRL(ratio3);
    }

    real  ratio_max = max(max(ratio1, ratio2), ratio3);

    if (ratio_max <= 1)
	return 0;
    else if (ratio_max == ratio1)
	return 1;
    else if (ratio_max == ratio2)
	return 2;
    else
	return 3;
}

// hier3:  Perform integration of a hirarchical three-body system,
//         initializing from a specified state.

local void hier3(initial_state3 & init,
		 real cpu_time_check,
		 real dt_out,         // diagnostic output interval
		 real dt_snap,        // snapshot output interval
		 real snap_cube_size, // limit output to particles within cube
		 real dt_print,	      // print output interval
		 real t_end,
		 bool quiet = false)
{
    check_init(init);
    sdyn3 * b = init_to_sdyn3(init);

    if (!b) err_exit("Unable to initialize system...");

    // Save some initial data:

    sdyn3_to_system(b, init.system);

    real e_init = energy(b);
    real cpu_init = cpu_time();
    real cpu = cpu_init;
    
    // Evolve the system until the interaction is over or a collision occurs.

    int n_kepler = 0;

    do {

	int n_stars = 0;
	for_all_daughters(sdyn3, b, bb) n_stars++;

	tree3_evolve(b, CHECK_INTERVAL, dt_out, dt_snap, snap_cube_size,
		     init.eta, cpu_time_check,
		     dt_print, print_elements);

	// Stop immediately if hierarchy has changed.

	if (outer_component(b) != 3) break;

	// Check the CPU time.  Note that the printed CPU time is the
	// time since this routine was entered.

	if (cpu_time() - cpu > cpu_time_check) {

	    if (!quiet) {

		if (cpu == cpu_init)
		    cerr << endl; // (May be waiting in mid-line!)

		int p = cerr.precision(STD_PRECISION);
		cerr << "hier3:  CPU time = " << cpu_time() - cpu_init
		     << "  time = " << b->get_system_time()
		     << "  n_steps = " << b->get_n_steps()
		     << endl;
		cerr << "           n_osc = " << b->get_n_ssd_osc()
		     << "  n_kepler = " << n_kepler
		     << endl << flush;
		cerr.precision(p);
	    }
	    cpu = cpu_time();
	}

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

    } while (b->get_system_time() < t_end &&
	     !extend_or_end_scatter(b, init));	// this is the same termination
						// condition as for scattering
						// experiments -- sufficient,
						// but probably unnecessarily
						// stringent.

    // In all cases, the final system should be a properly configured
    // 1-, 2-, or 3-body system.  Print it for the last time, to reflect
    // final state, if the "-D" command-line option was specified.

    if (dt_snap < VERY_LARGE_NUMBER && system_in_cube(b, snap_cube_size)) {
	put_node(b);
	cout << flush;
    }

    int p = cerr.precision(STD_PRECISION);
    cerr << "t = " << b->get_system_time()
	 << "  t/t_end = " << b->get_system_time()/t_end
	 << "  outer = " << outer_component(b)
	 << endl;
    cerr.precision(p);

    if (!quiet)	outer_component(b, true);

    // Delete the 3-body sytem.

    sdyn3 * bi = b->get_oldest_daughter();
    while (bi) {
	sdyn3 * tmp = bi->get_younger_sister();
	delete bi;
	bi = tmp;
    }
    delete b;
}

main(int argc, char **argv)
{
    initial_state3 init;
    make_hier_init(init);

    int  seed 	    = 0;    	// seed for random number generator
    int n_rand      = 0;        // number of times to invoke the generator
                                // before starting for real
    int  n_experiments = 1;     // default: only one run
    real dt_out     =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_print   =        	// print time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int planar_flag = 0;
    int psi_flag = 0;
    real psi = 0;

    real outer_peri = 0;
    bool peri_flag = false;

    bool quiet = false;

    bool X_flag = false, Y_flag = false;    // flags for using Eggleton
					    // and Kiseleva parameters

    real X_KE = 10, Y_KE = 4;		    // set above flags true to make
					    // these the defaults

    real T_stab = 100;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string
	= "a:A:c:C:d:D:e:E:g:m:M:n:N:o:pPq:QR:s:S:t:T:x:X:y:Y:z:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'a': init.a_out = atof(poptarg);
	    	      X_flag = false;
		      break;
	    case 'A': init.eta = atof(poptarg);
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
	    	      peri_flag = false;
		      break;
	    case 'E': init.e_out = atof(poptarg);
		      peri_flag = false;
		      Y_flag = false;
		      break;
	    case 'g': init.tidal_tol_factor = atof(poptarg);
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
	    case 'q': outer_peri = atof(poptarg);
		      peri_flag = true;
	    	      Y_flag = false;
		      break;
	    case 'Q': quiet = !quiet;
		      break;
	    case 'R': init.r_stop = atof(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 'S': init.r_stop = atof(poptarg);
		      break;
	    case 't': dt_print = atof(poptarg);
		      if (dt_print < 0) dt_print = pow(2.0, dt_print);
		      break;
	    case 'T': T_stab = pow(10.0, atof(poptarg));
		      break;
	    case 'x': init.r1 = atof(poptarg);
		      break;
	    case 'X': X_KE =atof(poptarg);
		      X_flag = true;
		      break; 
	    case 'y': init.r2 = atof(poptarg);
		      break;
	    case 'Y': Y_KE =atof(poptarg);
	    	      peri_flag = false;
		      Y_flag = true;
		      break; 
	    case 'z': init.r3 = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
	}            

    if (init.m2 > 1) {
	cerr << "hier3: init.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
    }

    if (X_flag) {

	// Convert period ratio X_KE to outer semi-major axis.
	// Inner semi-major axis = 1 by convention.

	if (X_KE <= 1) err_exit("X_KE <= 1");

	init.a_out = pow((1+init.m3)*X_KE*X_KE, 1.0/3.0);
    }

    if (Y_flag) {

	// Convert periastron ratio Y_KE to outer eccentricity.
	// Note that Y_KE is outer periastron / inner apastron.

	if (Y_KE <= 1) err_exit("Y_KE <= 1");

	init.e_out = 1 - (init.sma*(1+init.ecc) * Y_KE) / init.a_out;
    }

    if (peri_flag) {

	// Convert outer periastron to eccentricity.

	if (outer_peri <= 0) err_exit("outer_peri <= 0");

	init.e_out = 1 - outer_peri / init.a_out;;
    }

    cpu_init();
    int random_seed = srandinter(seed, n_rand);

    for (int i = 0; i < n_experiments; i++) {

	if (!quiet) {
	    if (n_experiments > 1) cerr << i+1 << ": ";

	    cerr << "Random seed = " << get_initial_seed()
	 	 << "  n_rand = " << get_n_rand() << endl << flush;
	}

	init.id = get_initial_seed() + get_n_rand();

	// Phase angles will be such that the inner binary is at
	// periastron, the outer binary at apastron.  Orbital
	// orientation will be random, by default.

	// In the planar case, it may be desireable to
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

	hier3(init, cpu_time_check,
	      dt_out, dt_snap, snap_cube_size, dt_print,
	      2*M_PI*T_stab*sqrt(pow(init.a_out, 3)/(1+init.m3)),
	      quiet);
    }
}

#endif
