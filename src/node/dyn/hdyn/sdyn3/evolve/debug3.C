
// debug.C: Perform three-body scattering experiments, with diagnostic output.
//
// Input is an initial state, output is intermediate and final states.
// No reference is made to scatter profiles in this file.

// Starlab application:  scatter3.

#include "scatter3.h"

local real energy(sdyn3 * b)
{
    real k = 0, u = 0;

    for_all_daughters(sdyn3, b, bi) {
	k += bi->get_mass() * bi->get_vel() * bi->get_vel();
	for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	  u -= bi->get_mass() * bj->get_mass()
	       / abs(bi->get_pos() - bj->get_pos());
    }
    return 0.5*k + u;
}

local void print_debug(sdyn3 *s)

// Print out information on a scattering experiment.  Control is passed
// to this function at time intervals of dt_print.  Note that, for a flyby,
// the first two sdyn3s below s are the components of the binary.

{
/*
    int p = cout.precision(LOW_PRECISION);
    cout.setf(ios::scientific, ios::floatfield);

    sdyn3 *d1, *d2, *d3;

    d1 = s->get_oldest_daughter();
    d2 = d1->get_younger_sister();
    d3 = d2->get_younger_sister();

    real m1 = d1->get_mass();
    real m2 = d2->get_mass();
    real m3 = d3->get_mass();

    kepler k;

    k.set_time(s->get_time());
    k.set_total_mass(d1->get_mass() + d2->get_mass());
    k.set_rel_pos(d2->get_pos() - d1->get_pos());
    k.set_rel_vel(d2->get_vel() - d1->get_vel());
    k.initialize_from_pos_and_vel();
 */

    cout << "ssd = " << s->get_ssd()
	 << "  n_osc = " << s->get_n_ssd_osc();
    cout << "  T = " << s->get_time() + s->get_time_offset();

//    cout << "  CPU = " << cpu_time() << endl;
    cout << "  E = " << energy(s) << endl;

    cout.precision(p);
}

// The main program is identical to scatter3, except that it takes an
// extra command-line flag (dt_print) and passes print_debug as an
// argument to the scatter3 function.

main(int argc, char **argv)
    {

    initial_state3 init;
    make_standard_init(init);

    int  seed 	    = 0;    	// seed for random number generator
    int  n_rand = 0;
    int  n_experiments = 1;     // default: only one run
    real dt_out     =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_print    =       	// diagnostic print time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    bool  b_flag = FALSE;
    bool  q_flag = FALSE;

    extern char *poptarg;
    int c;
    char* param_string = "A:bc:C:d:D:e:L:m:M:nN::qp:r:s:U:v:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'b': b_flag = 1 - b_flag;
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      break;
	    case 'D': dt_snap = atof(poptarg);
		      break;
	    case 'e': init.ecc = atof(poptarg);
		      break;
	    case 'L': init.r_init_min = atof(poptarg);
		      break;
	    case 'm': init.m2 = atof(poptarg);
		      break;
	    case 'M': init.m3 = atof(poptarg);
		      break;
            case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'n': n_experiments = atoi(poptarg);
		      break;
	    case 'p': dt_print = atof(poptarg);
		      break;
	    case 'q': q_flag = 1 - q_flag;
		      break;
	    case 'r': init.rho = atof(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
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
		      exit(1);
	    }            

    if (init.m2 > 1)
	{
	cerr << "sigma3: initial.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
	}

    int random_seed = srandinter(seed, n_rand);

    while (n_experiments--)
	{
	int initial_seed = get_initial_seed();
	int n_rand = get_n_rand();

	randomize_angles(init.phase);

	intermediate_state3 inter;
	final_state3 final;

	scatter3(init, inter, final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size,
		 dt_print, print_debug);

	cerr << "result:  " << state_string(inter.descriptor) << " "
	     << state_string(final.descriptor)
	     << "        Random seed = " << initial_seed
	     << "  n_rand = " << n_rand << endl;

	if (!q_flag)
	    {
	    print_initial(cerr, init, b_flag);
	    print_intermediate(cerr, inter, b_flag);
	    print_final(cerr, final, b_flag);
	    }
	}
    }

