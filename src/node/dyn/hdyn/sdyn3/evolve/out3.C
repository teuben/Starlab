
// out3.C: Perform a series of three-body scattering experiments
//	   and print out the result of each.

// Starlab application:  scatter3.

#include "scatter3.h"

local void print(sdyn3* b)
    {
    fprintf(stdout, "%.5f ", b->get_time() + b->get_time_offset());
    for_all_daughters(sdyn3, b, bb)
	for (int k = 0; k < 3; k++) fprintf(stdout, "%.5f ", bb->get_pos()[k]);
    fprintf(stdout, "\n");
    }

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
    real dt_print    =       	// print time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int planar_flag = 0;

    bool  b_flag = FALSE;
    bool  q_flag = FALSE;
    bool  Q_flag = FALSE;

    extern char *poptarg;
    int c;
    char* param_string = "A:bc:C:d:D:e:L:m:M:n:N:pPqQr:R:s:S:U:v:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
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
	    case 'D': dt_print = atof(poptarg);
		      break;
	    case 'e': init.ecc = atof(poptarg);
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
	    case 'R': init.r_stop = init.r_init_min
				  = init.r_init_max
				  = atof(poptarg);
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
		      exit(1);
	    }            

    if (Q_flag) q_flag = TRUE;

    if (init.m2 > 1)
	{
	cerr << "sigma3: init.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
	}

    cpu_init();
    int random_seed = srandinter(seed, n_rand);

    while (n_experiments--)
	{
	cerr << "Random seed = " << get_initial_seed()
	     << "  n_rand = " << get_n_rand() << flush;

	randomize_angles(init.phase);

	if (planar_flag == 1)
	    init.phase.cos_theta = 1;	// Planar prograde
	else if (planar_flag == -1)
	    init.phase.cos_theta = -1;	// Planar retrograde

	intermediate_state3 inter;
	final_state3 final;

	real cpu = cpu_time();	
	scatter3(init, inter, final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size,
		 dt_print, print);
	cpu = cpu_time() - cpu;

	cerr << ":  ";
	print_scatter3_outcome(inter, final, cerr);

	if (Q_flag) print_scatter3_summary(inter, final, cpu, cerr);

	if (!q_flag) print_scatter3_report(init, inter, final, b_flag, cerr);

	}
    }
