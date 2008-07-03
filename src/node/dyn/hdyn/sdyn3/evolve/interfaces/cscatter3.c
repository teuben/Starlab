
/*
 * cscatter3.c:  C interface to three-body scattering experiments.
 *		 The "c_" routines are C-callable links to Starlab
 *		 C++ functions.
 */

#define   C_ONLY
#include "scatter3.h"
#include "c_interface.h"

void main(int argc, char **argv)
{
    int  seed 	    = 0;    	/* seed for random number generator */
    int n_rand      = 0;        /* number of times to invoke the generator */
                                /* before starting for real */
    int  n_experiments = 1;     /* default: only one run */
    real dt_out     =       	/* output time interval */
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	/* output time interval */
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int  planar_flag = 0;
    bool psi_flag = FALSE;
    real psi = 0;

    bool b_flag = FALSE;
    bool q_flag = FALSE;
    bool Q_flag = FALSE;

    extern char *poptarg;
    int c;
    const char *param_string = "A:bc:C:d:D:e:g:L:m:M:n:N:o:pPqQr:R:s:S:U:v:x:y:z:";

    int random_seed;
    real cpu;
    int i;

    initial_state3 init;
    intermediate_state3 inter;
    final_state3 final;

    /* Initialize the initial state structure. */

    c_make_standard_init(&init);

    while ((c = c_pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'b': b_flag = 1 - b_flag;
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg); /* (in hours) */
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      break;
	    case 'D': dt_snap = atof(poptarg);
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
		      psi_flag = TRUE;
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
            case '?': /*fprintf(stderr, "Usage: %s\n", param_string);*/
		      exit(1);
	}            

    if (Q_flag) q_flag = TRUE;

    if (init.m2 > 1) {
/*	fprintf(stderr, "cscatter3:  init.m2 = %.5f > 1\n", init.m2);*/
	exit(1);
    }

    c_cpu_init();
    random_seed = c_srandinter(seed, n_rand);

    for (i = 0; i < n_experiments; i++) {

	if (n_experiments > 1) fprintf(stderr, "%d: ", i+1);
	c_print_initial_random_parameters();

	c_initialize_angles(&init, planar_flag, psi_flag, psi);

	cpu = c_cpu_time();	
	c_scatter3(&init, &inter, &final, cpu_time_check,
		   dt_out, dt_snap, snap_cube_size);
	cpu = c_cpu_time() - cpu;

	c_print_scatter3_info(&init, &inter, &final,
			      Q_flag, q_flag, b_flag, cpu);
    }
}

