
//// tree3_evolve:  simpler interface to low_n3_evolve with offset time,
////                dynamic timesteps, and no softening.
////
//// Options:      -A    specify accuracy parameter [0.02]
////               -c    specify CPU time check interval (s) [3600]
////               -C    specify cube size for snapshot output [10]
////               -d    specify log output interval [none]
////               -D    specify snapshot output interval [none]
////               -q    quiet output [verbose]
////               -t    specify time span of integration [1]

//   Note: low_n3_evolve() uses flat trees

// Starlab library function.

#include "scatter3.h"

#ifndef TOOLBOX

// Integration parameters:

#define SOFTENING	     0
#define X_FLAG		     1
#define TIMESTEP_CRITERION  "dynamic_timestep"
#define S_FLAG		     1
#define N_ITER    	     1
#define	N_MAX		    -1

void tree3_evolve(sdyn3 * b,          // sdyn3 array
		 real delta_t,        // time span of the integration
		 real dt_out,         // output time interval
		 real dt_snap,        // snapshot output interval
		 real snap_cube_size, // limit output to particles within cube
		 real eta,            // time step parameter
		 real cpu_time_check,
		 real dt_print,       // external print interval
		 sdyn3_print_fp p)    // pointer to external print function
{
    // Use offset times within the integrator to avoid roundoff.
    // Reset before returning.

    real t_offset = b->get_time();

    b->begin_offset_time(t_offset);
    for_all_daughters(sdyn3, b, bb)
	bb->begin_offset_time(t_offset);

    low_n3_evolve(b, delta_t, dt_out, dt_snap, snap_cube_size,
		  SOFTENING, eta, X_FLAG, TIMESTEP_CRITERION,
		  S_FLAG, N_ITER, N_MAX,
		  cpu_time_check, dt_print, p);

    b->end_offset_time();
    for_all_daughters(sdyn3, b, bbb)
	bbb->end_offset_time();

}

#else

main(int argc, char **argv)
{
    sdyn3* b;            // pointer to the nbody system
    
    real  delta_t = 10;  // time span of the integration
    real  eta = 0.02;    // time step parameter (for fixed time step,
                         //   equal to the time step size; for variable
                         //   time step, a multiplication factor)
    real  dt_out = VERY_LARGE_NUMBER;
                         // output time interval
    real  dt_snap = VERY_LARGE_NUMBER;
                         // snap output interval
    real  snap_cube_size = 10;

    real cpu_time_check = 3600;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:qt:";

    bool  a_flag = FALSE;
    bool  d_flag = FALSE;
    bool  D_flag = FALSE;
    bool  q_flag = FALSE;
    bool  t_flag = FALSE;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'A': a_flag = TRUE;
		      eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = atof(poptarg);
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': d_flag = TRUE;
		      dt_out = atof(poptarg);
		      break;
	    case 'D': D_flag = TRUE;
		      dt_snap = atof(poptarg);
		      break;
	    case 'q': q_flag = TRUE;
		      break;
	    case 't': t_flag = TRUE;
		      delta_t = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
		      exit(1);
	    }            

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag) cerr << "default delta_t = " << delta_t << "\n";
	if (!a_flag) cerr << "default eta = " << eta << "\n";
	if (!d_flag) cerr << "default dt_out = " << dt_out << "\n";
    }

    if (!D_flag) dt_snap = delta_t;

    b = get_sdyn3(cin);
    
    b->log_history(argc, argv);

    tree3_evolve(b, delta_t, dt_out,
		 dt_snap, snap_cube_size,
		 eta, cpu_time_check);
}

#endif
