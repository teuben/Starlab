
//// leapfrog:  Leapfrog integrator for flat N-body dyn systems.
////
//// Options:   -a    specify timestep [0.02]
////            -c    add a comment to the output snapshot [false]
////            -d    specify log output time interval [0.25]
////            -D    specify snapshot output time interval [none]
////            -e    specify softening parameter [0.05]
////            -q    echo parameters to cerr [true]
////            -t    specify time span of integration [1]
////            -x    terminate precisely at specified end time [no]

#include "dyn.h"

#ifdef TOOLBOX

local void predict_step(dyn * b,          // n-body system pointer
			real dt)          // timestep
{
    if (b->get_oldest_daughter() !=NULL) {
	for (dyn * bb = b->get_oldest_daughter(); bb != NULL;
	    bb = bb->get_younger_sister()) {
	    predict_step(bb, dt);
	}
    } else {
	b->inc_vel( 0.5 * dt * b->get_acc() );
	b->inc_pos( dt * b->get_vel() );
    }
}

local void correct_step(dyn * b,          // n-body system pointer
			real dt)          // timestep
{
    if (b->get_oldest_daughter() !=NULL) {
	for (dyn * bb = b->get_oldest_daughter(); bb != NULL;
	    bb = bb->get_younger_sister()) {
	    correct_step(bb,dt);
	}
    } else {
	b->inc_vel( 0.5 * dt * b->get_acc() );
    }
}

local void step(real& t,        // time                         
		dyn* b,         // dyn array                   
		real dt,        // time step of the integration 
		real eps)       // softening length             
{
    t += dt;
    
    predict_step(b, dt);
    b->calculate_acceleration(b, eps*eps);
    correct_step(b, dt);
}

local void evolve(real& t,        // time                         
		  dyn* b,         // dyn array                   
		  real delta_t,   // time span of the integration 
		  real dt,        // (fixed) integration time step
		  real dt_out,    // output time interval
		  real dt_snap,   // snapshot output interval
		  real eps,       // softening length             
		  int x_flag)     // exact-time termination flag
{
    real t_end = t + delta_t;      // final time, at the end of the integration
    real t_out = t + dt_out;       // time of next diagnostic output
    real t_snap = t + dt_snap;     // time of next snapshot;
    int steps = 0;

    
    b->calculate_acceleration(b, eps*eps);

    print_recalculated_energies(b, 0, eps*eps, 0);

    while (t < t_end) {

        int termination_flag;

        if (t + dt < t_end)
	    termination_flag = 0;
	else
	    {
	    termination_flag = 1;
	    if (x_flag)
		dt = t_end - t;
	    }

        step(t, b, dt, eps);
	steps++;

	// Check for (trivial) output to cerr.

	if (t >= t_out) {
	    cerr << "Time = " << t << "  steps = " << steps << endl;
	    print_recalculated_energies(b, 1, eps*eps, 0);
	    t_out += dt_out;
	}

	// Output a snapshot to cout at the scheduled time, or at end of run.

	if (termination_flag == 1) {
	    put_node(b);
	    cout << flush;
	    break;                          // end the run
	}

	if (t >= t_snap) {
	    put_node(b);             	    // do not synchronize all particles
	    cout << flush;
	    t_snap += dt_snap;              // and continue the run
	}
    }
}

main(int argc, char **argv)
{
    dyn* b;              // pointer to the nbody system
    real  t = 0;         // time

    real  delta_t = 1;   // time span of the integration
    real  dt_out = .25;  // output time interval
    real  dt = 0.02;     // time step of the integration
    real  dt_snap;       // snap output interval
    real  eps = 0.05;    // softening length 	       	   
    char  *comment;

    bool  a_flag = FALSE;
    bool  c_flag = FALSE;
    bool  d_flag = FALSE;
    bool  D_flag = FALSE;
    bool  e_flag = FALSE;
    bool  q_flag = FALSE;
    bool  t_flag = FALSE;
    bool  x_flag = FALSE;   // if true: termination at the exact time of
                            //          of the final output, by
                            //          adjustment of the last time step;
                            // if false: no adjustment of the last time step,
                            //           as a consequence the time of final
                            //           output might be slightly later than
                            //           the time specified.

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:c:d:D:e:qt:x";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'a': a_flag = TRUE;
		      dt = atof(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'd': d_flag = TRUE;
		      dt_out = atof(poptarg);
		      break;
	    case 'D': D_flag = TRUE;
		      dt_snap = atof(poptarg);
		      break;
	    case 'e': e_flag = TRUE;
		      eps = atof(poptarg);
		      break;
	    case 'q': q_flag = TRUE;
		      break;
	    case 't': t_flag = TRUE;
		      delta_t = atof(poptarg);
		      break;
	    case 'x': x_flag = TRUE;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag) cerr << "default delta_t = " << delta_t << "\n";
	if (!a_flag) cerr << "default dt = " << dt << "\n";
	if (!d_flag) cerr << "default dt_out = " << dt_out << "\n";
	if (!e_flag) cerr << "default eps = " << eps << "\n";
	if (!x_flag) cerr << "default termination: not at exact t_end" << "\n";
    }

    if (!D_flag) dt_snap = delta_t;

    b = get_dyn();
    
    if (c_flag == TRUE) b->log_comment(comment);
    b->log_history(argc, argv);

    evolve(t, b, delta_t, dt, dt_out, dt_snap, eps, x_flag);
    rmtree(b);
}

#endif
