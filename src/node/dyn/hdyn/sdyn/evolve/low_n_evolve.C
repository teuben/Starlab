
//// low_n_evolve:  time symmetrized Hermite integrator for low-N systems.
////
//// Options:      -A    specify accuracy parameter [0.02]
////               -c    specify CPU time check interval (s) [3600]
////               -C    specify cube size for snapshot output [10]
////               -d    specify log output interval [none]
////               -D    specify snapshot output interval [none]
////               -e    specify softening parameter [0.0]
////               -f    specify constant or dynamic timesteps [constant]
////               -n    specify number of iterations in symmetrization [0]
////               -q    quiet output [verbose]
////               -s    symmetric timestep [false]
////               -t    specify time span of integration [1]
////               -x    terminate at exact end time [no]
////               -z    specify maximum number of steps to take [unspecified]

#include "sdyn.h"

#ifndef TOOLBOX

void predict_step(sdyn * b,          // n-body system pointer
		  real dt)         // timestep
    {
    if(b->get_oldest_daughter() !=NULL)
        for_all_daughters(sdyn, b, bb)
	    predict_step(bb,dt);
    else
        b->taylor_pred_new_pos_and_vel(dt);
    }

void correct_step(sdyn * b,          // n-body system pointer
		  real new_dt,       // new timestep
		  real prev_new_dt)  // previous new timestep
    {
    if(b->get_oldest_daughter() !=NULL)
        for_all_daughters(sdyn, b, bb)
	    correct_step(bb, new_dt, prev_new_dt);
    else
	{
	b->correct_new_acc_and_jerk(new_dt, prev_new_dt);
        b->correct_new_pos_and_vel(new_dt);
	}
    }

typedef  real  (*tfp)(sdyn *, real);

void step(sdyn * b,        // sdyn array                   
	  real & t,        // time                         
	  real eps,        // softening length             
	  real eta,        // time step parameter
	  real dt,         // time step of the integration 
	  real max_dt,     // maximum time step (till end of the integration)
	  int & end_flag,  // to flag that integration end is reached
	  int & coll_flag,  // to flag that collision has occured
	  tfp the_tfp,     // timestep function pointer
	  int n_iter,      // number of iterations
          int x_flag,      // exact-time termination flag
	  int s_flag)      // symmetric timestep ?
{
    int collision_flag = 0;

    predict_step(b, dt);

    real new_dt = dt;
    for (int i = 0; i <= n_iter; i++) {

	real prev_new_dt = new_dt;
        b->calculate_new_acc_and_jerk_from_new(b, eps*eps, n_iter - i,  // hack
					       collision_flag);
	if (i < n_iter && s_flag) {

	    real end_point_dt = the_tfp(b, eta);

	    new_dt = 0.5 * (dt + end_point_dt);
	    new_dt = dt + 0.5 * (end_point_dt - dt) * (new_dt/prev_new_dt);

	    if (x_flag) {

		if (new_dt >= max_dt) {
		    end_flag = 1;
		    new_dt = max_dt;
		} else
		    end_flag = 0;
	    }
	}
	correct_step(b, new_dt, prev_new_dt);
    }

    for_all_daughters(sdyn, b, bb) bb->store_new_into_old();
    t += new_dt;

    if (collision_flag) end_flag = 1;
    coll_flag = collision_flag;
}

void initialize(sdyn * b,        // sdyn array                   
		real eps)        // softening length             
{
    predict_step(b, 0);

    int collision_flag;
    b->calculate_new_acc_and_jerk_from_new(b, eps*eps, 1, collision_flag);

    for_all_daughters(sdyn, b, bb) bb->store_new_into_old();
}

real  constant_timestep(sdyn * b, real eta)
{
    if (b == NULL) eta = 0; // To keep the HP compiler happy.
    return  eta;
}

real  dynamic_timestep(sdyn * b, real eta)
    {
    real global_min_encounter_time_sq = VERY_LARGE_NUMBER;
    real global_min_free_fall_time_sq = VERY_LARGE_NUMBER;

    for_all_daughters(sdyn, b, bb)
	{
	global_min_encounter_time_sq =
	    min(global_min_encounter_time_sq, bb->get_min_encounter_time_sq());
	global_min_free_fall_time_sq =
	    min(global_min_free_fall_time_sq, bb->get_min_free_fall_time_sq());
	}

    return  eta *
	sqrt(min(global_min_encounter_time_sq, global_min_free_fall_time_sq));
    }

tfp timestep_function_ptr(char * timestep_name)
    {
    if (streq(timestep_name, "constant_timestep"))
	return constant_timestep;
    else if (streq(timestep_name, "dynamic_timestep"))
	return dynamic_timestep;
    else
	{
	cerr << "timestep_function_ptr: no timestep function implemented"
	     << " with name `" << timestep_name << "'" << endl;
	exit(1);
	}
    return NULL; // To keep HP g++ happy!
    }

real calculate_energy(sdyn * b, real & ekin, real & epot)
    {
    ekin = epot = 0;
    for_all_daughters(sdyn, b, bb)
	{
	epot += bb->get_mass() * bb->get_pot();
	ekin += 0.5 * bb->get_mass() * (bb->get_vel() * bb->get_vel());
	}
    epot *= 0.5;

    return ekin + epot;
    }

real calculate_energy_from_scratch(sdyn * b, real & ekin, real & epot)
{
    ekin = epot = 0;

    for_all_daughters(sdyn, b, bi) {
	ekin += bi->get_mass() * bi->get_vel() * bi->get_vel();
	real depot = 0;
	for (sdyn* bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	    depot -= bj->get_mass() / abs(bi->get_pos() - bj->get_pos());
	epot += bi->get_mass() * depot;
    }
    ekin *= 0.5;

    return ekin + epot;
}

void  start_up(sdyn * b, real & n_steps)
    {
    if(b->get_oldest_daughter() !=NULL)
	{
	n_steps = b->get_n_steps();
	if (n_steps == 0)
	    b->prepare_root();

	for_all_daughters(sdyn, b, bb)
	    start_up(bb, n_steps);
	}
    else
	{
	if (n_steps == 0)
	    b->prepare_branch();
	}
    }

void  clean_up(sdyn * b, real n)
    {
    b->set_n_steps(n);
    }

// system_in_cube: return TRUE if most of the (top-level) nodes in the
//                 system are contained within the specified cube.

#define MOST 0.75

local bool system_in_cube(sdyn* b, real cube_size) {

    int n = 0;
    for_all_daughters(sdyn, b, bi) n++;

    int n_in = 0;
    for_all_daughters(sdyn, b, bb) {
	int in = 1;
	for (int k = 0; k < 3; k++)
	    if (abs(bb->get_pos()[k]) > cube_size) in = 0;
	n_in += in;
    }

    if (n_in >= MOST * n)
	return TRUE;
    else
	return FALSE;
}

local void copy_node_partial(sdyn*b, sdyn* c)
{
    c->set_index(b->get_index());
    c->set_mass(b->get_mass());
    c->set_time(b->get_time());
    c->set_pos(b->get_pos());
    c->set_vel(b->get_vel());
}

local sdyn* copy_flat_tree(sdyn* b)
// Make a partial copy of the flat tree with b as root, and return a
// pointer to it.  Intention is for use with xstarplot, so only 
// copy pointers, index, mass, time, pos, and vel.
{
    sdyn* root = new sdyn;

    root->set_parent(NULL);
    root->set_elder_sister(NULL);
    root->set_younger_sister(NULL);
    copy_node_partial(b, root);

    sdyn* s = NULL;
    for_all_daughters(sdyn, b, bi) {
	sdyn* ci = new sdyn;
	ci->set_parent(root);
	ci->set_oldest_daughter(NULL);
	ci->set_younger_sister(NULL);
	ci->set_elder_sister(s);
	if (s)
	    s->set_younger_sister(ci);
	else
	    root->set_oldest_daughter(ci);
	copy_node_partial(bi, ci);
	s = ci;
    }

    return root;
}

local void delete_node(sdyn* b)
// Recursively delete node b and its descendents.
{
    sdyn* bi = b->get_oldest_daughter();
    while (bi) {
	sdyn* tmp = bi->get_younger_sister();
	delete_node(bi);
	bi = tmp;
    }
    delete b;
}

// pp: Recursively pretty-print a node:

void pp(sdyn* b, ostream & s, int level = 0) {

    s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    b->pretty_print_node(s);
    s << " \t"<< b->get_mass() << " \t"
      << b->get_pos() << "   "
      << b->get_vel() <<endl;

    for_all_daughters(sdyn, b, daughter) pp(daughter, s, level + 1);	
}


#define N_STEP_CHECK 1000 // Interval for checking CPU time

bool low_n_evolve(sdyn * b,       // sdyn array
		  real delta_t,   // time span of the integration
		  real dt_out,    // output time interval
		  real dt_snap,   // snapshot output interval
		  real snap_cube_size,
		  real eps,       // softening length 
		  real eta,       // time step parameter
		  int x_flag,     // exact-time termination flag
		  char * timestep_name,
		  int s_flag,     // symmetric timestep ?
		  int n_iter,     // number of iterations
		  real n_max,     // if > 0: max. number of integration steps
		  real cpu_time_check,
		  real dt_print,  // external print interval
		  sdyn_print_fp   // pointer to external print function
		       print)
    {

    bool terminate = false;

    real t = b->get_time();
    real tr = t + (real)b->get_time_offset();  // total real time in N-body units

    real t_end = t + delta_t;      // final time, at the end of the integration
    real t_out = t + dt_out;       // time of next diagnostic output
    real t_snap = t + dt_snap;     // time of next snapshot;
    real t_print = t + dt_print;   // time of next printout;

    real n_steps;
    int count_steps = 0;
    real cpu_init = cpu_time();
    real cpu_save = cpu_init;

    tfp  the_tfp = timestep_function_ptr(timestep_name);

    start_up(b, n_steps);
    initialize(b, eps);

    cerr.precision(16);

    real ekin, epot;
    calculate_energy(b, ekin, epot);

    if (t_out <= t_end && n_steps == 0)
	{
	  cerr << "Time = " << t << " " << tr << "  n_steps = " << n_steps
	    << "  Etot = " << ekin + epot << endl;
	}

    if (b->get_n_steps() == 0)             // should be better interfaced
	{                                  // with start_up ; to be done
	b->set_e_tot_init(ekin + epot);
	b->clear_de_tot_abs_max();
	}

    int end_flag = 0;
    int coll_flag = 0;
    while (t < t_end && !end_flag)
	{
        real max_dt = t_end - t;
        real dt = the_tfp(b, eta);

        end_flag = 0;
        if (dt > max_dt && x_flag)
	    {
	    end_flag = 1;
	    dt = max_dt;
	    }

        step(b, t, eps, eta, dt, max_dt, end_flag, coll_flag, 
	     the_tfp, n_iter,
	     x_flag, s_flag);
	b->set_time(t);                  // should be prettified some time soon
	for_all_daughters(sdyn, b, bi)
	    bi->set_time(t);

	n_steps += 1;
	count_steps++;

	if(coll_flag>=1) {
	  //cerr << "Collision has occured: " << coll_flag << endl;
	  terminate = false;
	  return terminate;
	}

//      Check for (trivial) output to cerr...

	if (t >= t_out && n_steps==0)
	    {
	    calculate_energy(b, ekin, epot);
	    cerr << "Time = " << t << " " << tr << "  n_steps = " << n_steps
		 << "  Etot = " << ekin + epot << endl;
	    t_out += dt_out;
	    }

//      ...and (not-so-trivial) output handled elsewhere.

	if (t >= t_print && print != NULL) 
	    {
	    (*print)(b);
	    t_print += dt_print;
	    }

//      Output a snapshot to cout at the scheduled time.

        if(n_max > 0 && n_steps >= n_max)
	    end_flag = 1;

	if (t >= t_snap && system_in_cube(b, snap_cube_size)) {

	    // Looks like we can't just make a tree and flatten it
	    // again before continuing--better to work with a copy.

/*	    sdyn* c = copy_flat_tree(b);

	    bool dynamics = FALSE;
	    bool stability = FALSE;
	    int k_max = 2;
	    make_tree(c, dynamics, stability, k_max, 0);
	    put_node(cout, *c);
	    delete_node(c);
*/
	    put_node(cout, *b);

	    cout << flush; 
	    t_snap += dt_snap;
	}

	// Check the number of steps and the CPU time every N_STEP_CHECK steps.
	// Note that the printed CPU time is the time since this routine was
	// entered.

	if (count_steps >= N_STEP_CHECK) {
	    count_steps = 0;

	    // added (SPZ: 12 Oct 2000)
	    if (b->get_n_steps() > MAX_N_STEPS) {
	      //return;
	      cerr << "Terminate after " << b->get_n_steps() 
		   << " steps" << endl;
	      terminate = true;
	      }

	    if (cpu_time() - cpu_save > abs(cpu_time_check)) {
		cpu_save = cpu_time();
		calculate_energy(b, ekin, epot);
		cerr.precision(6);
		cerr << "low_n3_evolve:  CPU time = " << cpu_save - cpu_init;
		cerr.precision(13);
		cerr << "  time = " << t;
		cerr.precision(6);
		cerr << "  offset = " << b->get_time_offset() << endl;
		cerr << "                n_steps = " << n_steps
		     << "  Etot = " << ekin + epot
		     << "  dt = " << dt
		     << endl << flush;

		if (cpu_time_check < 0) return terminate;

		}
	    }
	}
 
   clean_up(b, n_steps);       // too late for snapshot?

    return terminate;

    }

#else

main(int argc, char **argv)
{
    sdyn* b;             // pointer to the nbody system
    int   n_iter = 1;	 // number of iterations (0: explicit; >=1: implicit)
    real  n_max = -1;    // if > 0: maximum number of integration steps
    
    real  delta_t = 10;   // time span of the integration
    real  eta = 0.02;    // time step parameter (for fixed time step,
                         //   equal to the time step size; for variable
                         //   time step, a multiplication factor)
    real  dt_out = VERY_LARGE_NUMBER;
                         // output time interval
    real  dt_snap = VERY_LARGE_NUMBER;
                         // snap output interval
    real snap_cube_size = 10;

    real  cpu_time_check = 3600;

    real  eps = 0.0;     // softening length 	       	   
    char  *timestep_name = "dynamic_timestep";

    bool  a_flag = FALSE;
    bool  d_flag = FALSE;
    bool  D_flag = FALSE;
    bool  e_flag = FALSE;
    bool  f_flag = FALSE;
    bool  n_flag = FALSE;
    bool  q_flag = FALSE;
    bool  r_flag = FALSE;
    bool  s_flag = TRUE;    // symmetric timestep ?
    bool  t_flag = FALSE;
    bool  x_flag = TRUE;    // if true: termination at the exact time of
                            //          of the final output, by
                            //          adjustment of the last time step;
                            // if false: no adjustment of the last time step,
                            //           as a consequence the time of final
                            //           output might be slightly later than
                            //           the time specified.
    bool  z_flag = FALSE;     // to specify termination after n_max steps

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:e:f:n:qr:st:xz:";

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
	    case 'e': e_flag = TRUE;
		      eps = atof(poptarg);
		      break;
	    case 'f': f_flag = TRUE;
		      timestep_name = poptarg;
		      break;
	    case 'n': n_flag = TRUE;
		      n_iter = atoi(poptarg);
		      break;
	    case 'q': q_flag = TRUE;
		      break;
	    case 's': s_flag = FALSE;
		      break;
	    case 't': t_flag = TRUE;
		      delta_t = atof(poptarg);
		      break;
	    case 'x': x_flag = FALSE;
		      break;
	    case 'z': z_flag = TRUE;
		      n_max = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
	    }            

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag) cerr << "default delta_t = " << delta_t << "\n";
	if (!a_flag) cerr << "default eta = " << eta << "\n";
	if (!d_flag) cerr << "default dt_out = " << dt_out << "\n";
	if (!e_flag) cerr << "default eps = " << eps << "\n";
	if (!f_flag) cerr << "default timestep_name = " << timestep_name
	                  << "\n";
	if (!n_flag) cerr << "default n_iter = " << n_iter << "\n";
	if (!s_flag) cerr << "s_flag (symmetric timestep ?) = FALSE" << "\n";
	if (x_flag) cerr << "default termination: at exact t_end" << "\n";
	if (!z_flag) cerr << "default n_max = " << n_max << "\n";
    }

    if (!D_flag) dt_snap = delta_t; // Guarantee output at end

    b = get_sdyn(cin);
    b->log_history(argc, argv);

    cpu_init();
    low_n_evolve(b, delta_t, dt_out, dt_snap, snap_cube_size,
		 eps, eta, x_flag,
		 timestep_name, s_flag, n_iter, n_max,
		 cpu_time_check);

}

#endif
