
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// flat_hermite.C: hermite scheme with simple (flat-tree) _dyn_ structure
//
// J. Makino    92-12-04:  Seems to work... Energy conservation looks good.
// S. McMillan  93-04-02:  Cleaned up and modified command-line interface.
// J. Makino    96-07-xx:  Seems not to work... Energy conservation looks bad.
// S. McMillan  99-06-xx:  Seems to work... Energy conservation looks good.

#include "_dyn_.h"

local void accumulate_energies(_dyn_ * b,
			       real & epot, real & ekin, real & etot)
{
    // Recursion is one level deep for flat trees...

    if (b->get_oldest_daughter()) {

	for_all_daughters(_dyn_, b, bb)
	    accumulate_energies(bb, epot, ekin, etot);	// recursion...

    } else {

	real mi = b->get_mass();
	epot += 0.5 * mi * b->get_pot();
	vector vel = b->get_vel();
	ekin += 0.5 * mi * vel * vel;

    }
    etot = ekin + epot;
}

local void compute_energies(_dyn_ * b,
			    real & epot, real & ekin, real & etot)
{
    epot = ekin = etot = 0;
    accumulate_energies(b, epot, ekin, etot);
}

local void print_energies(_dyn_ * b, int mode)
{
    real epot, ekin, etot;
    static real old_etot = 0;

    compute_energies(b, epot, ekin, etot);
    cerr << "    Energies:  " << ekin << "  " << epot << "  " << etot;

    if (mode) {
	real de = (etot - old_etot) / etot;
	cerr << ";  de = " << de << endl;
    } else {
	old_etot = etot;
	cerr << endl;
    }
}

local _dyn_ *particle_to_move(_dyn_ * b, real & tmin)
{
    // Scheduler returns a pointer to a single particle.  Since steps
    // are forced to be powers of 2 (see _dyn_ev.C), it would be better
    // if we returned a list of nodes to move, as in hdyn/evolve/kira.

    _dyn_ *bmin = NULL;

    if (b->get_oldest_daughter()) {

	for_all_daughters(_dyn_, b, bb) {
	    real ttmp = VERY_LARGE_NUMBER;
	    _dyn_ *btmp = particle_to_move(bb, ttmp);	// recursion...
	    if (ttmp < tmin) {
		tmin = ttmp;
		bmin = btmp;
	    }
	}

    } else {

	if (tmin > b->get_next_time()) {
	    bmin = b;
	    tmin = b->get_next_time();
	}

    }

    return bmin;
}

local void flat_initialize_system_phase1(_dyn_ * b,  real t)
{
    if (b->get_oldest_daughter())
	for_all_daughters(_dyn_, b, d)
	    flat_initialize_system_phase1(d, t);	// recursion...

    b->set_time(t);
    b->set_system_time(t);
    b->init_pred();
}

local void flat_initialize_system_phase2(_dyn_ * b, _dyn_ * root, real eps,
					 real eta, real dtout)
{
    if (b->get_oldest_daughter() != NULL) {

	for_all_daughters(_dyn_, b, bb)
	    flat_initialize_system_phase2(bb, root,
					  eps, eta, dtout); // recursion...

    } else {

	b->clear_interaction();
	b->flat_calculate_acc_and_jerk(root, eps * eps);
	b->store_old_force();
	b->flat_set_first_timestep(0.5 * eta, dtout);

    }
}

local void evolve_system(_dyn_ * b,	// root node
			 real delta_t,	// time span of the integration 
			 real eta,	// accuracy parameter 
			 real dtout,	// output time interval
			 real dt_snap,	// snapshot output interval
			 real eps,	// softening length             
			 real snap_cube_size)
{
    real t = b->get_time();	// current time
    real t_end = t + delta_t;	// final time, at the end of the integration
    real t_next = t + dtout;	// time of next output
    real t_snap = t + dt_snap;	// time of next snapshot;

    int steps = 0;

    flat_initialize_system_phase1(b, t);
    flat_initialize_system_phase2(b, b, eps, eta, dtout);

    // Initial output (to cerr, note) only if later output is anticipated.

    if (t_next <= t_end) {
	cerr << "Time = " << t << ",  steps = " << steps << endl;
	print_energies(b, 0);
    }

    // Main integration loop.

    while (t <= t_end) {

	real ttmp = 1e100;
	_dyn_ *bi = particle_to_move(b, ttmp);

	// Check for output and termination before taking the step.

	if (ttmp > t_next) {
	    cerr << "Time = " << t << ",  steps = " << steps << endl;
	    print_energies(b, 1);
	    t_next += dtout;
	}

	// Output a snapshot to cout at the scheduled time, or at end of run.
	// Use snap_cube_size to force all particles in the cube to have C.M.
	// at rest at the origin.

	if (ttmp > t_snap || ttmp > t_end) {
	    put_node(b);
	    cout << flush;
	    t_snap += dt_snap;
	}
	if (ttmp > t_end) break;

	b->set_system_time(t=ttmp);
	predict_loworder_all(b, t);

	bi->clear_interaction();

	bi->flat_calculate_acc_and_jerk(b, eps * eps);

	bi->flat_correct();
	bi->flat_update(eta, dtout);
	bi->store_old_force();

	steps++;
    }
}

main(int argc, char **argv)
{
    _dyn_ *b;			// root node

    real delta_t = 10;		// time span of the integration
    real dtout = .25;		// output interval--make a power of 0.5
    real dt_snap;		// snap output interval
    real eps = 0.05;		// softening length               
    real eta = 0.05;		// time step parameter

    char *comment;		// comment string

    real snap_cube_size = VERY_LARGE_NUMBER;

    bool a_flag = FALSE;
    bool c_flag = FALSE;
    bool d_flag = FALSE;
    bool D_flag = FALSE;
    bool e_flag = FALSE;
    bool q_flag = FALSE;
    bool t_flag = FALSE;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:c:C:d:D:e:qt:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)

	switch (c) {
	    case 'a':	a_flag = TRUE;
	    		eta = atof(poptarg);
			break;
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'C':	snap_cube_size = atof(poptarg);
			break;
	    case 'd':	d_flag = TRUE;
	    		dtout = atof(poptarg);
			break;
	    case 'D':	D_flag = TRUE;
			dt_snap = atof(poptarg);
			break;
	    case 'e':	e_flag = TRUE;
			eps = atof(poptarg);
			break;
	    case 'q':	q_flag = TRUE;
			break;
	    case 't':	t_flag = TRUE;
			delta_t = atof(poptarg);
			break;
            case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			exit(1);
	}

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag)
	    cerr << "default delta_t = " << delta_t << endl;
	if (!d_flag)
	    cerr << "default dtout = " << dtout << endl;
	if (!a_flag)
	    cerr << "default eta = " << eta << endl;
	if (!e_flag)
	    cerr << "default eps = " << eps << endl;
    }
    if (!D_flag)
	dt_snap = delta_t;

    b = get__dyn_();

    if (c_flag == TRUE)
	b->log_comment(comment);
    b->log_history(argc, argv);

    evolve_system(b, delta_t, eta, dtout, dt_snap, eps, snap_cube_size);
}
