
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// hermite.C
// hermite scheme with simple (flat-tree) hdyn structure
//
// J. Makino    92-12-04:  Seems to work... Energy conservation looks good.
// S. McMillan  93-04-02:  Cleaned up and modified command-line interface.
// J. Makino    96-07-xx:  Seems not to work... Energy conservation looks bad.

void hdyn::flat_set_first_timestep(real eta_for_firststep, real max_step_size)
{
    real a1scaled2 = jerk * jerk;
    real a2 = acc * acc;
    real newstep;

    if (a1scaled2 > 0.0) {
	newstep = eta_for_firststep * sqrt(a2 / a1scaled2);
    } else {
	newstep = eta_for_firststep * 0.125;
    }
    newstep = adjust_number_to_power(newstep, max_step_size);
    timestep = newstep;
}

real flat_new_timestep(vector & at3,	// third order term
		       vector & bt2,	// 2nd order term
		       vector & jerk,	// 1st order term
		       vector & acc,	// 0th order term
		       real timestep,	// old timestep
		       real time,	// present time
		       real eta,	// accuracy parameter
		       real dtmax)	// maximum stepsize
{
    real a3scaled2 = 36 * at3 * at3;
    real a2scaled2 = 4 * bt2 * bt2;
    real a1scaled2 = timestep * jerk * jerk;
    real a2 = acc * acc;

// simple criterion:
    real newstep = eta * sqrt(sqrt(a2 / a2scaled2)) * timestep;

    if (newstep < timestep) {
	return timestep * 0.5;
    } else if (newstep < 2 * timestep) {
	return timestep;
    } else if (fmod(time, 2 * timestep) != 0.0 || 2 * timestep > dtmax) {
	return timestep;
    } else {
	return timestep * 2;
    }
}

void hdyn::flat_update(const real eta, const real dtmax)
{
    vector at3 = 2 * (old_acc - acc) + timestep * (old_jerk + jerk);
    vector bt2 = -3 * (old_acc - acc) - timestep * (2 * old_jerk + jerk);
    time = time + timestep;
    timestep = flat_new_timestep(at3, bt2, jerk, acc,
				 timestep, time, eta, dtmax);
}

void accumulate_energies(hdyn * b, real & epot, real & ekin, real & etot)
{
    if (b->get_oldest_daughter() != NULL) {
	for (hdyn * bb = b->get_oldest_daughter(); bb != NULL;
	     bb = bb->get_younger_sister()) {
	    accumulate_energies(bb, epot, ekin, etot);
	}
    } else {
	real mi = b->get_mass();
	epot += 0.5 * mi * b->get_pot();
	vector vel = b->get_vel();
	ekin += 0.5 * mi * vel * vel;
    }
    etot = ekin + epot;
}

void compute_energies(hdyn * b, real & epot, real & ekin, real & etot)
{
    epot = ekin = etot = 0;
    accumulate_energies(b, epot, ekin, etot);
}

void print_energies(hdyn * b, int mode)
{
    real epot;
    real ekin;
    real etot;
    static real old_etot;
    compute_energies(b, epot, ekin, etot);
    cerr << "    Energies:  " << ekin << "  " << epot << "  " << etot;
    if (mode) {
	real de = (etot - old_etot) / etot;
	cerr << ";  de = " << de << "\n";
    } else {
	old_etot = etot;
	cerr << "\n";
    }
    hdyn *dummy = b;		// to make the compiler happy

}

hdyn *particle_to_move(hdyn * b, real & tmin)
{
    hdyn *bmin = NULL;
    if (b->get_oldest_daughter() != NULL) {
	for (hdyn * bb = b->get_oldest_daughter(); bb != NULL;
	     bb = bb->get_younger_sister()) {
	    hdyn *btmp;
	    real ttmp = 1e100;
	    btmp = particle_to_move(bb, ttmp);
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

void flat_initialize_system_phase2(hdyn * b, hdyn * root, real eps,
				   real eta, real dtout)
{
    if (b->get_oldest_daughter() != NULL) {
	for (hdyn * bb = b->get_oldest_daughter(); bb != NULL;
	     bb = bb->get_younger_sister()) {
	    flat_initialize_system_phase2(bb, root, eps, eta, dtout);
	}
    } else {
	b->clear_interaction();
	b->flat_calculate_acc_and_jerk(root, eps * eps);
	b->store_old_force();
	b->flat_set_first_timestep(0.5 * eta, dtout);
    }
}

local void shift_cm(hdyn * b, real snap_cube_size)
{
    real mass = 0;
    vector cmpos = vector(0, 0, 0);

    for_all_daughters(hdyn, b, bb) {
	vector x = bb->get_pos();
	if (abs(x[0]) < snap_cube_size
	    && abs(x[1]) < snap_cube_size
	    && abs(x[2]) < snap_cube_size) {
	    mass += bb->get_mass();
	    cmpos += bb->get_mass() * bb->get_pos();
	}
    }

    if (mass > 0) {
	cmpos /= mass;
	for_all_daughters(hdyn, b, bbb) bbb->inc_pos(-cmpos);
    }
}

void evolve_system(hdyn * b,		// hdyn array                   
		   real delta_t,	// time span of the integration 
		   real eta,		// accuracy parameter 
		   real dtout,		// output time interval
		   real dt_snap,	// snapshot output interval
		   real eps,		// softening length             
		   real snap_cube_size)
{
    real t = b->get_time();	// current time
    real t_end = t + delta_t;	// final time, at the end of the integration
    real t_next = t + dtout;	// time of next output
    real t_snap = t + dt_snap;	// time of next snapshot;

    int steps = 0;

    initialize_system_phase1(b, t);
    flat_initialize_system_phase2(b, b, eps, eta, dtout);

    // Initial output (to cerr, note) only if later output is anticipated.

    if (t_next <= t_end) {
	cerr << "Time = " << t << ",  steps = " << steps << "\n";
	print_energies(b, 0);
    }

    while (t <= t_end) {
	hdyn *bi;
	real ttmp = 1e100;
	bi = particle_to_move(b, ttmp);

	// Check for output and termination before taking the step.

	if (ttmp > t_next) {
	    cerr << "Time = " << t << ",  steps = " << steps << "\n";
	    print_energies(b, 1);
	    t_next += dtout;
	}

	// Output a snapshot to cout at the scheduled time, or at end of run.
	// Use snap_cube_size to force all particles in the cube to have C.M.
	// at rest at the origin.

	if (ttmp > t_snap || ttmp > t_end) {
	    if (snap_cube_size < VERY_LARGE_NUMBER)
		shift_cm(b, snap_cube_size);
	    put_node(cout, *b);
	    cout << flush;
	    t_snap += dt_snap;
	}
	if (ttmp > t_end)
	    break;

	t = ttmp;
	predict_loworder_all(b, t);
	bi->clear_interaction();
	bi->flat_calculate_acc_and_jerk(b, eps * eps);
	bi->correct();
	bi->flat_update(eta, dtout);
	bi->store_old_force();
	steps++;
    }
}

void print_usage_and_exit()
{
    cerr <<
	 "usage: hermite -t # -a # -d # -D # -e # "
	 << "[-c \"...\"]      "
	 << "for t (time span),\n a (accuracy parameter ), "
	 << "d (output interval), D (snapshot output interval), "
	 << "e (softening length),\n";
    exit(1);
}

main(int argc, char **argv)
{
    hdyn *b;			// hdyn root node

    real delta_t = 10;		// time span of the integration
    real dtout = .25;		// output interval--make a power of 0.5
    real dt_snap;		// snap output interval
    real eps = 0.05;		// softening length               
    real eta = 0.05;		// time step parameter

    char *comment;		// comment string

    real snap_cube_size = VERY_LARGE_NUMBER;

    extern char *poptarg;
    int pgetopt(int, char **, char *), c;

    bool a_flag = FALSE;
    bool c_flag = FALSE;
    bool d_flag = FALSE;
    bool D_flag = FALSE;
    bool e_flag = FALSE;
    bool q_flag = FALSE;
    bool t_flag = FALSE;

    while ((c = pgetopt(argc, argv, "a:c:C:d:D:e:qt:")) != -1) {
	switch (c) {
	    case 'a':
		a_flag = TRUE;
		eta = atof(poptarg);
		break;
	    case 'c':
		c_flag = TRUE;
		comment = poptarg;
		break;
	    case 'C':
		snap_cube_size = atof(poptarg);
		break;
	    case 'd':
		d_flag = TRUE;
		dtout = atof(poptarg);
		break;
	    case 'D':
		D_flag = TRUE;
		dt_snap = atof(poptarg);
		break;
	    case 'e':
		e_flag = TRUE;
		eps = atof(poptarg);
		break;
	    case 'q':
		q_flag = TRUE;
		break;
	    case 't':
		t_flag = TRUE;
		delta_t = atof(poptarg);
		break;
	    case '?':
		print_usage_and_exit();
	}
    }

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag)
	    cerr << "default delta_t = " << delta_t << "\n";
	if (!d_flag)
	    cerr << "default dtout = " << dtout << "\n";
	if (!a_flag)
	    cerr << "default eta = " << eta << "\n";
	if (!e_flag)
	    cerr << "default eps = " << eps << "\n";
    }
    if (!D_flag)
	dt_snap = delta_t;

    b = get_hdyn(cin);

    if (c_flag == TRUE)
	b->log_comment(comment);
    b->log_history(argc, argv);

    evolve_system(b, delta_t, eta, dtout, dt_snap, eps, snap_cube_size);
}
