
//// makesphere: construct a simple homogeneous sphere, with
////
////              M = 1, T/U = -1/2, E = -1/4.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -i    number the particles sequentially [don't number]
////              -n    specify number of particles [no default]
////              -o    echo value of random seed [don't echo]
////              -s    specify random seed [random from system clock]
////              -u    leave unscaled [scale to E=-1/4, M = 1, R = 1]
////
////  If the "-u" flag is set, the particles are left unscaled,
////  uniformly distributed in a sphere with unit radius, with
////  all velocity components uniformly distributed in [-1, 1]
////  and all masses equal to 1/n.

#include "dyn.h"

#ifndef TOOLBOX

// Don't make local -- used elsewhere.

void  makesphere(dyn * root, int n,
	       int u_flag)			// default = false
{
    real radius, costheta, sintheta, phi;
    dyn  * bi;

    root->set_mass(1);
    real pmass = 1.0 / n;

    for (bi = root->get_oldest_daughter(); bi != NULL;
         bi = bi->get_younger_sister()) {
	bi->set_mass(pmass);

	radius = pow(randinter(0, 1), 1.0/3.0);
	costheta = randinter(-1.0, 1.0);
	sintheta = 1 - costheta*costheta;
	if (sintheta > 0) sintheta = sqrt(sintheta);
	phi = randinter(0.0, TWO_PI);
        bi->set_pos(vec(radius * sintheta * cos(phi),
			   radius * sintheta * sin(phi),
			   radius * costheta));

        bi->set_vel(vec(randinter(-1,1),
			   randinter(-1,1),
			   randinter(-1,1)));
    }

//  Transform to center-of-mass coordinates and optionally
//  scale to standard parameters.

    root->to_com();
    putrq(root->get_log_story(), "initial_mass", 1.0);

    if (!u_flag && n > 1) {

        real kinetic, potential;

	// Note: scale_* operates on internal energies.

	get_top_level_energies(root, 0.0, potential, kinetic);
	scale_virial(root, -0.5, potential, kinetic);	// scales kinetic
	real energy = kinetic + potential;
	scale_energy(root, -0.25, energy);		// scales energy
	putrq(root->get_log_story(), "initial_total_energy", -0.25);
	putrq(root->get_log_story(), "initial_rvirial", 1.0);
	putrq(root->get_dyn_story(), "total_energy", -0.25);
    }
}

#else

#define  SEED_STRING_LENGTH  60

main(int argc, char ** argv) {
    int  i;
    int  n;
    int  input_seed, actual_seed;
    bool c_flag = false;
    bool C_flag = false;
    int  i_flag = FALSE;
    int  n_flag = FALSE;
    int  o_flag = FALSE;
    int  s_flag = FALSE;
    int  u_flag = FALSE;

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:Cin:os:u";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'n': n_flag = true;
		      n = atoi(poptarg);
		      break;
	    case 'o': o_flag = true;
                      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'u': u_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            
    
    if (!n_flag) {
        cerr << "makesphere: must specify the number # of";
	cerr << " particles with -n#\n";
	exit(1);
    }
    
    if (n < 1) {
        cerr << "makesphere: n < 1 not allowed\n";
	exit(1);
    }

    dyn *b, *by, *bo;
    b = new dyn();
    bo = new dyn();
    if (i_flag) bo->set_label(1);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    for (i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag) by->set_label(i+1);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    if (C_flag) b->set_col_output(true);
    if (c_flag) b->log_comment(comment);

    b->log_history(argc, argv);

    if (!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makesphere: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    if (n > 0) makesphere(b, n, u_flag);

    put_dyn(b);
    rmtree(b);
}

#endif

