
//// makepowerlaw: construct a system of test particles drawn from a
////               power-law mass distribution, in virial equilibrium
////               with the field due to that distribution.  The mass
////               distribution is of the form:
////
////                     M(r)  =  A r^x		(0 < x < 3)
////
////               except that r is modified to include a softening scale
////               parameter:
////
////                     r  -->  sqrt(r^2 + a^2)
////
//// We assume that x > 0.  Particle parameters will be chosen so that the
//// total mass of the N-body system is 1, independent of the actual mass
//// of the background distribution.  For now, we are interested only in
//// test particles.
////
//// Options:     -A    specify the coefficient A [1]
////              -a/R  specify scale [1]
////              -c    add a comment to the output snapshot [false]
////              -i    number the particles sequentially [don't number]
////              -n    specify number of particles [no default]
////              -o    echo value of random seed [don't echo]
////              -r    maximum radius of the N-body system [1]
////              -s    specify random seed [random from system clock]
////              -x    specify exponent [1]
////
//// The output snapshot will have the external field already enabled, and
//// will contain a flag to ignore internal interactions.

#include "dyn.h"

#ifdef TOOLBOX

local void  makepowerlaw(dyn * root, int n,
			 real A, real a, real x, real r_max)
{
    real xi = 1/x, pmass = 1./n, armax = pow(a/r_max, x),
	 vc = sqrt(3*A*pow(a, x-1)/(3-x));
    root->set_mass(1);

    for_all_daughters(dyn, root, bi) {

	bi->set_mass(pmass);

	// Radii are distributed between 0 and r_max, roughly uniformly
	// in r^x.  Power-law is OK if a << r_max; otherwise, better to
	// use a "core-halo" approximation to m(r):
	//
	//	m(r)  ~  (r/a)^3 (a/r_max)^x	(r < a)
	//		 (r/r_max)^x		(r > a)

	real radius;
	if (randinter(0, 1) < armax)
	    radius = a*pow(randinter(0, 1), 1./3);
	else
	    radius = r_max*pow(randinter(armax, 1), xi);

	real costheta = randinter(-1.0, 1.0);
	real sintheta = 1 - costheta*costheta;
	if (sintheta > 0) sintheta = sqrt(sintheta);
	real phi = randinter(0.0, TWO_PI);

        bi->set_pos(radius*vector(sintheta * cos(phi),
				  sintheta * sin(phi),
				  costheta));

	// Choose velocities to ensure equilibrium in the external field.

	real vrms = vc;
	if (radius > a) vrms *= pow(radius/a, (x-1)/2);

        bi->set_vel(vrms*vector(randinter(-1,1),
				randinter(-1,1),
				randinter(-1,1)));
    }

    putrq(root->get_log_story(), "initial_mass", 1.0);
}

#define  SEED_STRING_LENGTH  60

#define  FALSE  0
#define  TRUE   1

main(int argc, char ** argv) {
    int  i;
    int  n;
    int  input_seed, actual_seed;
    int  c_flag = FALSE;
    int  i_flag = FALSE;
    int  n_flag = FALSE;
    int  o_flag = FALSE;
    int  s_flag = FALSE;

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    real coeff = 1, scale = 1, exponent = 1;
    vector center = 0;

    real r_max = 1;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "A:a:c:in:oR:r:s:x:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'A': coeff = atof(poptarg);
		      break;
	    case 'a':
	    case 'R':
		      scale = atof(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = TRUE;
		      break;
	    case 'n': n_flag = TRUE;
		      n = atoi(poptarg);
		      break;
	    case 'o': o_flag = TRUE;
                      break;
	    case 'r': r_max = atof(poptarg);
		      break;
	    case 's': s_flag = TRUE;
		      input_seed = atoi(poptarg);
		      break;
	    case 'x': exponent = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    if (n_flag == FALSE) {
        cerr << "makepowerlaw: must specify the number # of";
	cerr << " particles with -n#" << endl;
	exit(1);
    }
    
    if (n < 1) {
        cerr << "makepowerlaw: n < 1 not allowed" << endl;
	exit(1);
    }

    if (exponent <= 0 || exponent >= 3) {
        cerr << "makepowerlaw: 0 < x < 3 required" << endl;
	exit(1);
    }

    // Create a linked list of (labeled) nodes.

    dyn *b, *by, *bo;

    b = new dyn();
    if (i_flag) b->set_label("root");

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

    // Add comments, etc.

    if (c_flag == TRUE) b->log_comment(comment);

    b->log_history(argc, argv);

    if (s_flag == FALSE) input_seed = 0;
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makepowerlaw: random seed = " << actual_seed << endl;

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    // Create dyn data.

    makepowerlaw(b, n, coeff, scale, exponent, r_max);

    // Flag actions for use by kira.

    putrq(b->get_log_story(), "kira_pl_coeff", coeff);
    putrq(b->get_log_story(), "kira_pl_exponent", exponent);
    putrq(b->get_log_story(), "kira_pl_scale", scale);
    putvq(b->get_log_story(), "kira_pl_center", center);

    putiq(b->get_log_story(), "ignore_internal", 1);

    put_node(cout, *b);
}

#endif

/* end of: makepowerlaw.c */
