
//// sigma3:  Determine all cross-sections for three-body scattering.
////          Use scatter3 in a Monte-Carlo fashion to compute cross-
////          sections of interest given the overall parameters of the
////          binary-single star interaction.
////
////          Units: G = 1, binary mass = 1, binary semi-major axis = 1.
////
////
//// Options:   -A    specify accuracy parameter [0.05]
////            -c    specify CPU time check, in hours [1]
////            -C    specify snap cube size [10]
////            -d    specify maximum trial density [1]
////            -D    specify snap output interval [none]
////            -e    specify initial binary eccentricity [thermal]
////            -g    specify tidal tolerance [1.e-6]
////            -I    output intermediate cross sections [true]
////            -m    specify secondary mass (binary mass = 1) [0.5]
////            -M    specify incomer mass [0.5]
////            -N    specify random number count [0]
////            -o    specify outer orbit orientation [random]
////            -p    print raw counts [false]
////            -q    minimal output [false]
////            -Q    intermediate amount of output [false]
////            -s    specify random seed [taken from system clock]
////            -v    specify incomer velocity at infinity [0 ==> Etot = 0]
////            -V    maximize output [false]
////            -x    specify primary radius [0]
////            -y    specify secondary radius [0]
////            -z    specify incomer radius [0]

// Starlab application:  get_sigma3.

#include "sigma3.h"

#ifdef TOOLBOX

main(int argc, char **argv)
{
    int  debug  = 0;
    int  seed 	= 0;
    int  n_rand = 0;
    real max_trial_density = 1.0;

    real cpu_time_check = 3600; // One check per CPU hour!
    real dt_snap = VERY_LARGE_NUMBER;
    real snap_cube_size = 10;
    int scatter_summary_level = 0;

    bool print_counts = FALSE;
    bool intermediate_sigma = TRUE;

    // Find which version we are running:

    bool pvm = false;
    if (strstr(argv[0], ".pvm")) {

#ifndef HAS_PVM					// Compile-time check.
	err_exit("PVM not available");
#endif
	if (getenv("PVM_ROOT") == NULL)
	    err_exit("PVM not available");	// Run time check...

	pvm = true;
    }

    check_help();

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    const char *param_string = "A:c:C:d:D:e:g:Im:M:N:pqQs:v:V:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'A': prof.eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': if (!pvm) 
			  snap_cube_size = atof(poptarg);
	    	      else
			  cerr << "\"-C\" option disallowed in PVM mode\n";
		      break;
	    case 'd': max_trial_density = atof(poptarg);
		      break;
	    case 'D': if (!pvm) {
			  dt_snap = atof(poptarg);
	       		  scatter_summary_level = 2;  // Undo with later "-q/Q"
	    	      } else
			  cerr << "\"-D\" option disallowed in PVM mode\n";
		      break;
	    case 'e': prof.ecc = atof(poptarg);
		      prof.ecc_flag = 1;
		      break;
	    case 'g': prof.tidal_tol_factor = atof(poptarg);
		      break;
	    case 'I': intermediate_sigma = 1 - intermediate_sigma;
		      break;
	    case 'm': prof.m2 = atof(poptarg);
		      break;
	    case 'M': prof.m3 = atof(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'p': print_counts = 1 - print_counts;
		      break;
	    case 'q': if (scatter_summary_level > 0)
		          scatter_summary_level = 0;
		      else
			  scatter_summary_level = 1;
		      break;
	    case 'Q': if (scatter_summary_level > 0)
		          scatter_summary_level = 0;
		      else
			  scatter_summary_level = 2;
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 'v': prof.v_inf = atof(poptarg);
		      break;
	    case 'V': debug = atoi(poptarg);
		      break;
	    case 'x': prof.r1 = atof(poptarg);
		      break;
	    case 'y': prof.r2 = atof(poptarg);
		      break;
	    case 'z': prof.r3 = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
	}

    // Debugging:    |debug| = 0 ==> no debugging output (default)
    //               |debug| = 1 ==> output only at end
    //               |debug| = 2 ==> output at end of each top-level iteration
    //        	     |debug| = 3 ==> output after each sub-trial
    //
    //		      debug  < 0 ==> simple statistics also (default: off)

    cpu_init();

    sigma_out out;
    int first_seed = srandinter(seed, n_rand);

    cerr << "random seed = " << first_seed << endl;
    print_profile(cerr, prof);

    // Note that dt_snap and snap_cube_size should probably be discarded
    // in the parallel implementation.

    int k = 0;
    real init_dens = max_trial_density;

    if (intermediate_sigma) {
	int min_dens = (pvm ? 16 : 4);
	while (init_dens > min_dens) {
	    k++;
	    init_dens /= 4;
	}
    }

    // Determine the cross sections:

    for (real max_dens = init_dens; k >= 0; k--, max_dens *= 4) {
	get_sigma3(max_dens, prof, out,
		   debug, cpu_time_check,
		   dt_snap, snap_cube_size,
		   scatter_summary_level);
	print_sigma3(out, prof.v_inf * prof.v_inf);
	if (print_counts) print_all_sigma3_counts(out, cerr);
    }

}

#endif
