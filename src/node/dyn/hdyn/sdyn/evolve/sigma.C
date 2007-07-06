//// sigma:   Determine all cross-sections for N-body scattering.
////          Use scatter in a Monte-Carlo fashion to compute cross-
////          sections of interest given the overall parameters of the
////          binary-single star interaction.
////
////          Units: G = 1, binary mass = 1, binary semi-major axis = 1.
////
////
//// Options:   -A    specify accuracy parameter [0.02]
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
////            -t    specify time span of integration [infty]
////            -v    specify incomer velocity at infinity [0 ==> Etot = 0]
////            -V    maximize output [false]
////            -x    specify primary radius [0]
////            -y    specify secondary radius [0]
////            -z    specify incomer radius [0]

// Starlab application:  get_sigma.

#include "sigma_MPI.h"

#ifdef TOOLBOX

main(int argc, char **argv) {

  sigma_input input;

  // identical binary collision
  //  char* default_init  
  //       = "-M 1 -v 1 -rm 3 -t -r1 0 -r2 0 -q 1 -p -a 1 -q 1 -r1 0 -r2 0";
  char *default_init; 
  //  default_init 
  //    = "-M 0.66 -rm 3 -v 0.0358 -t -q 0.43 -r1 0.0322 -r2 0.0217 -p -a 5.1190 -q 0.84 -r1 0.0278 -r2 0.0172";

   // = "M 0.92 -rm 3 -v 0.031631 -t -q 0.56 -r1 0.0638 -r2 0.0383 -p -a 9.2571 -q 1.0 -r1 0.0446 -r2 0.0510";



  //  strcpy(&input.init_string[0], default_init);

  //    check_help();
  
  real  delta_t = VERY_LARGE_NUMBER;       // time span of the integration
  real  dt_out = VERY_LARGE_NUMBER;       // time output interval
  
  extern char *poptarg;
  int c;
  char* param_string = "A:c:C:d:D:e:g:Ii:m:M:N:pqQs:t:v:V:";
  
  while ((c = pgetopt(argc, argv, param_string,
		  "$Revision$", _SRC_)) != -1)
    switch(c) {
    case 'A': input.eta = atof(poptarg);
      break;
    case 'c': input.cpu_time_check = 3600*atof(poptarg);
      // (Specify in hours)
      break;
    case 'C': input.snap_cube_size = atof(poptarg);
      break;
    case 'd': input.max_trial_density = atof(poptarg);
      break;
    case 'D': input.dt_out = atof(poptarg);
#if 0
    case 'D': if (!pvm) {
      input.dt_snap = atof(poptarg);
    } else
      cerr << "\"-D\" option disallowed in PVM mode\n";
      break;
#endif
      //	    case 'e': input.ecc = atof(poptarg);
      //		      input.ecc_flag = 1;
      //		      break;
    case 'g': input.tidal_tol_factor = atof(poptarg);
      break;
      //case 'I': intermediate_sigma = 1 - intermediate_sigma;
      //break;
      //    case 'i': strcpy(input.init_string, poptarg);
      //    case 'i': strcpy(&default_init[0], poptarg);
    case 'i': default_init = poptarg;
      break;
    case 'M': input.pmass = atof(poptarg);
      break;
    case 'm': input.pmass = atof(poptarg);
      break;
    case 'N': input.n_rand = atoi(poptarg);
      break;
      // case 'p': print_counts = 1 - print_counts;
      // break;
#if 0
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
#endif
    case 's': input.seed = atoi(poptarg);
      break;
    case 't': input.delta_t = atof(poptarg);
      break;
    case 'v': input.v_inf = atof(poptarg);
      break;
    case 'V': input.debug = 1 - input.debug;
              input.verbose = atoi(poptarg);
      break;
    case '?': params_to_usage(cerr, argv[0], param_string);
      //		      get_help();
    }

    strcpy(&input.init_string[0], default_init);
    
    execute_sigma_experiment(input);
    
}

#endif

