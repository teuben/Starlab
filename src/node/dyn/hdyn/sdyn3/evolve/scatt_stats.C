
// scatt_stats.C: Determine all cross-sections for three-body scattering.

// Starlab application:  get_sigma3.

#include "sigma3.h"

// Global statistics data:

#define NA 10
#define NE 10

static int n = 0;
static int counter[NA][NE][N_RHO_ZONE_MAX];

static real sigma[NA][NE];
static real sigma_err_sq[NA][NE];

local void initialize_stats()
{
    for (int ia = 0; ia < NA; ia++)
	for (int je = 0; je < NE; je++)
	    for (int kr = 0; kr < N_RHO_ZONE_MAX; kr++)
		counter[ia][je][kr] = 0;
}

// accumulate_stats: Gather statistics supplied by get_sigma3. Customize here
//                   to the specific application at hand.

local void accumulate_stats(scatter_profile& prof,
			    initial_state3& init,
			    intermediate_state3& inter,
			    final_state3& final,
			    int rho_bin,
			    sigma_out& out)
{
    if (final.descriptor == exchange_1 || final.descriptor == exchange_2) {

//	cerr << "acc_stats: " << state_string(inter.descriptor)
//	     << " " << state_string(final.descriptor) << endl;

	if (prof.r3 > 0) {

	    real ratio = final.sma / prof.r3;

	    int ia = (int) (log(ratio) / log(2.0));
	    if (ia >= NA) ia = NA - 1;

	    int je = (int) (10 * final.ecc);
	    if (je >= NE) je = NE - 1;

	    n++;
	    counter[ia][je][rho_bin]++;

//	    cerr << "           a/R = " << ratio
//		 << "  e = " << final.ecc << endl;
	}
    }
}

// normalize_counts: convert raw counts into cross-sections and errors

local void normalize_counts(sigma_out& out)
{
    for (int ia = 0; ia < NA; ia++)
	for (int je = 0; je < NE; je++) {

	    sigma[ia][je] = 0;
	    sigma_err_sq[ia][je] = 0;

	    real rho_sq_min, rho_sq_max;
	    int kr;
	    for (kr = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
		 kr <= out.i_max;
		 kr++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

		real factor = (rho_sq_max - rho_sq_min) / out.trials_per_zone;
		sigma[ia][je] += factor * counter[ia][je][kr];
		sigma_err_sq[ia][je] += factor * factor * counter[ia][je][kr];
	    }
	}
} 

local void print_stats(scatter_profile& prof, sigma_out& out)
{
    normalize_counts(out);
    real v2 = prof.v_inf * prof.v_inf;

    int ia, je;

    cerr << "\nDifferential cross sections:\n";

    cerr << "\nRaw counts (column = log2(sma/r3), row = 10*ecc, total = "
	 << n << "):\n\n";
    for (ia = 0; ia < NA; ia++) {
	for (je = 0; je < NE; je++) {
	    int total = 0;
	    for (int kr = 0; kr <= out.i_max; kr++)
		total += counter[ia][je][kr];
	    fprintf(stderr, " %7d", total);
	}
	cerr << endl;
    }

    cerr << "\nv^2 sigma / (pi a^2)\n\n";
    for (ia = 0; ia < NA; ia++) {
	for (je = 0; je < NE; je++)
	    fprintf(stderr, " %7.3f", v2 * sigma[ia][je]);
	cerr << endl;
    }

    cerr << "\nv^2 (sigma error) / (pi a^2)\n\n";
    for (ia = 0; ia < NA; ia++) {
	for (je = 0; je < NE; je++)
	    fprintf(stderr, " %7.3f", v2 * sqrt(sigma_err_sq[ia][je]));
	cerr << endl;
    }
}

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

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:e:m:M:N:pqQs:v:V:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'A': prof.eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': max_trial_density = atof(poptarg);
		      break;
	    case 'D': dt_snap = atof(poptarg);
		      scatter_summary_level = 2;  // Undo with later "-q/Q"
		      break;
	    case 'e': prof.ecc = atof(poptarg);
		      prof.ecc_flag = 1;
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
		      exit(1);
	    }

    // Debugging:    |debug| = 0 ==> no debugging output
    //               |debug| = 1 ==> output only at end (default)
    //               |debug| = 2 ==> output at end of each top-level iteration
    //        	     |debug| = 3 ==> output after each sub-trial
    //
    //		      debug  < 0 ==> simple statistics also (default: off)

    cpu_init();

    sigma_out out;
    int first_seed = srandinter(seed, n_rand);

    cerr << "Random seed = " << first_seed << endl;
    print_profile(cerr, prof);

    initialize_stats();

    get_sigma3(max_trial_density, prof, out,
	       debug, cpu_time_check,
	       dt_snap, snap_cube_size,
	       scatter_summary_level,
	       accumulate_stats);

    print_sigma3(out, prof.v_inf * prof.v_inf);
    if (print_counts) print_all_sigma3_counts(out, cerr);

    print_stats(prof, out);
}
