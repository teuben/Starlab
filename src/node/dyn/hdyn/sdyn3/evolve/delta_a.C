
// delta_a.C:	Generate statistics on changes in binary semi-major axis
//		following encounters.

// Starlab application:  get_sigma3.

#include "sigma3.h"

#define DEBUG 0

#define MIN_OSC 5	    // Defines a "strong" resonance

// Global statistics data:

#define NBIN 5		    // defined bins:	0: direct exchange
			    //			1: hier_res
			    //			2: strong dem_res
			    //			3: weak dem_res pres
			    //			4: weak dem_res exch

#define NR N_RHO_ZONE_MAX   // for brevity (N_RHO_ZONE_MAX defined in sigma3.h)

// Subdivide hits by bin and radial zone:

static int  n_total[NBIN][NR];
static real da_total[NBIN][NR];

// Keep track of individual scatterings in each bin & zone:

#define NMAX 25000	    // should be adequate...

static int  nscat[NBIN][NR];
static real dascat[NBIN][NR][NMAX];

local void initialize_stats()
{
    for (int ibin = 0; ibin < NBIN; ibin++)
	for (int kr = 0; kr < NR; kr++) {
	    n_total[ibin][kr] = 0;
	    da_total[ibin][kr] = 0;
	    nscat[ibin][kr] = 0;
	}
}

// accumulate_stats: Gather statistics supplied by get_sigma3.
//		     Control will be transferred here after each scattering;
//		     all relevent state information and the rho zone used
//		     are provided by get_sigma3.
//		     Customize here to the specific application at hand.

local void accumulate_stats(scatter_profile& prof,
			    initial_state3& init,
			    intermediate_state3& inter,
			    final_state3& final,
			    int rho_zone,
			    sigma_out& out)
{
    if (inter.n_stars == 3 &&
	(inter.descriptor == hierarchical_resonance
	 || inter.descriptor == democratic_resonance
	 || final.descriptor == exchange_1
	 || final.descriptor == exchange_2)) {

	real da = final.sma/init.sma - 1;

	if (DEBUG) {
	    int p = cerr.precision(STD_PRECISION);
	    cerr << "accumulate_stats: " << state_string(inter.descriptor)
		 << " " << state_string(final.descriptor)
		 << ",  n_osc = " << inter.n_osc << endl;
	    cerr << "  rho_zone = " << rho_zone;
	    cerr << "  escaper = " << final.escaper;
	    cerr << "  mass = " << final.system[final.escaper-1].mass;
	    cerr << "  da/a = " << da << endl;
	    cerr.precision(p);
	}

	int ibin = 0;
	if (inter.descriptor == hierarchical_resonance)
	    ibin = 1;
	else if (inter.descriptor == democratic_resonance) {
	    if (inter.n_osc >= MIN_OSC)
		ibin = 2;
	    else {
		if (final.descriptor == preservation)
		    ibin = 3;
		else
		    ibin = 4;
	    }
	}

	// Store statistical information.

	n_total[ibin][rho_zone]++;
	da_total[ibin][rho_zone] += da;
	dascat[ibin][rho_zone][nscat[ibin][rho_zone]++] = da;
    }
}

local int realcomp(const void * a, const void * b)	// For use by qsort
{
    if (*((real *) a) > *((real *) b))
	return 1;
    else if (*((real *) a) < *((real *) b))
	return -1;
    else
	return 0;
}

// print_stats:  Convert raw counts into final data and print them out.
//	         Weights are provided by a library routine to minimise
//		 duplication in functionality and to allow the internal
//		 workings of the package to change without affecting
//		 user-written code.

local void print_stats(scatter_profile& prof, sigma_out& out)
{
    char* s[NBIN] = {"non-resonant exchange       ",
		     "hierarchical resonance      ",
		     "strong democratic resonance ",
		     "weak democratic preservation",
		     "weak democratic exchange    "};

    for (int ibin = 0; ibin < NBIN; ibin++) {

	real sigma = 0;
	real da = 0;
	int nt = 0;

	cerr << endl << s[ibin] << endl;

	real w_min = VERY_LARGE_NUMBER, w_max = 0;

	for (int kr = 0; kr < out.i_max; kr++) {

	    real w = zone_weight(out, kr);  // = zone_area / trials_per_zone

	    w_min = min(w_min, w);
	    w_max = max(w_max, w);

	    // May as well determine the cross-section too, as a check.

	    sigma += w * n_total[ibin][kr];
	    da += w * da_total[ibin][kr];	// Note: da_total includes
						//       a factor of n_total
	    nt +=nscat[ibin][kr] ;

	    cerr << "    radial zone " << kr << "  (weight = " << w
		 << "):  n_total = " << n_total[ibin][kr]
		 << ",  <da/a> = " << da_total[ibin][kr]
		     			/ max(1, n_total[ibin][kr])
		 << endl;

	    if (nscat[ibin][kr] > 3) {

		cerr << "        da/a quartiles:";

		qsort((void*)&dascat[ibin][kr][0], (size_t)nscat[ibin][kr],
		      sizeof(real), realcomp);
		real factor = 0.25*nscat[ibin][kr];
		for (int i = 1; i <= 3; i++)
		    cerr << "  " << dascat[ibin][kr][(int)(factor*i)];

		cerr << endl;

		if (DEBUG)
		    for (int i = 0; i < nscat[ibin][kr]; i++)
			cerr << "            " << dascat[ibin][kr][i] << endl;

	    }
	}

	if (sigma > 0) da /= sigma;

	cerr << "<da/a> for " << s[ibin] << " = " << da
	     << ",  v^2 sigma = " << prof.v_inf * prof.v_inf * sigma << endl;

	if (nt > 3) {

	    // Combine all radial zones.

	    cerr << "da/a quartiles:";

	    real* tempda = new real[out.i_max*nt*(int)(w_max/w_min + 0.1)];
	    nt = 0;

	    // Note quick and dirty treatment of zone weights...

	    for (int kr = 0; kr <= out.i_max; kr++)
		for (int i = 0; i < nscat[ibin][kr]; i++)
		    for (int ww = 0;
			 ww < (int)(zone_weight(out, kr)/w_min + 0.1); ww++)
			tempda[nt++] = dascat[ibin][kr][i];

	    qsort((void*)tempda, (size_t)nt, sizeof(real), realcomp);
	    real factor = 0.25*nt;
	    for (int i = 1; i <= 3; i++)
		cerr << "  " << tempda[(int)(factor*i)];

	    cerr << endl;

	    if (DEBUG)
		for (int i = 0; i < nt; i++)
		    cerr << "        " << tempda[i] << endl;

	}
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Probably no particular need to modify anything below this line.

//----------------------------------------------------------------------
//----------------------------------------------------------------------

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

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    const char *param_string = "A:c:C:D:d:e:m:M:N:pqQs:V:v:x:y:z:";

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

    initialize_stats();

    // Generic scattering experiment call, but an with additional hook to
    // transfer control to accumulate_stats (above) after each scattering.

    get_sigma3(max_trial_density, prof, out,
	       debug, cpu_time_check,
	       dt_snap, snap_cube_size,
	       scatter_summary_level,
	       accumulate_stats);

    // Out is the standard final descriptor of the scattering experiment.
    // It contains useful statistical and diagnostic information about
    // the details of what went on, but it may not be needed in all cases.

    // Standard "sigma3" output:

    print_sigma3(out, prof.v_inf * prof.v_inf);
    if (print_counts) print_all_sigma3_counts(out, cerr);

    // User-specific output:

    for (int i = 0; i < 72; i++) cerr << "-";
    cerr << endl;

    print_stats(prof, out);
}
