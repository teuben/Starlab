
// all_sig.C:	Generate cross sections and statistics on point-mass
//		encounters for identical binary components and
//		arbitrary incoming mass.
//
//		Determine:	(1) Total cross-sections, broken down by type:
//					0 flybys (energy exchange only)
//					1 exchange
//					2 hierarchical resonance
//					3 "strong" resonance (n_ocs >= MIN_OSC)
//					4 "weak" resonant preservation
//					5 "weak" resonant exchange
//				(2) Differential cross-sections in
//					energy exchange dE/E
//					final a and e (2-D),
//				    both broken down by type

// Starlab application:  get_sigma3.

#include "sigma3.h"

#define MIN_OSC 4

// Global statistics data:

#define NTYPE		6	// number of defined outcome types

#define DE_MIN		0.01	// starting point for dE/E binning
#define EBINS_PER_DECADE  5	// number of bins per decade in dE/E
#define N_DECADES	3	// number of decades spanned in dE/E
#define NE	EBINS_PER_DECADE * N_DECADES	// number of bins in dE/E

#define A_MIN_DEF	0.17677669529663688110	// = sqrt(2)/8
real a_min = A_MIN_DEF;
				// starting point for final sma binning
#define ABINS_PER_DOUBLING  2	// number of bins per doubling in sma
#define NSMA		10	// number of bins in sma 

#define NECC		5	// number of bins in final eccentricity^2

#define NR	N_RHO_ZONE_MAX	// (for brevity)

// Counters:

static int n_total[NTYPE][NR];
static int ne[NE+1][NTYPE][NR];
static int nae[NSMA][NECC][NTYPE][NR];

static int total_trials[NR];	// Keep track of progress locally.
static int total_total;

static bool print_counts = FALSE;
static int output_frequency = -1;

// Initialize all counters:

local void initialize_stats()
{
    for (int type = 0; type < NTYPE; type++)
	for (int kr = 0; kr < NR; kr++) {
	    n_total[type][kr] = 0;
	    for (int je = 0; je < NE; je++) ne[je][type][kr] = 0;
	    for (int ja = 0; ja < NSMA; ja++)
		for (int ke = 0; ke < NECC; ke++)
		    nae[ja][ke][type][kr] = 0;
	}

    for (int kr = 0; kr < NR; kr++) total_trials[kr] = 0;
    total_total = 0;
}

// print_stats:  Convert raw counts into final data and print them out.
//	         Weights are provided by a library routine to minimise
//		 duplication in functionality and to allow the internal
//		 workings of the package to change without affecting
//		 user-written code.

local void print_stats(scatter_profile& prof, sigma_out& out)
{
    char* s[NTYPE] = {"non-resonant preservations   ",
		      "non-resonant exchanges       ",
		      "hierarchical resonances      ",
		      "strong democratic resonances ",
		      "weak democratic preservations",
		      "weak democratic exchanges    "};

    real v2 = prof.v_inf * prof.v_inf;

    cerr << "\n--------------- CROSS-SECTIONS BY TYPE ---------------\n";

    for (int type = 0; type < NTYPE; type++) {

	cerr << endl << "**** " << s[type] << "\n\n";

	real sigma, error, dsig[NE], derr[NE],
	     dsig2[NSMA][NECC+1], derr2[NSMA][NECC+1];
	int je, ja, ke;

	sigma = error = 0;

	for (je = 0; je < NE; je++) dsig[je] = derr[je] = 0;

	for (ja = 0; ja < NSMA; ja++)
	    for (ke = 0; ke <= NECC; ke++) dsig2[ja][ke] = derr2[ja][ke] = 0;

	for (int kr = 0; kr <= out.i_max; kr++) {

	    real w = zone_weight(out, kr);  // = zone_area / trials_per_zone
	    real w2 = w * w;

	    sigma += w * n_total[type][kr];
	    error += w2 * n_total[type][kr];

	    for (je = 0; je < NE; je++) {
		dsig[je] += w * ne[je][type][kr];
		derr[je] += w2 * ne[je][type][kr];
	    }

	    for (ja = 0; ja < NSMA; ja++)
		for (ke = 0; ke < NECC; ke++) {
		    dsig2[ja][ke] += w * nae[ja][ke][type][kr];
		    derr2[ja][ke] += w2 * nae[ja][ke][type][kr];
		    dsig2[ja][NECC] += w * nae[ja][ke][type][kr];
		    derr2[ja][NECC] += w2 * nae[ja][ke][type][kr];
		}
	}

	cerr << "v^2 sigma = " << v2 * sigma << endl;

	if (sigma > 0) {

	    cerr << "v^2 error = " << v2 * sqrt(error) << endl;

	    cerr << "\nv^2 dsig(E)/dE (dE_min = " << DE_MIN
		 <<  ", " << EBINS_PER_DECADE << " bins per decade):\n\n";
	    int je1, je2;
	    real e_fac = pow(10.0, 1.0/EBINS_PER_DECADE);
	    for (je1 = 0, je = 0; je1 < N_DECADES; je1++) {
		real e1 = DE_MIN, e2 = DE_MIN;
		for (je2 = 0; je2 < EBINS_PER_DECADE; je2++, je++) {
		    e2 *= e_fac;
		    fprintf(stderr, " %11.7f", v2*dsig[je]/(e2-e1));
		    e1 = e2;
		}
		cerr << endl;
	    }

	    cerr << "\nv^2 derr(E)/dE:\n\n";
	    for (je1 = 0, je = 0; je1 < N_DECADES; je1++) {
		real e1 = DE_MIN, e2 = DE_MIN;
		for (je2 = 0; je2 < EBINS_PER_DECADE; je2++, je++) {
		    e2 *= e_fac;
		    fprintf(stderr, " %11.7f", v2*sqrt(derr[je])/(e2-e1));
		    e1 = e2;
		}
		cerr << endl;
	    }

	    cerr << "\n" << NECC
		 << " v^2 dsig2(a,e) / v^2 dsig2(a)" << endl
		 << "(a -->, a_min = " << a_min
		 << ", " << ABINS_PER_DOUBLING
		 << " bins per doubling; linear in e^2):\n\n";

	    for (ke = 0; ke <= NECC; ke++) {
		for (ja = 0; ja < NSMA; ja++) {

		    real average = v2*dsig2[ja][NECC] / NECC;
		    real dsig = v2*dsig2[ja][ke];
		    if (ke < NECC && average > 0) dsig /= average;

		    fprintf(stderr, "%7.3f", min(99.999, dsig));
		}

		if (ke == 0)      fprintf(stderr, "    e = 0");
		if (ke == NECC-1) fprintf(stderr, "    e = 1\n");
		if (ke == NECC)   fprintf(stderr, "    total");
		cerr << endl;
	    }

	    cerr << "\n" << NECC << " v^2 derr2(a,e) / v^2 dsig2(a):\n\n";

	    for (ke = 0; ke <= NECC; ke++) {
		for (ja = 0; ja < NSMA; ja++) {

		    real average = v2*dsig2[ja][NECC] / NECC;
		    real derr = v2*sqrt(derr2[ja][ke]);
		    if (ke < NECC && average > 0) derr /= average;
		    
		    fprintf(stderr, "%7.3f", min(99.999, derr));
		}
		if (ke == NECC-1) cerr << endl;
		cerr << endl;
	    }
	}
    }
    cerr << endl;
}

local real log2(real x)
{
    return log(x)/log(2.0);
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
    // Divide into types and determine binning:

    if (inter.n_stars == 3
	&& inter.descriptor != unknown_intermediate
	&& final.descriptor != error
	&& final.descriptor != stopped
	&& final.descriptor != unknown_final) {

	// Type:

	int type = -1;
	if (inter.descriptor == non_resonance) {
	    if (final.descriptor == preservation)
		type = 0;
	    else
		type = 1;
	} else if (inter.descriptor == hierarchical_resonance)
	    type = 2;
	else if (inter.descriptor == democratic_resonance) {
	    if (inter.n_osc >= MIN_OSC)
		type = 3;
	    else {
		if (final.descriptor == preservation)
		    type = 4;
		else
		    type = 5;
	    }
	}

	// Energy change (de = dE/E):

	real m1 = 1 - init.m2;
	real m2 = init.m2;

	real e_init = 0.5*m1*m2 / init.sma;

	if (final.escaper == 1)
	    m1 = init.m3;
	else if (final.escaper == 2)
	    m2 = init.m3;
	real e_final = 0.5*m1*m2 / final.sma;

	real de = e_final/e_init - 1;

	// Use logarithmic binning in dE/E (starting at DE_MIN).

	int je = -1;
	if (de >= DE_MIN)
	    je = min(NE, (int) (EBINS_PER_DECADE * log10(de/DE_MIN)));

	// Use log_2 binning in final sma.

	int ja = (int) (ABINS_PER_DOUBLING * log2(final.sma/a_min));

	// Use linear binning in final ecc^2.

	int ke = (int) (NECC * final.ecc * final.ecc);

	// Update counters:

	if (type >= 0) {
	    n_total[type][rho_zone]++;
	    if (je >= 0) ne[je][type][rho_zone]++;
	    if (ja >= 0 && ja < NSMA && ke >= 0 && ke < NECC)
		nae[ja][ke][type][rho_zone]++;
	}

	// Diagnostics:
/*
	int p = cerr.precision(STD_PRECISION);
	cerr << "accumulate_stats: " << state_string(inter.descriptor)
	     << " " << state_string(final.descriptor)
	     << ",  n_osc = " << inter.n_osc << endl;
	cerr << "  rho_zone = " << rho_zone;
	cerr << ",  je, ja, ke = " << je << "  " << ja << "  " << ke << endl;
	cerr << "  escaper = " << final.escaper;
	cerr << "  mass = " << final.system[final.escaper-1].mass;
	cerr << "  e_init = " << e_init;
	cerr << "  de/e = " << de << endl;
	cerr.precision(p);
*/
    }

    // Finally, check to see if any intermediate output is required.

    if (output_frequency >= 0) {

	total_trials[rho_zone]++;
	total_total++;

	if (total_total >= output_frequency) {

	    // Check that all zones have the same number of trials.

	    int kr;
	    for (kr = 1; kr <= out.i_max; kr++)
		if (total_trials[kr] != total_trials[0]) return;

	    // Reset local counter and print cross-sections.

	    total_total = 0;

	    if (rho_zone != out.i_max)
		cerr << "\naccumulate_stats:  Warning: not at end of row!\n";

	    // We can only reach this point when we have just completed
	    // a "row" of trials (so we should have rho_zone = out.i_max).
	    // However, control will be transferred from the scatter driver
	    // here *before* out.trials_per_zone is updated, so we must
	    // temporarily do it here to make the intermediate cross
	    // sections come out right.

	    // NOTE: (1) The "trial density" is updated at a higher level
	    //		 within sigma3, so it will be wrong in these
	    //		 intermediate tables.
	    //	     (2) If we happen to get output just as the outer zone
	    //		 has been expanded and is being built up, the
	    //		 cross sections will be *wrong* (but counts are OK)
	    //		 because, in that rare case, trials_per_zone should
	    //		 not be incremented here...

	    out.trials_per_zone++;	// <-----

	    counts_to_sigma(out);
	    print_sigma3(out, prof.v_inf * prof.v_inf);
	    if (print_counts) print_all_sigma3_counts(out, cerr);

	    // User-specific output:

	    print_stats(prof, out);

	    out.trials_per_zone--;	// <-----

	}	
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// No particular need to modify anything below this line.

//----------------------------------------------------------------------
//----------------------------------------------------------------------

main(int argc, char **argv)
    {
    int  debug  = 1;
    int  seed 	= 0;
    int  n_rand = 0;
    real max_trial_density = 1.0;

    real cpu_time_check = 3600; // One check per CPU hour!
    real dt_snap = VERY_LARGE_NUMBER;
    real snap_cube_size = 10;
    int scatter_summary_level = 0;

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    const char *param_string = "a:A:c:C:d:D:e:m:M:N:o:pqQs:v:V:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'a': a_min = atof(poptarg);
		      break;
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
	    case 'o': output_frequency = atoi(poptarg);
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

    print_stats(prof, out);
}
