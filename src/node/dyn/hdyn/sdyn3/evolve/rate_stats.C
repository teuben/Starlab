
// rate_stats.C: Determine Maxwellian averaged cross-sections for 3-body
//		 scattering.  Input in dimensionless units only.
//
//		 Steve McMillan, Piet Hut, Fred Rasio (July 1994)

// Starlab application:  get_rate3.

#ifdef TOOLBOX
#include "sigma3.h"

#define USAGE "usage: rate_stats \
[-A #] [-d #] [-e #] [-m #] [-M #] [-n #] [-r[123] # ] [-s #] [-v #] \n"

//------------------------------------------------------------------------

// Local statistics data:

#define NA	10
#define NE	10
#define N1	50

static int n_hit = 0;
static int counter[NA][NE][N_RHO_ZONE_MAX];	// Raw counts
static int ecounter[N1][N_RHO_ZONE_MAX];
static int vcounter[N1][N_RHO_ZONE_MAX];

static real sigma[NA][NE];			// Maxwellian averaged sigma
static real sigma_err_sq[NA][NE];
static real esigma[N1];
static real vsigma[N1];
static real v_scale = 1;

local void initialize_local_stats()	// Initialize overall integrals
{
    for (int ia = 0; ia < NA; ia++)
	for (int je = 0; je < NE; je++)
	    sigma[ia][je] = sigma_err_sq[ia][je] = 0;

    for (int i1 = 0; i1 < N1; i1++)
	esigma[i1] = vsigma[i1] = 0;
}

local void initialize_local_counters()	// Initialize counters for specific v
{
    int kr;

    for (int ia = 0; ia < NA; ia++)
	for (int je = 0; je < NE; je++)
	    for (kr = 0; kr < N_RHO_ZONE_MAX; kr++)
		counter[ia][je][kr] = 0;

    for (int i1 = 0; i1 < N1; i1++)
	for (kr = 0; kr < N_RHO_ZONE_MAX; kr++)
	    ecounter[i1][kr] = vcounter[i1][kr] = 0;
}

local real v_recoil(initial_state3& init, final_state3& final)
{
    // Determine the recoil velocity of the final binary.

    real m_total = 1 + init.m3;
    real m2i = init.m2;
    real m1i = 1 - m2i;
    real m3i = init.m3;

    // Find the indices of the remaining bound stars in the "system" array:

    int i1f = 0, i2f = 1;
    if (final.system[0].index == final.escaper) {
	i1f = 2;
    } else if (final.system[1].index == final.escaper) {
	i2f = 2;
    }

    real m1f = final.system[i1f].mass;
    real m2f = final.system[i2f].mass;
    real m3f = m_total - m1f - m2f;

    real de = m1f*m2f/final.sma - m1i*m2i/init.sma;	// N.B. init.sma = 1.
    real vinff = sqrt(max( (2*m_total*de + (m1i+m2i)*m3i*init.v_inf*init.v_inf)
			     / ((m1f+m2f)*m3f), 0.0));

    return m3f*vinff/m_total;
}

// accumulate_local_stats: Gather statistics supplied by get_sigma3.
//                         Customize here to the specific application at hand.

local void accumulate_local_stats(scatter_profile& prof,
				  initial_state3& init,
				  intermediate_state3& inter,
				  final_state3& final,
				  int rho_bin,
				  sigma_out& out)
{
    if (final.descriptor == exchange_1 || final.descriptor == exchange_2) {

	if (prof.r3 > 0) {

	    real ratio = final.sma / prof.r3;

	    int ia = (int) (log(ratio) / log(2.0));
	    if (ia >= NA) ia = NA - 1;

	    int je = (int) (NE * final.ecc);
	    if (je >= NE) je = NE - 1;

	    counter[ia][je][rho_bin]++;
	    n_hit++;

	    // Higher resolution ecentricity:

	    int i1 = (int) (N1 * final.ecc);
	    if (i1 >= N1) i1 = N1 - 1;
	    ecounter[i1][rho_bin]++;

	    // Recoil velocity:

	    i1 = (int) (N1 * v_recoil(init, final) / v_scale);
	    if (i1 >= N1) i1 = N1 - 1;
	    vcounter[i1][rho_bin]++;
	}
    }
}

local void integrate_local_stats(sigma_out& out, real weight)
{
    for (int ia = 0; ia < NA; ia++)
	for (int je = 0; je < NE; je++) {

	    // Normalize the current counts to dsigma and dsigma_err_eq...

	    real dsigma = 0;
	    real dsigma_err_sq = 0;

	    // Standard normalization, from sigma3:

	    real rho_sq_min, rho_sq_max;
	    int kr;
	    for (kr = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
		 kr <= out.i_max;
		 kr++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

		real factor = (rho_sq_max - rho_sq_min) / out.trials_per_zone;
		dsigma += factor * counter[ia][je][kr];
		dsigma_err_sq += factor * factor * counter[ia][je][kr];
	    }

	    // ...then update the integrals (add errors in quadrature).
	    // Integration is the same as in get_rate3 below.

	    sigma[ia][je] += dsigma * weight;
	    sigma_err_sq[ia][je] += dsigma_err_sq * weight * weight;
	}

    for (int i1 = 0; i1 < N1; i1++) {

	// Normalize the current counts to desigma and dvsigma...

	real desigma = 0;
	real dvsigma = 0;

	real rho_sq_min, rho_sq_max;
	int kr;
	for (kr = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
	     kr <= out.i_max;
	     kr++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

	    real factor = (rho_sq_max - rho_sq_min) / out.trials_per_zone;
	    desigma += factor * ecounter[i1][kr];
	    dvsigma += factor * vcounter[i1][kr];
	}

	// Integration:

	esigma[i1] += desigma * weight;
	vsigma[i1] += dvsigma * weight;
    }
}

local void print_local_stats()
{
    int ia, je, i1;

    cerr << "\nDifferential cross sections:\n";

    cerr << "\n<v sigma> / (pi a^2 v_rel_th)\n\n";
    for (ia = 0; ia < NA; ia++) {
	for (je = 0; je < NE; je++)
	    fprintf(stderr, " %7.3f", sigma[ia][je]);
	cerr << endl;
    }

    cerr << "\n(<v sigma> error) / (pi a^2 v_rel_th)\n\n";
    for (ia = 0; ia < NA; ia++) {
	for (je = 0; je < NE; je++)
	    fprintf(stderr, " %7.3f", sqrt(sigma_err_sq[ia][je]));
	cerr << endl;
    }

    cerr << "\n<v sigmae> / (pi a^2 v_rel_th)\n\n";
    for (i1 = 0; i1 < N1; i1++) fprintf(stderr, " %7.3f", esigma[i1]);
    cerr << endl;

    cerr << "\n<v sigmav> / (pi a^2 v_rel_th)\n\n";
    for (i1 = 0; i1 < N1; i1++) fprintf(stderr, " %7.3f", vsigma[i1]);
    cerr << endl;
}

//------------------------------------------------------------------------

local real stat_weight(real v, real v_th_3d) // $v^3 f(v)$
{
    // Return the normalized Maxwellian distribution function at velocity v.

    real x = v/v_th_3d;
    return 3 * sqrt(6/PI) * x * x * x * exp(-1.5*x*x);

    // Normalization is such that $\int v^2 f(v) dv = 1$
    //			      and $\int v^4 f(v) dv = v_{th,3d}^2$
}

local int get_rate3(scatter_profile& prof, real v_rel_th,
		    int nv, real max_dens,
		    real total_sigma[][N_FINAL], real total_err_sq[][N_FINAL])
{
    real v_max = 2 * v_rel_th;
    real dv = v_max / nv;

    int ii, jj;

    // Initialize "standard" counters...

    int total_trials = 0;

    for (ii = 0; ii < N_INTER; ii++)
	for (jj = 0; jj < N_FINAL; jj++) {
	    total_sigma[ii][jj] = 0;
	    total_err_sq[ii][jj] = 0;
	}

    // ...and extra ones:

    initialize_local_stats();

    // Perform the integration using the trapezoid rule in v_inf
    // (skip the end-points at 0 and 2 v_max):

    for (int i = 1; i < nv; i++) {

	prof.v_inf = i * dv;

	initialize_local_counters();

	sigma_out out;
	get_sigma3(max_dens, prof, out,
		   0,			// debug
		   VERY_LARGE_NUMBER,	// cpu_time_check
		   VERY_LARGE_NUMBER,	// dt_snap
		   VERY_LARGE_NUMBER,	// snap_cube_size
		   0,			// scatter_summary_level
		   accumulate_local_stats);

	total_trials += out.total_trials;

	// Accumulate the standard integrals...

	real weight = stat_weight(prof.v_inf, v_rel_th) * dv / v_rel_th;

	for (ii = 0; ii < N_INTER; ii++)
	    for (jj = 0; jj < N_FINAL; jj++) {

		total_sigma[ii][jj] += out.sigma[ii][jj] * weight;

		// Just add errors in quadrature (sqrt taken at end):

		total_err_sq[ii][jj] += out.sigma_err_sq[ii][jj]
		                        * weight * weight;
		}

	// ...and the integrals for differential cross-sections:

	integrate_local_stats(out, weight);

	// (Note that the unit of sigma here is pi a^2 v_rel_th.)
    }

    return total_trials;
}

local real total_coll(real total_sigma[][N_FINAL])
{
    real sigma_coll = 0;	

    for (int i_int = 0; i_int < N_INTER; i_int++)
	sigma_coll += total_sigma[i_int][merger_binary_1]
			    + total_sigma[i_int][merger_binary_2]
			    + total_sigma[i_int][merger_binary_3]
			    + total_sigma[i_int][merger_escape_1]
			    + total_sigma[i_int][merger_escape_2]
			    + total_sigma[i_int][merger_escape_3]
			    + total_sigma[i_int][triple_merger];
    return sigma_coll;
}

main(int argc, char **argv) {

    scatter_profile prof;
    make_standard_profile(prof);

    // Default parameters:

    int  seed 	= 0;
    prof.m2 = prof.m3 = 0.5;
    prof.r1 = prof.r2 = prof.r3 = 0;
    real v_rel_th = 1;			// 3-D rms relative velocity

    real max_dens = 64;			// Best joint choice of parameters?
    int nv = 6;

    // Parse the command line:

    int i = 0;
    while (++i < argc) if (argv[i][0] == '-')
	switch (argv[i][1]) {
	    case 'A': prof.eta = atof(argv[++i]);
		      break;
	    case 'd': max_dens = atoi(argv[++i]);
		      break;
	    case 'e': prof.ecc = atof(argv[++i]);
		      prof.ecc_flag = 1;
		      break;
	    case 'm': prof.m2 = atof(argv[++i]);
		      break;
		      break;
	    case 'M': prof.m3 = atof(argv[++i]);
		      break;
	    case 'n': nv = atoi(argv[++i]);
		      break;
	    case 'r': switch(argv[i][2]) {
			  case '1':	prof.r1 = atof(argv[++i]);
					break;
			  case '2':	prof.r2 = atof(argv[++i]);
					break;
			  case '3':	prof.r3 = atof(argv[++i]);
					break;
		      }
		      if (prof.r3 < 0)
			  err_exit("r3 < 0"); // To avoid bus error in
		      			      // optimized g++ 2.3.3!
		      break;
	    case 's': seed = atoi(argv[++i]);
		      break;
	    case 'v': v_rel_th = atof(argv[++i]);
		      break;
            default:  cerr << USAGE;
		      exit(1);
	}

    if (prof.r3 <= 0) err_exit("rate_stats: must specify r3 > 0.");

    int first_seed = srandinter(seed);
    cerr << "Thermally averaged reaction rates, random seed = "
	 << first_seed << endl;

    // Output information on the scatter profile:

    prof.v_inf = -1;
    print_profile(cerr, prof);

    cerr << "binary eccentricity = ";
    if (prof.ecc_flag)
	cerr << prof.ecc;
    else
	cerr << "thermal [f(e) = 2e]";

    cerr << "   v_rel_th = " << v_rel_th << endl;

    real total_sigma[N_INTER][N_FINAL];
    real total_err_sq[N_INTER][N_FINAL];

    // Calculate the Maxwellian-averaged cross sections:

    v_scale = 4 * v_rel_th;
    int total_trials = get_rate3(prof, v_rel_th, nv, max_dens,
				 total_sigma, total_err_sq);

    // Display the results:

    cerr << "\nmax_dens = " << max_dens
	 << "  nv = " << nv
	 << "  total_trials = " << total_trials
	 << "  (" << n_hit << " hits)\n\n";

    cerr << "<v sigma> (unit = pi a^2 v_rel_th):\n\n";
    print_sigma3_array(total_sigma, 1);

    cerr << "<v sigma_err> (unit = pi a^2 v_rel_th):\n\n";
    print_sigma3_err_array(total_err_sq, 1);

    // More detail:

    cerr << "<v sigma> and <v sigma_err> (details):\n\n";
    print_sigma3_nonmergers(total_sigma, total_err_sq, 1);

    if (total_coll(total_sigma) > 0) 
	print_sigma3_mergers(total_sigma, total_err_sq, 1);

    // Finally, print the differential cross sections:

    print_local_stats();
}
#endif
