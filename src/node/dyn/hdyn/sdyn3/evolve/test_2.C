
// test_2.C: Determine merger cross-sections for two-body collisions.

// Starlab application:  get_sigma3.

#include "sigma3.h"
#define DEBUG 0

#define Grav	6.673e-8
#define Msun	1.989e33
#define AU	1.486e13
#define Rsun	6.960e10

#define Kms	1e5
#define SECONDS_IN_YEAR 3.16e7
#define CM_IN_PC 3.086e18

#define N_HE   10
#define MIN_HE 0.24	// MIN_HE and MAX_HE determine limits on mass for
#define MAX_HE 0.39	// binning purposes only.

#define N_MASS 24
#define M_MIN  0.0	// M_MIN and M_MAX determine limits on mass for binning
#define M_MAX  2.4	// purposes only.  Actual limits are m_lower, m_upper.

static real m_lower  = 0.1;
static real m_upper  = 0.8;
static real exponent = 1.35;	// Convention: Salpeter = +1.35!
static real m_unit;

// Interpolation from real stellar models:

#define N_INT 9
static real m_int[N_INT] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8};
static real r_int[N_INT] = {0.1, 0.2, 0.3, 0.4, 0.51, 0.63, 0.75, 1.09, 1.87};
static real y_int[N_INT] = {0.240, 0.240, 0.243, 0.248, 0.254,
			    0.273, 0.313, 0.347, 0.385};

static int n_scatt = 0;
static int n_merge = 0;
static int temp_counter[N_MASS][N_HE][N_RHO_ZONE_MAX];

static int total_counter[N_MASS][N_HE][N_RHO_ZONE_MAX];
static real sigma[N_MASS][N_HE];
static real sigma_err_sq[N_MASS][N_HE];


local int get_index(real m)
{
    for (int i = 0; i < N_INT; i++) if (m_int[i] > m) return i-1;
    return N_INT - 1;
}


local real get_quantity(real m, real* q)
{
    int i = get_index(m);

    if (i < 0)
	return q[0];
    else if (i >= N_INT-1)
	return q[N_INT-1];

    return q[i] + (m - m_int[i]) * (q[i+1] - q[i]) / (m_int[i+1] - m_int[i]);
}


local real radius(real m)	// Mass-radius relation (units: solar)
{
    return get_quantity(m, r_int);
}


local real helium(real m)	// Helium fraction (mass unit = solar mass)
{
    return get_quantity(m, y_int);
}


// get_mass: Return a random mass between ml and mu, distributed as a
//	     power law with exponent x (i.e. dN/dm \propto m^x).

local real get_mass(real ml, real mu, real x)
{
    real r = randinter(0,1);

    if (x == -1) 
	return ml * exp( r*log(mu/ml) );
    else
	return ml * pow( 1 + r * (pow(mu/ml, 1+x) - 1), 1.0/(1+x));
}


local void zero_temp_counter()
{
    for (int i = 0; i < N_MASS; i++)
	for (int j = 0; j < N_HE; j++)
	    for (int k = 0; k < N_RHO_ZONE_MAX; k++) temp_counter[i][j][k] = 0;
}


local void initialize_stats()
{
    for (int i = 0; i < N_MASS; i++)
	for (int j = 0; j < N_HE; j++) {
	    sigma[i][j] = sigma_err_sq[i][j] = 0;
	    for (int k = 0; k < N_RHO_ZONE_MAX; k++) total_counter[i][j][k] = 0;
	}
}


// accumulate_stats: Gather statistics supplied by get_sigma3. Customize here
//                   to the specific application at hand.

local void accumulate_stats(scatter_profile& prof,
			    initial_state3& init,
			    intermediate_state3& inter,
			    final_state3& final,
			    int rho_bin)
{

    n_scatt++;

    real mmerge = 0;
    real ymerge;

    if (final.descriptor == merger_binary_1
	|| final.descriptor == merger_escape_1) {	 // 2 and 3 merged

	mmerge = prof.m2 + prof.m3;
	ymerge = (prof.m2 * helium(prof.m2*m_unit)
		  + prof.m3 * helium(prof.m3*m_unit)) / mmerge;

    } else if (final.descriptor == merger_binary_2
	       || final.descriptor == merger_escape_2) { // 1 and 3 merged

	real m1 = 1 - prof.m2;
	mmerge = m1 + prof.m3;
	ymerge = (m1 * helium(m1*m_unit)
		  + prof.m3 * helium(prof.m3*m_unit)) / mmerge;

    } else if (final.descriptor == merger_binary_3
	       || final.descriptor == merger_escape_3) { // 1 and 2 merged

	real m1 = 1 - prof.m2;
	mmerge = 1;
	ymerge = m1 * helium(m1*m_unit) + prof.m2 * helium(prof.m2*m_unit);

    } else if (final.descriptor == triple_merger) {	 // 1, 2 and 3 merged

	real m2 = prof.m2;
	real m1 = 1 - m2;
	mmerge = 1 + prof.m3;
	ymerge = (m1 * helium(m1*m_unit) + m2 * helium(m2*m_unit)
		  + prof.m3 * helium(prof.m3*m_unit)) / mmerge;
    }

    if (mmerge > 0) {
	mmerge *= m_unit;

	int i = (int) (N_MASS * (mmerge - M_MIN) / (M_MAX - M_MIN));
	if (i < 0) i = 0;
	if (i >= N_MASS) i = N_MASS - 1;

	int j = (int) (N_HE * (ymerge - MIN_HE) / (MAX_HE - MIN_HE));
	if (j < 0) j = 0;
	if (j >= N_HE) j = N_HE - 1;

	n_merge++;

	temp_counter[i][j][rho_bin]++;
    }
}

// normalize_and_store_counts: convert raw counts into cross-sections and
//			       errors, and add them into the total

local void normalize_and_store_counts(sigma_out& out, real weight)
{
    for (int i = 0; i < N_MASS; i++)
	for (int j = 0; j < N_HE; j++) {
	    real rho_sq_min, rho_sq_max;
	    int k;
	    for (k = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
		 k <= out.i_max;
		 k++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

		real factor = weight * (rho_sq_max - rho_sq_min)
		                     / out.trials_per_zone;
		sigma[i][j] += factor * temp_counter[i][j][k];
		sigma_err_sq[i][j] += factor * factor * temp_counter[i][j][k];
		total_counter[i][j][k] += temp_counter[i][j][k];
	    }
	}
} 

local void print_stats(sigma_out& out)
{
    int i, j;

    cerr << "\nDifferential cross sections (n_scatt = " << n_scatt << "):\n";

    cerr << "\nRaw counts (row = 10*mass, column = scaled He abundance"
	 << " (total = " << n_merge << "):\n\n";
    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE; j++) {
	    int total = 0;
	    for (int k = 0; k <= out.i_max; k++) 
		total += total_counter[i][j][k];
	    fprintf(stderr, " %7d", total);
	}
	cerr << endl;
    }

    cerr << "\nsigma / (pi a^2):\n\n";
    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE; j++)
	    fprintf(stderr, " %7.3f", sigma[i][j]);
	cerr << endl;
    }

    cerr << "\n(sigma error) / (pi a^2):\n\n";
    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE; j++)
	    fprintf(stderr, " %7.3f", sqrt(sigma_err_sq[i][j]));
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

    real sma = 1.0;	// Binary semi-major axis in A.U.
    real v_inf = 10.0;	// Relative velocity at infinity in km/s

    int n_sigma3 = 1;
    bool sig_print = FALSE;
    bool stat_print = 0;

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    const char *param_string ="a:A:c:C:D:d:e:l:L:n:N:p:PqQs:u:U:v:x:"; 

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'a': sma = atof(poptarg);
		      break;
	    case 'A': prof.eta = atof(poptarg);
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': max_trial_density = atof(poptarg);  // Max density for a
		      break;				  // single experiment
	    case 'D': dt_snap = atof(poptarg);
		      scatter_summary_level = 2;  // Undo with later "-q/Q"
		      break;
	    case 'e': prof.ecc = atof(poptarg);
		      prof.ecc_flag = 1;
		      break;
	    case 'l':
	    case 'L': m_lower = atof(poptarg);
		      break;
	    case 'n': n_sigma3 = atoi(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'p': stat_print = atoi(poptarg);
		      if (stat_print <= 0) stat_print = 1;
		      break;
	    case 'P': sig_print = 1 - sig_print;
		      stat_print = 1;
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
	    case 'u':
	    case 'U': m_upper = atof(poptarg);
		      break;
	    case 'v': v_inf = atof(poptarg);
		      break;
	    case 'V': debug = atoi(poptarg);
		      break;
	    case 'x': exponent = atof(poptarg);
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

    cerr << "\n\t*** "
	 << "Binary/single-star scattering with stellar collisions."
	 << " ***\n\n";
    cerr << "Binary semi-major-axis a = " << sma
	 << " A.U.  Velocity at infinity = " << v_inf << " km/s.\n";
    cerr << "Mass range = " << m_lower << "-" << m_upper
	 << " solar masses,  power-law exponent = " << exponent
	 << " (Salpeter = 1.35)\n";
    cerr << "Target trial density = " << max_trial_density
	 << ",  total number of mass combinations = " << n_sigma3 << endl;

    int first_seed = srandinter(seed, n_rand);
    cerr << "Random seed = " << first_seed << endl;

    sigma_out out;
    initialize_stats();

    // Loop over mass combinations, obtaining cross-sections for each.

    for (int i = 0; i < n_sigma3; i++) {

	// Get masses in solar masses, randomly from the mass function.

	real m1 = get_mass(m_lower, m_upper, -exponent - 1);
	real m2 = get_mass(m_lower, m_upper, -exponent - 1);
	real m3 = get_mass(m_lower, m_upper, -exponent - 1);

	// Scale masses, v_inf and radii:

	m_unit = m1 + m2;
	real r_unit = sma;
	real v_unit = sqrt( (Grav * Msun / AU)
			   * m1 * m2 * (m1 + m2 + m3) / ((m1 + m2) * m3)
			   / sma ) / Kms;

	prof.m2 = m2 / m_unit;
	prof.m3 = m3 / m_unit;
	prof.v_inf = v_inf / v_unit;

	prof.r1 = radius(m1) * Rsun / (sma * AU);
	prof.r2 = radius(m2) * Rsun / (sma * AU);
	prof.r3 = radius(m3) * Rsun / (sma * AU);

	int p = cerr.precision(STD_PRECISION);
	cerr << "\nMass combination #" << i+1 << ":\n";
	cerr << "    m1 = " << m1 << "  m2 = " << m2
	     << "  m3 = " << m3 << " solar\n";
	cerr << "    r1 = " << radius(m1) << "  r2 = " << radius(m2)
	     << "  r3 = " << radius(m3) << " solar\n";
	cerr << "    a = " << sma << " A.U. = "
	     << sma * AU / Rsun << " solar,";
	cerr << "  v_inf = " << v_inf << " km/s = "
	     << prof.v_inf << " scaled\n";

	// Get cross-sections. Note that the primary use of get_sigma3
	// here is to ensure that properly sampled data are sent to
	// routine accumulate_stats.

	zero_temp_counter();
	get_sigma3(max_trial_density, prof, out,
		   debug, cpu_time_check,
		   dt_snap, snap_cube_size,
		   scatter_summary_level,
		   accumulate_stats);

	// We are assuming here that all mass combinations carry
	// equal statistical weight.

	normalize_and_store_counts(out, 1.0/n_sigma3);

	// Optional diagnostic/intermediate output:

	if (sig_print) print_sigma3(out, prof.v_inf * prof.v_inf);
	if (stat_print > 0 && (i+1)%stat_print == 0) print_stats(out);

	cerr.precision(p);
    }

    print_stats(out);
}
