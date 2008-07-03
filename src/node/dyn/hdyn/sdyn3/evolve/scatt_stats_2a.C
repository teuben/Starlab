
// scatt_stats_2a.C: Determine merger cross-sections for 2-body encounters
//		     Steve McMillan, Charles Bailyn (Sep 1994)

// Starlab application:  no library function.

#include "stdinc.h"

#define Grav	6.673e-8
#define Msun	1.989e33
#define AU	1.486e13
#define Rsun	6.960e10

#define Kms	1e5
#define SECONDS_IN_YEAR 3.16e7
#define CM_IN_PC 3.086e18

// Bin final results in:
//	(1) total mass and helium abundance
//	(2) total mass and helium mass.

#define N_MASS 48
#define M_MIN  0.0	// M_MIN and M_MAX determine limits on mass for binning
#define M_MAX  2.4	// purposes only.  Actual limits are m_lower, m_upper.

#define N_HE   30
#define MIN_HE 0.24	// MIN_HE and MAX_HE determine limits on mass for
#define MAX_HE 0.39	// binning purposes only.

#define N_HEMASS 25
#define M_HEMIN  0.0	// M_HEMIN and M_HEMAX determine limits on helium mass
#define M_HEMAX  1.0	// for binning purposes only.

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

static int n_scatt;
static real sigma_1[N_MASS][N_HE];
static real sigma_2[N_MASS][N_HEMASS];
static real sigma_err_sq_1[N_MASS][N_HE];
static real sigma_err_sq_2[N_MASS][N_HEMASS];

// ------------------------------------------------------------------------

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

// ------------------------------------------------------------------------

// Zero out accumulators:

local void initialize_stats()
{
    n_scatt = 0;
    for (int i = 0; i < N_MASS; i++) {
	for (int j = 0; j < N_HE; j++) sigma_1[i][j] = sigma_err_sq_1[i][j] = 0;
	for (j = 0; j < N_HEMASS; j++) sigma_2[i][j] = sigma_err_sq_2[i][j] = 0;
    }
}

local void update_sigma2(real m1, real m2, real r1, real r2, real v_inf)
{
    n_scatt++;

    real mmerge = m1 + m2;
    real rmin = r1 + r2;
    v_inf *= Kms;

    // Express sig in solar units (PI Rsun^2):

    real sigma = rmin * rmin * (1 + 2 * Grav * mmerge * Msun
				         / (rmin * Rsun * v_inf * v_inf));

    real ymerge = (m1 * helium(m1) + m2 * helium(m2)) / mmerge;
    real mhemerge = mmerge * ymerge;

    // Bin the new cross section:

    int i = (int) (N_MASS * (mmerge - M_MIN) / (M_MAX - M_MIN));
    if (i < 0) i = 0;
    if (i >= N_MASS) i = N_MASS - 1;

    int j = (int) (N_HE * (ymerge - MIN_HE) / (MAX_HE - MIN_HE));
    if (j < 0) j = 0;
    if (j >= N_HE) j = N_HE - 1;

    sigma_1[i][j] += sigma;
    sigma_err_sq_1[i][j] += sigma*sigma;

    j = (int) (N_HEMASS * (mhemerge - M_HEMIN) / (M_HEMAX - M_HEMIN));
    if (j < 0) j = 0;
    if (j >= N_HEMASS) j = N_HEMASS - 1;

    sigma_2[i][j] += sigma;
    sigma_err_sq_2[i][j] += sigma*sigma;
}

#define SCALE 100

local void normalize_and_print_stats()
{
    if (n_scatt <= 0) {
	cerr << "\nn_scatt = 0!\n";
	return;
    }
    int i, j;

    real sum1 = 0, sum2 = 0;
    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE;     j++) sum1 += sigma_1[i][j];
	for (j = 0; j < N_HEMASS; j++) sum2 += sigma_2[i][j];
    }

    cerr << "\nNormalized differential cross sections "
	 << "(row = 10*mass, column = scaled He):\n\n";

    cerr << endl << SCALE << " sigma_1 (total = " << sum1/n_scatt
	 << " pi Rsun^2):\n\n";

    sum1 /= SCALE;

    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE; j++)
	    fprintf(stderr, " %7.3f", sigma_1[i][j]/sum1);
	cerr << endl;
    }

    cerr << endl << SCALE << " sigma_2 (total = " << sum2/n_scatt
	 << " pi Rsun^2):\n\n";

    sum2 /= SCALE;

    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HEMASS; j++)
	    fprintf(stderr, " %7.3f", sigma_2[i][j]/sum2);
	cerr << endl;
    }

    cerr << "\nsigma_1 error:\n\n";

    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HE; j++)
	    fprintf(stderr, " %7.3f", sqrt(sigma_err_sq_1[i][j])/sum1);
	cerr << endl;
    }

    cerr << "\nsigma_2 error:\n\n";

    for (i = 0; i < N_MASS; i++) {
	for (j = 0; j < N_HEMASS; j++)
	    fprintf(stderr, " %7.3f", sqrt(sigma_err_sq_2[i][j])/sum2);
	cerr << endl;
    }
}

// ------------------------------------------------------------------------

main(int argc, char **argv)
    {
    int  seed 	= 0;
    int  n_rand = 0;

    real v_inf = 10.0;	// Relative velocity at infinity in km/s
    int  n_sigma2 = 100;

    extern char *poptarg;
    int c;
    const char *param_string = "l:L:n:N:s:u:U:v:x:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'l':
	    case 'L': m_lower = atof(poptarg);
		      break;
	    case 'n': n_sigma2 = atoi(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 'u':
	    case 'U': m_upper = atof(poptarg);
		      break;
	    case 'v': v_inf = atof(poptarg);
		      break;
	    case 'x': exponent = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
        }

    cpu_init();

    cerr << "\nVelocity at infinity = " << v_inf << " km/s.\n";
    cerr << "Mass range = " << m_lower << "-" << m_upper
	 << " solar masses, power-law exponent = " << exponent
	 << " (Salpeter = 1.35)\n";
    cerr << "Total number of mass combinations = " << n_sigma2 << endl;

    int first_seed = srandinter(seed, n_rand);
    cerr << "Random seed = " << first_seed << endl;

    initialize_stats();

    // Loop over mass combinations, obtaining cross-sections for each.

    for (int i = 0; i < n_sigma2; i++) {

	// Get masses in solar masses, randomly from the mass function.

	real m1 = get_mass(m_lower, m_upper, -exponent - 1);
	real m2 = get_mass(m_lower, m_upper, -exponent - 1);

	// Get radii in solar radii.

	real r1 = radius(m1);
	real r2 = radius(m2);

	// Update cross-sections, assuming that all mass combinations
	// carry equal statistical weight.

	update_sigma2(m1, m2, r1, r2, v_inf);
    }

    normalize_and_print_stats();
}
