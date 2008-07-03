
// new_stats.C: Determine merger statistics for 3-body scattering.
//		Steve McMillan, Charles Bailyn, Alison Procter (1995).

//		This version doesn't compute cross sections.  Instead,
//		it prints out raw data in a form that can be summed
//		after the fact.

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

static real r_unit;		// Initial binary semi-major axis in A.U.
static real m_unit;		// Initial total binary mass, in solar masses
static real v_unit;		// Initial v_crit, in km/s

static int n_merge;

#define MAX_SAVE 100000

static intermediate_descriptor3	save_inter[MAX_SAVE];
static final_descriptor3	save_final[MAX_SAVE];
static int			save_bin[MAX_SAVE];
static real			save_ecc_0[MAX_SAVE];
static real			save_sma[MAX_SAVE];
static real			save_ecc_1[MAX_SAVE];

// ------------------------------------------------------------------------

static real m_lower  = 0.1;
static real m_upper  = 0.8;
static real exponent = 1.35;	// Convention: Salpeter = +1.35!

// Interpolation from real stellar models:

static int metal = 1;		// Default is solar metallicity

// Low metallicity, for globular clusters:

#define N_INT_0 9
static real m_int_0[N_INT_0]
		= {0.1, 0.2, 0.3, 0.4, 0.5,  0.6,  0.7,  0.75, 0.8};
static real r_int_0[N_INT_0]
		= {0.1, 0.2, 0.3, 0.4, 0.51, 0.63, 0.75, 1.09, 1.87};

// Solar metallicity, for open clusters:

#define N_INT_1 13
static real m_int_1[N_INT_1]
		= {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
		   0.8, 0.9, 1.0, 1.1, 1.2, 1.25};
static real r_int_1[N_INT_1]
		= {0.140, 0.216, 0.292, 0.377, 0.465, 0.555, 0.644,
		   0.733, 0.836, 0.977, 1.182, 1.466, 1.767};

local int get_index(real mass, real* m, int n)
{
    for (int i = 0; i < n; i++) if (m[i] > mass) return i - 1;
    return n - 1;
}

local real get_quantity(real mass, real* m, real* q, int n)
{
    int i = get_index(mass, m, n);

    if (i < 0)
	return q[0];
    else if (i >= n)
	return q[n-1];

    return q[i] + (mass - m[i]) * (q[i+1] - q[i]) / (m[i+1] - m[i]);
}

local real radius(real mass)	// Mass-radius relation (units: solar)
{
    real *r, *m;
    int i, n;

    // Choose which interpolation arrays to use:

    if (metal == 0) {
	m = m_int_0;
	r = r_int_0;
	n = N_INT_0;
    } else {
	m = m_int_1;
	r = r_int_1;
	n = N_INT_1;
    }

    return get_quantity(mass, m, r, n);
}

// ------------------------------------------------------------------------

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

// save_stats: Save essential raw data supplied by get_sigma3.

local void save_stats(scatter_profile& prof,
		      initial_state3& init,
		      intermediate_state3& inter,
		      final_state3& final,
		      int rho_bin,
		      sigma_out& out)
{
    if (   final.descriptor == merger_binary_1
	|| final.descriptor == merger_escape_1
	|| final.descriptor == merger_binary_2
	|| final.descriptor == merger_escape_2
	|| final.descriptor == merger_binary_3
	|| final.descriptor == merger_escape_3
	|| final.descriptor == triple_merger  ) {

	// Save data only in case of merger.

	save_ecc_0[n_merge] = init.ecc;

	save_inter[n_merge] = inter.descriptor;
	save_final[n_merge] = final.descriptor;
	save_bin[n_merge] = rho_bin;

	if (   final.descriptor == merger_binary_1
	    || final.descriptor == merger_binary_2
	    || final.descriptor == merger_binary_3 ) {
	    save_sma[n_merge] = final.sma;
	    save_ecc_1[n_merge] = final.ecc;
	}

	n_merge++;
    }
}

// dump_stats: print the data saved by save_stats.

void dump_stats(scatter_profile& prof, sigma_out& out, int verbose)
{
    if (n_merge > 0) {

	int i;

	// Header information:

	if (verbose) {

	    int p = cerr.precision(4);		// nonstandard precision
	    cerr << "\n    max_zone = " << out.i_max
		 << "  rho_max = " << out.rho_max
		 << endl;

	    cerr << "    trials per zone = " << out.trials_per_zone;
	    cerr << ", total hits per zone:";
	    for (i = 0; i <= out.i_max; i++) cerr << " " << out.n_hit[i];
	    cerr << endl;

	    cerr << "    (v/v_c)^2 * zone weights:";
	    for (i = 0; i <= out.i_max; i++)
		cerr << " " << prof.v_inf*prof.v_inf * zone_weight(out, i);
	    cerr << endl;

	    cerr << endl <<
	 "  total merging cross sections [(v/v_c)^2 * sigma / (pi a^2)]:\n\n";
	    print_sigma3_mergers(out.sigma, out.sigma_err_sq,
			     prof.v_inf*prof.v_inf);
	    cerr.precision(p);

	} else {

	    fprintf(stderr, "\n%s%s%s%s%s%s%s%s%s%s%s\n",
		    "    m1 ",
		    "    m2 ",
		    "    m3 ",
		    "     a_i ",
		    "   e_i ",
		    " bin",
		    " int",
		    " fin",
		    "  d(sigma)",
		    "    a_f ",
		    "   e_f");

	}

	// Details, merger by merger:

	int p = cerr.precision(STD_PRECISION);

	for (i = 0; i < n_merge; i++) {

	    if (verbose) {

		if (i > 0) cerr << endl;
		cerr << "  merger #" << i << endl;

		cerr << "    initial a = "
		     << r_unit << "  e = " << save_ecc_0[i];
		cerr << "   masses: " << m_unit*(1 - prof.m2)
	    		              << " " <<  m_unit*prof.m2
			              << " " <<  m_unit*prof.m3
		     << endl;

		cerr << "    outcome: " << state_string(save_inter[i])
		     << " " << state_string(save_final[i])
		     << endl;

		cerr << "    rho_bin = " << save_bin[i]
		     << "  d(sigma) [AU^2] = "
		     << zone_weight(out, save_bin[i])*PI*r_unit*r_unit
		     << endl;

		if (   save_final[i] == merger_binary_1
		    || save_final[i] == merger_binary_2
		    || save_final[i] == merger_binary_3 )
		    cerr << "    final binary parameters:  a = "
			 << r_unit*save_sma[i]
			 << "  e = " << save_ecc_1[i] << endl;

	    } else {

		fprintf(stderr,
			"%7.3f%7.3f%7.3f%9.3f%7.3f %3d %3d %3d%9.3f",
			m_unit*(1 - prof.m2), m_unit*prof.m2, m_unit*prof.m3,
			r_unit, save_ecc_0[i], save_bin[i],
			save_inter[i], save_final[i],
			zone_weight(out, save_bin[i])*PI*r_unit*r_unit);
		if (   save_final[i] == merger_binary_1
		    || save_final[i] == merger_binary_2
		    || save_final[i] == merger_binary_3 )
		    fprintf(stderr, "%9.3f%7.3f",
			    r_unit*save_sma[i], save_ecc_1[i]);
		fprintf(stderr, "\n");

	    }
	}
	cerr.precision(p);
    }
}

// ------------------------------------------------------------------------

main(int argc, char **argv)
{
    int  verbose  = 1;

    int  debug  = 0;
    int  seed 	= 0;
    int  n_rand = 0;
    real max_trial_density = 1.0;

    real cpu_time_check = 3600; // One check per CPU hour!
    real dt_snap = VERY_LARGE_NUMBER;
    real snap_cube_size = 10;
    int scatter_summary_level = 0;

    real sma = 1.0;	// Initial binary semi-major axis, in A.U.
    real v_inf = 10.0;	// Relative velocity at infinity in km/s

    int n_sigma3 = 1;
    bool sig_print = FALSE;
    bool stat_print = 0;

    scatter_profile prof;
    make_standard_profile(prof);

    extern char *poptarg;
    int c;
    const char *param_string = "a:A:c:C:D:d:e:El:L:Mn:N:p:PqQs:u:U:v:Vx:";

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
	    case 'e': prof.ecc = atof(poptarg);	  // Specify an eccentricity
		      prof.ecc_flag = 1;
		      break;
	    case 'E': prof.ecc_flag = 2;	  // Gaussian eccentricity
		      break;			  // distribution
	    case 'l':
	    case 'L': m_lower = atof(poptarg);
		      break;
	    case 'M': metal = 1 - metal;
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
	    case 'V': verbose = 1 - verbose;
		      break;
	    case 'x': exponent = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	    }

    cpu_init();
    int first_seed = srandinter(seed, n_rand);

    cerr << "\n\t*** "
	 << "Binary/single-star scattering with stellar collisions."
	 << " ***\n\n";
    cerr << "Binary semi-major-axis a = " << sma
	 << " AU, velocity at infinity = " << v_inf << " km/s\n";
    cerr << "Mass range: " << m_lower << " - " << m_upper
	 << " solar, power-law exponent = " << exponent
	 << " (Salpeter = 1.35)\n";
    cerr << "Target trial density = " << max_trial_density
	 << ", total number of mass combinations = " << n_sigma3 << endl;
    cerr << "Random seed = " << first_seed << endl;

    sigma_out out;

    // Loop over mass combinations, obtaining cross-sections for each.

    for (int i = 0; i < n_sigma3; i++) {

	// Get masses in solar masses, randomly from the mass function.

	real m1 = get_mass(m_lower, m_upper, -exponent - 1);
	real m2 = get_mass(m_lower, m_upper, -exponent - 1);
	real m3 = get_mass(m_lower, m_upper, -exponent - 1);

	// Scale masses, v_inf and radii:

	m_unit = m1 + m2;
	r_unit = sma;
	v_unit = sqrt( (Grav * Msun / AU)
		          * m1 * m2 * (m1 + m2 + m3) / ((m1 + m2) * m3)
		          / sma ) / Kms;

	prof.m2 = m2 / m_unit;
	prof.m3 = m3 / m_unit;
	prof.v_inf = v_inf / v_unit;

	prof.r1 = radius(m1) * Rsun / (sma * AU);
	prof.r2 = radius(m2) * Rsun / (sma * AU);
	prof.r3 = radius(m3) * Rsun / (sma * AU);

	if (verbose) {

	    cerr << endl;
	    for (int j = 0; j < 75; j++) cerr << '-';
	    cerr << endl;

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
	    cerr.precision(p);
	}

	// Note that the only use of get_sigma3 here is to ensure
	// that properly sampled data are sent to routine dump_stats.

	n_merge = 0;

	get_sigma3(max_trial_density, prof, out,
		   debug, cpu_time_check,
		   dt_snap, snap_cube_size,
		   scatter_summary_level,
		   save_stats);

	dump_stats(prof, out, verbose);
    }
}
