
//// rate3:  Determine Maxwellian averaged cross-sections for 3-body
////         scattering.  Use sigma3 repeatedly, then integrate over
////         a thermal velocity distribution.
////
////         Units: (i)  physical:  solar mass, solar radius for ri,
////                                A.U. for semi-major axis, km/s for v
////                (ii) dynamical: G = 1, binary mass = 1,
////                                binary semi-major axis = 1,
////
//// Options:   -a    specify binary semi-major axis (if using physical
////                     units) [1000 A.U.]
////            -A    specify accuracy parameter [0.05]
////            -c    specify CPU time check, in hours [1]
////            -C    specify snap cube size [10]
////            -d    use physical units [true]
////            -D    specify snap output interval [none]
////            -e    specify initial binary eccentricity [thermal]
////            -g    specify tidal tolerance [1.e-6]
////            -I    output intermediate cross sections [true]
////            -m    specify secondary mass (binary mass = 1) [0.5]
////            -m1   specify primary mass (solar units) [1]
////            -m2   specify secondary mass (solar units) [1]
////            -m3   specify incomer mass (solar units) [1]
////            -mv   specify mass to which v refers [1]
////            -M    specify incomer mass (binary mass = 1) [0.5]
////            -n    specify number of refinements over initial grid [1]
////            -q    minimal output [false]
////            -Q    intermediate amount of output [false]
////            -r1   specify primary radius (binary sma or solar) [0]
////            -r2   specify secondary radius [0]
////            -r3   specify incomer radius [0]
////            -s    specify random seed [taken from system clock]
////            -S    output individual cross-section summaries [false]
////            -v    specify rms velocity at infinity [10 km/s, or 1]

// Starlab library function, with driver application.

#include "sigma3.h"

#ifdef TOOLBOX

#define Grav	6.673e-8
#define Msun	1.989e33
#define Rsun	6.960e10
#define Kms	1e5
#define SECONDS_IN_YEAR 3.16e7
#define CM_IN_PC 3.086e18

// summarize_sigma: Output enough information to restart the current get_sigma.

local void summarize_sigma(real max_dens, scatter_profile & prof,
 			   real cpu_time_check,
                           real dt_snap, real snap_cube_size,
                           int scatter_summary_level)
{
    int p = cerr.precision(STD_PRECISION);
    cerr << "\nsigma3"
	 << " -s " << get_initial_seed()
	 << " -N " << get_n_rand()
	 << " -d " << max_dens;

    cerr.precision(16);				// Uglify the output
    cerr << " -m " << prof.m2
	 << " -M " << prof.m3
         << " -v " << prof.v_inf;
    if (prof.ecc_flag) cerr << " -e " << prof.ecc;
    cerr << " -x " << prof.r1
	 << " -y " << prof.r2
	 << " -z " << prof.r3
	 << " -A " << prof.eta
	 << " -g " << prof.tidal_tol_factor;

    cerr.precision(STD_PRECISION);
    cerr << " -c " << cpu_time_check/3600;	// (specified in hours!)

    if (dt_snap < VERY_LARGE_NUMBER) {
	cerr << " -D " << dt_snap
	     << " -C " << snap_cube_size;
    }

    if (scatter_summary_level == 1)
        cerr << " -q ";
    else if (scatter_summary_level == 2)
        cerr << " -Q ";

    cerr << endl;
    cerr.precision(p);
}

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
		    real total_sigma[][N_FINAL], real total_err_sq[][N_FINAL],
		    real cpu_time_check, real cpu_init, real& cpu_save,
		    real dt_snap, real snap_cube_size,
		    int scatter_summary_level, int sigma_summary_flag)
{
    real v_max = 2 * v_rel_th;
    real dv = v_max / nv;

    int ii, jj;

    // Initialize counters:

    int total_trials = 0;

    for (ii = 0; ii < N_INTER; ii++)
	for (jj = 0; jj < N_FINAL; jj++) {
	    total_sigma[ii][jj] = 0;
	    total_err_sq[ii][jj] = 0;
	}

    // Trapezoid rule in v_inf (skip end-points):

    for (int i = 1; i < nv; i++) {

	prof.v_inf = i * dv;

	real cpu_temp = cpu_time();

	if (sigma_summary_flag)
	    summarize_sigma(max_dens, prof,
			    cpu_time_check,
			    dt_snap, snap_cube_size,
			    scatter_summary_level);

	sigma_out out;
	get_sigma3(max_dens, prof, out,
		   0, cpu_time_check, dt_snap, snap_cube_size,
		   scatter_summary_level);

	int p = cerr.precision(STD_PRECISION);
	if (sigma_summary_flag) {
	    cerr << "\nv_inf = " << prof.v_inf << ":\n";
	    print_sigma3(out, prof.v_inf * prof.v_inf);
	} else
	    cerr << "\nv_inf = " << prof.v_inf
		 << ",  total_trials = " << out.total_trials
		 << ",  delta CPU time = " << cpu_time() - cpu_temp << " s"
		 << endl;
	cerr.precision(p);

	total_trials += out.total_trials;
	real weight = stat_weight(prof.v_inf, v_rel_th) * dv / v_rel_th;

	for (ii = 0; ii < N_INTER; ii++)
	    for (jj = 0; jj < N_FINAL; jj++) {

		total_sigma[ii][jj] += out.sigma[ii][jj] * weight;

		// Just add errors in quadrature (sqrt taken at end):

		total_err_sq[ii][jj] += out.sigma_err_sq[ii][jj]
		                        * weight * weight;
		}

	// (Unit of total here is pi a^2 v_rel_th.)

	// Check the CPU time.  Note that the printed CPU time is the
	// time since this routine was entered.

	if (cpu_time() - cpu_save > cpu_time_check) {
	    cpu_save = cpu_time();
	    cerr << "\nget_rate3:  CPU time = " << cpu_save - cpu_init
		 << "  nv = " << nv << "  i = " << i
		 << "  v_inf = " << prof.v_inf
		 << "  total_trials = " << total_trials
		 << endl << flush;
	}
    }

    return total_trials;
}

local int rescale_physical(real sma, real v_rel_th, real v_unit,
			   real total_sigma[][N_FINAL],
			   real total_err_sq[][N_FINAL])
{
    real scale_factor = PI * sma * sma * v_rel_th      		// to code units
	                   * Rsun * Rsun * v_unit * Kms         // to cm^3 / s
			   * SECONDS_IN_YEAR * 1e9	        // to cm^3 / Gyr
			   / pow(CM_IN_PC, 3);	        	// to pc^3 / Gyr

    real total_max = -1;
    real total_max_err = -1;
    for (int ii = 0; ii < N_INTER; ii++)
	for (int jj = 0; jj < N_FINAL; jj++) {
	    total_sigma[ii][jj] *= scale_factor;
	    total_err_sq[ii][jj] *= scale_factor*scale_factor;
	    if (ii + jj > 0) {
		// (Don't include the top-left element in the array...)
		total_max = max(total_max, total_sigma[ii][jj]);
		total_max_err = max(total_max_err, sqrt(total_err_sq[ii][jj]));
	    }
	}

    // Format the output to maximize significant digits:

    int i_total;
    if (total_max < 1)
	i_total = -1 - (int) (-log10(total_max));
    else
	i_total = (int) (log10(total_max));

    return i_total;
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

#define USAGE "usage: rate3 \
[-a #] [-A #] [-c #] [-C #] [-d] [-D] [-e #] [-g #] \
[-I] [-m[123] #] [-m #] [-mv #] [-M #] [-n #] \
[-q] [-Q] [-r[123] # ] [-s #] [-S] [-v #] \n"

main(int argc, char **argv)
{
    int  seed 	= 0;
    int  n_iter = 1;

    bool intermediate_sigma = TRUE;

    int  scatter_summary_level = 0;	// No low-level output.
    int  sigma_summary_flag = 0;
    bool physical_units_flag = 1;	// Specify input in physical units.

    real cpu_time_check = 3600;		// One check per CPU hour.
    real dt_snap = VERY_LARGE_NUMBER;
    real snap_cube_size = 10;

    char* indent = "12345678901234567890";

    // Default dimensionless parameters:

    real m = 0.5, M = 0.5;
    bool mM_flag = FALSE;		// To check for incorrect input.

    // Default physical parameters:

    real m1 = 1, m2 = 1, m3 = 1;
    bool m123_flag = FALSE;		// To check for incorrect input.
    real sma = 1000;			// 5 A. U. (approx)
    real v_rms, mv = 1;

    real r1 = 0, r2 = 0, r3 = 0;
    bool v_flag = FALSE;		// Defaults depend on context.

    // Units:	mass = 1 solar mass
    //		radius = 1 solar radius
    //		velocity = 1 km/s
    //		binary separation = 1 solar radius

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

    // Parse the command line (the hard way):

    int i = 0;
    while (++i < argc) if (argv[i][0] == '-')
	switch (argv[i][1]) {
	    case 'a': sma = atof(argv[++i]);
		      break;
	    case 'A': prof.eta = atof(argv[++i]);
		      break;
	    case 'c': cpu_time_check = 3600*atof(argv[++i]);  // (In hours)
		      break;
	    case 'C': if (!pvm) 
			  snap_cube_size = atof(argv[++i]);
	    	      else
			  cerr << "\"-C\" option disallowed in PVM mode\n";
		      break;
            case 'd': physical_units_flag = 1 - physical_units_flag;
		      break;
	    case 'D': if (!pvm) {
			  dt_snap = atof(argv[++i]);
	       		  scatter_summary_level = 2; // Undo with later "-q/Q/S"
			  sigma_summary_flag = 1;
	    	      } else
			  cerr << "\"-D\" option disallowed in PVM mode\n";
		      break;
	    case 'e': prof.ecc = atof(argv[++i]);
		      prof.ecc_flag = 1;
		      break;
	    case 'g': prof.tidal_tol_factor = atof(argv[++i]);
		      break;
	    case 'I': intermediate_sigma = 1 - intermediate_sigma;
		      break;
	    case 'm': switch(argv[i][2]) {
			  case '\0':	m = atof(argv[++i]);
					mM_flag = TRUE;
					break;
			  case '1':	m1 = atof(argv[++i]);
					m123_flag = TRUE;
					break;
			  case '2':	m2 = atof(argv[++i]);
					m123_flag = TRUE;
					break;
			  case '3':	m3 = atof(argv[++i]);
					m123_flag = TRUE;
					break;
			  case 'v':	mv = atof(argv[++i]);
					break;
                          default:      cerr << USAGE;
			                get_help();
		      }
		      break;
	    case 'M': M = atof(argv[++i]);
		      mM_flag = TRUE;
		      break;
	    case 'n': n_iter = atoi(argv[++i]);
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
	    case 'r': switch(argv[i][2]) {
			  case '1':	r1 = atof(argv[++i]);
					break;
			  case '2':	r2 = atof(argv[++i]);
					break;
			  case '3':	r3 = atof(argv[++i]);
					break;
		      }
		      if (r3 < 0) err_exit("r3 < 0"); // To avoid bus error in
		      				      // optimized g++ 2.3.3!
		      break;
	    case 's': seed = atoi(argv[++i]);
		      break;
	    case 'S': sigma_summary_flag = 1 - sigma_summary_flag;
		      break;
	    case 'v': v_rms = atof(argv[++i]);
		      v_flag = TRUE;
		      break;
            default:  cerr << USAGE;
		      get_help();
	}

    int first_seed = srandinter(seed);
    cerr << "Thermally averaged reaction rates, random seed = "
	 << first_seed << endl;

    real m_unit, r_unit, v_unit, v_rel_th;

    if (physical_units_flag) {

	// Convert to internal units:

	if (mM_flag) cerr << "Ignoring m/M settings." << endl;
	if (!v_flag) v_rms = 10; // Default

	m_unit = m1 + m2;
	r_unit = sma;
	v_unit = sqrt( (Grav * Msun / Rsun)
			   * m1 * m2 * (m1 + m2 + m3) / ((m1 + m2) * m3)
			   / sma ) / Kms;
	prof.m2 = m2 / m_unit;
	prof.m3 = m3 / m_unit;
	v_rel_th = sqrt( mv/(m1 + m2) + mv/m3 ) * v_rms / v_unit;

	indent = "                 ";
	cerr << "physical units:  m1 = " << m1 << "  m2 = " << m2
	     << "  m3 = " << m3 << "  mv = " << mv << "  (Msolar)\n";
	cerr << indent << "binary semi-major axis = " << sma << "  (Rsolar)\n";
	cerr << indent << "v_rms = " << v_rms
	     	       << "  v_crit = " << v_unit << "  (km/s, 3D)\n";

    } else {

	if (m123_flag) cerr << "Ignoring m[123] settings." << endl;
	if (!v_flag) v_rms = 1; // Default

	m_unit = 1;
	r_unit = 1;
	v_unit = 1;
	prof.m2 = m;
	prof.m3 = M;
	v_rel_th = v_rms;

	indent = "                      ";
	cerr << "dimensionless units:  m = " << m << "  M = " << M << endl;
    }

    prof.r1 = r1 / r_unit;
    prof.r2 = r2 / r_unit;
    prof.r3 = r3 / r_unit;

    cerr << indent << "r1 = " << r1 << "  r2 = " << r2
	 	   << "  r3 = " << r3;
    if (physical_units_flag) cerr << "  (Rsolar)";
    cerr << endl;

    cerr << indent << "binary eccentricity = ";
    if (prof.ecc_flag)
	cerr << prof.ecc << endl;
    else
	cerr << "thermal [f(e) = 2e]\n";

    cerr << indent << "v_rel_th / v_crit = " << v_rel_th << endl;

    prof.v_inf = -1;
    print_profile(cerr, prof);

    cpu_init();
    real cpu_init = cpu_time();
    real cpu_save = cpu_init;

    real max_dens, init_dens = 4;
    int nv, init_nv = 2;

    if (!intermediate_sigma) {
	for (; n_iter > 1; n_iter--, init_dens *= 4, init_nv *= 2);
	n_iter = 1;
    }

    // Loop over pairs of (max_dens, nv):

    for (max_dens = init_dens, nv = init_nv; n_iter > 0;
	 n_iter--, max_dens *= 4, nv *= 2) {	// Best joint choice of params?

	cerr.precision(STD_PRECISION);

	real total_sigma[N_INTER][N_FINAL];
	real total_err_sq[N_INTER][N_FINAL];

	cerr << endl << "max_dens = " << max_dens
	     << ",  nv = " << nv << endl;

        int total_trials = get_rate3(prof, v_rel_th, nv, max_dens,
				     total_sigma, total_err_sq,
				     cpu_time_check, cpu_init, cpu_save,
				     dt_snap, snap_cube_size,
				     scatter_summary_level,
				     sigma_summary_flag);

	cerr << "\ntotal_trials = " << total_trials << "\n\n";

	cerr << "<v sigma> (unit = pi a^2 v_rel_th):\n\n";
	print_sigma3_array(total_sigma, 1);

	cerr << "<v sigma_err> (unit = pi a^2 v_rel_th):\n\n";
	print_sigma3_err_array(total_err_sq, 1);

	// More detail:

	real sigma_coll = total_coll(total_sigma);

	cerr << "<v sigma> and <v sigma_err> (details):\n\n";
	print_sigma3_nonmergers(total_sigma, total_err_sq, 1);

	if (sigma_coll > 0) 
	    print_sigma3_mergers(total_sigma, total_err_sq, 1);

	if (physical_units_flag) {

	    int i_total = rescale_physical(sma, v_rel_th, v_unit,
					   total_sigma, total_err_sq);
	    real print_factor = pow(10., -i_total); 
	
	    cerr << "<v sigma> (unit = 10^"
		 << i_total << " pc^3/Gyr):\n\n";
	    print_sigma3_array(total_sigma, print_factor);

	    // Use the *same* print_factor for the errors as for sigma:

	    cerr << "<v sigma_err> (unit = 10^"
		 << i_total << " pc^3/Gyr):\n\n";
	    print_sigma3_err_array(total_err_sq, print_factor);

	    Details:

	    cerr << "<v sigma> and <v sigma_err> (details):\n\n";
	    print_sigma3_nonmergers(total_sigma, total_err_sq, print_factor);

	    if (sigma_coll > 0) 
		print_sigma3_mergers(total_sigma, total_err_sq, print_factor);
	}
    }
}

#endif
