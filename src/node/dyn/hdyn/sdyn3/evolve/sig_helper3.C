
// sig_helper3.C: Engine for determining cross-sections in
//                three-body scattering.

// Differences from the previous version:
//
//	more streamlined scan over impact parameter
//	decoupling of initialization, integration, and analysis functions
//	parallel operation, eventually

// Starlab library function.

#include "sigma3.h"
#include "pvm_scatt.h"

// Pretty-print a scatter profile:

void print_profile(ostream & s, scatter_profile & p, int prec)
{
    int pp = s.precision(prec);

    s << "scatter_profile:" << endl
      << "   m1 = " << 1 - p.m2 << "  m2 = " << p.m2 << "  m3 = " << p.m3
      << "  r1 = " << p.r1 << "  r2 = " << p.r2 << "  r3 = " << p.r3;

    if (p.v_inf >= 0) s << "  v_inf = " << p.v_inf;

    s << endl;

    if (p.rho_flag)
	s << "   rho_min = " << p.rho_min
	  << " rho_max = " << p.rho_max << endl;

    if (p.ecc_flag)
	s << "   binary_ecc = " << p.ecc << endl;

    if (p.phase_flag.cos_theta || p.phase_flag.phi || p.phase_flag.psi
	 || p.phase_flag.mean_anomaly) {

        s << "  phase = ";

	if (p.phase_flag.cos_theta)
	    s << acos(p.phase.cos_theta) << " ";
	else
	    s << "- ";

	if (p.phase_flag.phi)
	    s << p.phase.phi << " ";
	else
	    s << "- ";

	if (p.phase_flag.psi)
	    s << p.phase.psi << " ";
	else
	    s << "- ";

	if (p.phase_flag.mean_anomaly)
	    s << p.phase.mean_anomaly;
	else
	    s << "-";
	
	s << endl;
	}

    s.precision(pp);
}

// Set up a template profile:

void make_standard_profile(scatter_profile & prof)
{
    prof.m2 = 0.5;
    prof.m3 = 0.5;
    prof.r1 = 0;
    prof.r2 = 0;
    prof.r3 = 0;
    prof.v_inf = 1;
    prof.tidal_tol_factor = DEFAULT_TIDAL_TOL_FACTOR;
    prof.rho_flag = 0;					// random impact par.
    prof.ecc_flag = 0;					// random eccentricity
    prof.phase_flag.cos_theta = 0;			// random angles
    prof.phase_flag.phi = 0;
    prof.phase_flag.psi = 0;
    prof.phase_flag.mean_anomaly = 0;
    prof.ecc = 0;			                // for definiteness
    prof.rho_min = 0;
    prof.rho_max = 0;
    prof.phase.cos_theta = 0;
    prof.phase.phi = 0;
    prof.phase.psi = 0;
    prof.phase.mean_anomaly = 0;
    prof.eta = DEFAULT_ETA;
    prof.intermediate_target = unknown_intermediate;
    prof.final_target1 = unknown_final;
    prof.final_target2 = unknown_final;
}

// Turn a scatter profile into a default initial state:

void prof_to_init(scatter_profile & prof, initial_state3 & init)
{
    if (prof.m2 > 1) {
	cerr << "prof_to_init: prof.m2 = " << prof.m2 << " > 1" << endl;
	exit(1);
    }

    make_standard_init(init);

    init.m2 = prof.m2;
    init.m3 = prof.m3;
    init.r1 = prof.r1;
    init.r2 = prof.r2;
    init.r3 = prof.r3;
    init.sma  = 1;

    // Maximum eccentricity to stay well away (by factor ECC_TOL) from contact:

    real ecc_max = 1 - ECC_TOL * (init.r1 + init.r2) / init.sma;
    if (ecc_max < 0)
	err_exit("prof_to_init: radii too large to initialize binary");

    init.v_inf = prof.v_inf;

    if (prof.rho_flag)
      init.rho = sqrt(randinter(prof.rho_min*prof.rho_min,
				prof.rho_max*prof.rho_max));
    else
      init.rho = 0;                       // The impact parameter rho will be
                                          // initialized elsewhere if random
    init.tidal_tol_factor = prof.tidal_tol_factor;
    init.eta = prof.eta;

    randomize_angles(init.phase);
    if (prof.phase_flag.mean_anomaly)
	init.phase.mean_anomaly = prof.phase.mean_anomaly;
    if (prof.phase_flag.cos_theta)
	init.phase.cos_theta = prof.phase.cos_theta;
    if (prof.phase_flag.phi) init.phase.phi = prof.phase.phi;
    if (prof.phase_flag.psi) init.phase.psi = prof.phase.psi;

    //  Initialize the eccentricity last for ease of programming in sigma3...

    if (prof.ecc_flag == 0) {

	// Thermal distribution of eccentricities in the allowed range.

	init.ecc = sqrt(randinter(0, ecc_max*ecc_max));

    } else if (prof.ecc_flag == 1) {

	// A particular eccentricity was specified.

	init.ecc = prof.ecc;

    } else if (prof.ecc_flag == 2) {

	// Bizarre observed eccentricity distribution (Mathieu)...

	init.ecc = -1;
	while (init.ecc < 0 || init.ecc > 1) init.ecc = gausrand(0.3, 0.15);

    } else

	err_exit("prof_to_init: illegal eccentricity flag");

}

int n_coll(int n_hits[N_INTER][N_FINAL][N_RHO_ZONE_MAX], int i_max)
{
    int total = 0;	

    for (int i_int = 0; i_int < N_INTER; i_int++)
	for (int i_rho = 0; i_rho <= i_max; i_rho++)
	    total += n_hits[i_int][merger_binary_1][i_rho]
	                  + n_hits[i_int][merger_binary_2][i_rho]
	                  + n_hits[i_int][merger_binary_3][i_rho]
			  + n_hits[i_int][merger_escape_1][i_rho]
	                  + n_hits[i_int][merger_escape_2][i_rho]
	                  + n_hits[i_int][merger_escape_3][i_rho]
			  + n_hits[i_int][triple_merger][i_rho];
    return total;
}

#define output_header_1 \
  "   pres   exch    ion merg_b merg_e merg_3  error   stop unknown\n"
//     0      1       2     3      4      5      6	<-- n/s_total
#define output_header_2 \
  "   pres   exch   ion  merg_b merg_e merg_3  error   stop\n"

#define output_header_3 \
  "   pres  exch1  exch2    ion  error   stop\n"

#define output_header_4 \
  "  mrg_b1 mrg_b2 mrg_b3 mrg_e1 mrg_e2 mrg_e3 3coll\n"

static char *f_lab[] = {"non_res", "hier_res", "dem_res", "unknown"};
// (static here to keep Sun CC happy...)

// NOTE: All "print_sigma" functions start on the current line (i.e. they
//       do not skip a line at the start) and skip a line at the end.

char dummy_string[20];

void print_sigma3_counts(int n_hits[N_INTER][N_FINAL][N_RHO_ZONE_MAX],
			 int i_max)
{
    cerr << output_header_1;

    for (int i_int = 0; i_int < N_INTER; i_int++) {

	int n_total[N_FINAL - 5];
	for (int k = 0; k < N_FINAL - 5; k++) n_total[k] = 0;

	// Combine some categories:

	for (int i_rho = 0; i_rho <= i_max; i_rho++) {
	    n_total[0] += n_hits[i_int][preservation][i_rho];
	    n_total[1] += n_hits[i_int][exchange_1][i_rho]
	                  + n_hits[i_int][exchange_2][i_rho];
	    n_total[2] += n_hits[i_int][ionization][i_rho];
	    n_total[3] += n_hits[i_int][merger_binary_1][i_rho]
	                  + n_hits[i_int][merger_binary_2][i_rho]
	                  + n_hits[i_int][merger_binary_3][i_rho];
	    n_total[4] += n_hits[i_int][merger_escape_1][i_rho]
	                  + n_hits[i_int][merger_escape_2][i_rho]
	                  + n_hits[i_int][merger_escape_3][i_rho];
	    n_total[5] += n_hits[i_int][triple_merger][i_rho];
	    n_total[6] += n_hits[i_int][error][i_rho];
	    n_total[7] += n_hits[i_int][stopped][i_rho];
	    n_total[8] += n_hits[i_int][unknown_final][i_rho];
	}

	for (int i_fin = 0; i_fin < N_FINAL - 5; i_fin++) {

	    cerr << " ";

	    sprintf(dummy_string, "%d", n_total[i_fin]);

	    // Use exactly 6 chars (+ 1 space) unless integer is very large.

	    for (int i = 0; i < Starlab::max(0, 6 - (int)strlen(dummy_string)); i++)
		cerr << " ";
	    
	    cerr << dummy_string;
	}
	cerr << "  " << f_lab[i_int] << endl;
    }
    cerr << endl;
}

local void print_formatted(real x)
{
    cerr << " ";

    if ( (x < 1e6 && x > 0.01) || x <= 0) {

	sprintf(dummy_string, "%.3f", x);

	// Use exactly 7 characters unless integer is very large.

	for (int i = 0; i < Starlab::max(1, 6 - (int)strlen(dummy_string)); i++)
	    cerr << " ";

	// Truncate the string, if x large (max = 999999).

	if (strlen(dummy_string) > 6) dummy_string[6] = '\0';

	cerr << dummy_string;

    } else if (x < 0.01)
	fprintf(stderr, "%7.1e", x);	// 2 sig. fig. for small numbers
    else
	fprintf(stderr, "%8.2e", x);	// 3 sig. fig. for large numbers
}

local void print_error(real x)
{
    cerr << " ";

    if (x < 1) {
	if (x > 0.01 || x <= 0) {
	    sprintf(dummy_string, "%.3f", x);
	    fprintf(stderr, "+-%s", dummy_string+1);
	} else
	    fprintf(stderr, "+-%5.0e", x);
    } else if (x < 10)
	fprintf(stderr, "+-%.2f", x);
    else
	fprintf(stderr, "+-%6.1e", x);
}

#define N_SIGMA3_ARR (N_FINAL - 6)

void print_sigma3_array(real sigma[N_INTER][N_FINAL], real scale_factor,
			int sqrt_flag)
{
    cerr << output_header_2;

    for (int i_int = 0; i_int < N_INTER - 1; i_int++) {

	real s_total[N_SIGMA3_ARR];

	// Combine some categories (as in print_sigma3_counts):

	s_total[0] = sigma[i_int][preservation];
	s_total[1] = sigma[i_int][exchange_1]
	              + sigma[i_int][exchange_2];
	s_total[2] = sigma[i_int][ionization];
	s_total[3] = sigma[i_int][merger_binary_1]
	              + sigma[i_int][merger_binary_2]
	              + sigma[i_int][merger_binary_3];
	s_total[4] = sigma[i_int][merger_escape_1]
	              + sigma[i_int][merger_escape_2]
	              + sigma[i_int][merger_escape_3];
	s_total[5] = sigma[i_int][triple_merger];
	s_total[6] = sigma[i_int][error];
	s_total[7] = sigma[i_int][stopped];

	for (int i_fin = 0; i_fin < N_SIGMA3_ARR; i_fin++) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "       ");
	    else {
		real temp = (sqrt_flag == 0 ? s_total[i_fin]
			                    : sqrt(s_total[i_fin]));
		print_formatted(scale_factor * temp);
	    }
	}
	cerr << "         " << f_lab[i_int] << endl;
    }
    cerr << endl;
}

// print_sigma3_err_array: print errors, assuming that SQUARES are sent.

void print_sigma3_err_array(real sigma_err_sq[N_INTER][N_FINAL],
			    real scale_factor)
{
    print_sigma3_array(sigma_err_sq, scale_factor, 1);
}

#define N_SIGMA3_NM_ARR (N_FINAL - 1)

void print_sigma3_nonmergers(real sigma[N_INTER][N_FINAL],
			     real sigma_err_sq[N_INTER][N_FINAL],
			     real v2)
{
    cerr << output_header_3 << endl; 

    for (int i_int = 0; i_int < N_INTER - 1; i_int++) {

	// Cross-sections:

	for (int i_fin = 0; i_fin <  N_SIGMA3_NM_ARR; i_fin++) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "       ");
	    else if (i_fin <= 3 || i_fin >=  N_SIGMA3_NM_ARR - 2)
		print_formatted(v2 * sigma[i_int][i_fin]);
	}
	cerr << "                       " << f_lab[i_int] << endl;

	// Errors:

	for (int i_fin = 0; i_fin <  N_SIGMA3_NM_ARR; i_fin++) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "       ");
	    else if (i_fin <= 3 || i_fin >=  N_SIGMA3_NM_ARR - 2)
		print_error(v2 * sqrt(sigma_err_sq[i_int][i_fin]));
	}
	cerr << "\n\n";
    }
}

void print_sigma3_mergers(real sigma[N_INTER][N_FINAL],
			  real sigma_err_sq[N_INTER][N_FINAL],
			  real v2)
{
    cerr << output_header_4 << endl;

    for (int i_int = 0; i_int < N_INTER - 1; i_int++) {

	// Cross-sections:
	    
	for (int i_fin = 4; i_fin < 11; i_fin++) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "       ");
	    else
		print_formatted(v2 * sigma[i_int][i_fin]);
	}
	cerr << "                " << f_lab[i_int] << endl;
	    
	// Errors:
	    
	for (int i_fin = 4; i_fin < 11; i_fin++) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "       ");
	    else
		print_error(v2 * sqrt(sigma_err_sq[i_int][i_fin]));
	}
	cerr << "\n\n";
    }
}

local void create_temp(char* temp, int length, char* string)
{
    for (int k = 0; k < length; k++) temp[k] = ' ';
    strcpy(temp, string);
    for (int k = 0; k < length; k++) if (temp[k] == '\0') temp[k] = ' ';
    temp[length-1] = '\0';
}

void print_sigma3(sigma_out & out, real v2)
{
    int i;

    cerr << endl;
    for (i = 0; i < 72; i++) cerr << "-";
    cerr << endl;

    int p = cerr.precision(STD_PRECISION);
    cerr << "\n   rho_max = " << out.rho_max
	 << "  i_max = " << out.i_max << endl;
    cerr << "   central_trial_density = " << out.central_trial_density
	 << "  total_trials = " << out.total_trials << endl;

    cerr << "   average steps = " << out.total_steps / (real) out.total_trials
	 << "  maximum = " << out.max_steps << endl;
    cerr << "   timesteps    (bin factor = 10): ";
    for (i = 0; i < N_STEP_BIN; i++) cerr << " " << out.step_counter[i];
    cerr << endl;

    int coll_total = n_coll(out.n_hits, out.i_max);

    cerr << "   oscillations (bin factor =  2):" << endl;
    int jmax = (int)triple_merger;
    if (coll_total == 0) jmax = (int)ionization;

#   define N_TEMP   36
#   define N_INDENT 17

    char temp[N_TEMP];
    for (int k = 0; k < N_TEMP; k++) temp[k] = ' ';

    for (int j = 0; j <= jmax; j++) {
	create_temp(temp + N_INDENT, N_TEMP - N_INDENT,
		    state_string((final_descriptor3)j));
	cerr << temp;
        for (i = 0; i < N_OSC_BIN; i++) cerr << " " << out.osc_counter[i][j];
	cerr << endl;
    }

    create_temp(temp + N_INDENT, N_TEMP - N_INDENT, "total");
    cerr << temp;

    int osc_total[N_OSC_BIN];
    for (i = 0; i < N_OSC_BIN; i++) osc_total[i] = 0;
    for (int j = 0; j <= jmax; j++)
        for (i = 0; i < N_OSC_BIN; i++) osc_total[i] += out.osc_counter[i][j];
    for (i = 0; i < N_OSC_BIN; i++) cerr << " " << osc_total[i];
    cerr << endl;

    // Raw counts:

    cerr << endl << "raw counts:\n\n";
    print_sigma3_counts(out.n_hits, out.i_max);

    // Cross-sections:

    cerr << "v^2 sigma / (pi a^2):\n\n";
    print_sigma3_array(out.sigma, v2);

    // More detail and errors:

    // Non-mergers first (0-3, N_FINAL-2):

    print_sigma3_nonmergers(out.sigma, out.sigma_err_sq, v2);

    // Then mergers (4-10):

    if (coll_total > 0) 
	print_sigma3_mergers(out.sigma, out.sigma_err_sq, v2);

    cerr << flush;
    cerr.precision(p);
}

// print_counts: print out enough information to determine cross-sections
//		 and errors, and to combine them with other runs.

void print_all_sigma3_counts(sigma_out & out, ostream & s)
{
    int p = s.precision(STD_PRECISION);
    s << "Counts:"
	 << "\n  i_max = " << out.i_max
	 << "  trials_per_zone = " << out.trials_per_zone
	 << "\n  rho_sq_init = " << out.rho_sq_init
	 << "  RHO_SQ_FACTOR = " << RHO_SQ_FACTOR << endl;

    for (int i_int = 0; i_int < N_INTER; i_int++)
	for (int i_fin = 0; i_fin < N_FINAL; i_fin++)
	    for (int i_rho = 0;	i_rho <= out.i_max; i_rho++)
		s << out.n_hits[i_int][i_fin][i_rho] << " ";
    s << endl;
    cerr.precision(p);
} 

//------------------------------------------------------------------------

// Encapsulate the following functions -- efficiency is not an issue here.

real zone_area(sigma_out& out, int i)	// (without the PI)
{
    real rho_sq_min = 0, rho_sq_max;

    rho_sq_max = out.rho_sq_init * pow((real)RHO_SQ_FACTOR, i);
    if (i > 0) rho_sq_min = rho_sq_max / RHO_SQ_FACTOR;

    return rho_sq_max - rho_sq_min;
}

real zone_density(sigma_out& out, int i)
{
    return out.trials_per_zone / zone_area(out, i);
}

real zone_weight(sigma_out& out, int i)
{
    return zone_area(out, i) / out.trials_per_zone;
}

//------------------------------------------------------------------------

// counts_to_sigma: Convert raw counts into cross-sections and errors.
//		    This will update sigma from the current counts,
//		    making NO other changes in the sigma_out structure.

void counts_to_sigma(sigma_out & out)
{
    out.sigma_total = 0;
    out.sigma_total_err_sq = 0;

    for (int i_int = 0; i_int < N_INTER; i_int++)
	for (int i_fin = 0; i_fin < N_FINAL; i_fin++) {

	    out.sigma[i_int][i_fin] = 0;
	    out.sigma_err_sq[i_int][i_fin] = 0;

	    for (int i_rho = 0; i_rho <= out.i_max; i_rho++) {

//		if (out.n_hits[i_int][i_fin][i_rho] > 0)
//		    cerr<<"i_int = "<<i_int<<" i_fin = "<<i_fin
//			<<" i_rho = "<<i_rho<<" "
//			<<out.n_hits[i_int][i_fin][i_rho]<<endl;

		real weight = zone_weight(out, i_rho);

		real delta_sigma =  weight * out.n_hits[i_int][i_fin][i_rho];
		real delta_err_sq = weight * weight
		                     * out.n_hits[i_int][i_fin][i_rho];

		out.sigma[i_int][i_fin] += delta_sigma;
		out.sigma_err_sq[i_int][i_fin] += delta_err_sq;

		if (i_int > 0 || i_fin > 0) {
		    out.sigma_total += delta_sigma;
		    out.sigma_total_err_sq += delta_err_sq;
		}
	    }
	}
} 

// Elementary statistics:

real max_t = 0;					       // global variables!!
real max_err = 0;

local void print_status(real v_inf, sigma_out & out, int debug)
{
    int p = cerr.precision(STD_PRECISION);
    cerr.setf(ios::fixed, ios::floatfield);
     
    counts_to_sigma(out);
 
    cerr << "sigma[" << v_inf << "] = " << out.sigma_total
 	 << " +/- " << sqrt(out.sigma_total_err_sq)
	 << "  v^2 sigma = " << out.sigma_total*v_inf*v_inf
 	 << " +/- " << sqrt(out.sigma_total_err_sq)*v_inf*v_inf
	 << endl;
 
    if (debug < 0)
	cerr << "  n_total = "	<< out.total_trials
	     << "  max_t = "	<< max_t
	     << "  max_err = "	<< max_err
	     << "  n_hit_tot = "<< out.n_hit_tot
	     << endl;

    cerr.precision(p);
}

// summarize_scattering_initial: Print a line specifying a
//				 scattering experiment.

void summarize_scattering_initial(initial_state3 & init,
				  int n_rand,
				  real dt_snap,
				  real snap_cube_size)
{
    int p = cerr.precision(STD_PRECISION);

    cerr << "\nscatter3"
	 << " -s " << get_initial_seed()
	 << " -N " << n_rand;

    cerr.precision(16);				// Uglify the output
    cerr << " -m " << init.m2
	 << " -M " << init.m3
	 << " -v " << init.v_inf
	 << " -e " << init.ecc
	 << " -r " << init.rho
         << " -x " << init.r1
	 << " -y " << init.r2
	 << " -z " << init.r3
	 << " -A " << init.eta
	 << " -g " << init.tidal_tol_factor;

    cerr.precision(STD_PRECISION);
    if (dt_snap < VERY_LARGE_NUMBER) {
	cerr << " -D " << dt_snap
	     << " -C " << snap_cube_size;
    }

    cerr << endl;

    // print_initial(cerr, init);

    cerr.precision(p);
}

// summarize_scattering_final: Print a line specifying a scattering result.

void summarize_scattering_final(intermediate_state3 & inter,
				final_state3 & final,
				int level,
				real cpu)
{
    cerr << "outcome:  ";
    print_scatter3_outcome(inter, final, cerr);

    if (level > 1) print_scatter3_summary(inter, final, cpu, cerr);

//    print_final(cerr, final);
}

//=========================================================================
//
// Functions to initialize, perform, and analyze a scattering experiment.
//
//=========================================================================

// single_scatter_init: Initialize a scattering experiment in the specified
//                      rho range with the given v_inf.

void single_scatter_init(scatter_profile & prof,
			 real rho_sq_min, real rho_sq_max,
			 initial_state3 & init, int & n_rand,
			 int scatter_summary, real dt_snap, real snap_cube_size)
{
    n_rand = get_n_rand();

    // Initialize the scattering experiment:

    prof_to_init(prof, init);

    // Choose rho randomly from the allowed range:

    if (prof.rho_flag == 0)
	init.rho = sqrt(randinter(rho_sq_min, rho_sq_max));
    else
	init.rho = rho_sq_min;

    if (scatter_summary > 0)
	summarize_scattering_initial(init, n_rand, dt_snap, snap_cube_size);

}

// single_scatter: Perform a scattering experiment with the specified init.

int single_scatter(initial_state3 & init,
		   intermediate_state3 & inter,
		   final_state3 & final,
		   real cpu_time_check,
		   real dt_snap,
		   real snap_cube_size)
{
    // Perform a scattering experiment with the specified init, returning
    // result = false if the outcome is a non-resonant preservation.

    scatter3(init, inter, final, cpu_time_check,
	     VERY_LARGE_NUMBER, dt_snap, snap_cube_size);

    // Return a "hit" if the result was
    //
    //		(1) not a flyby (preservation non_resonance)
    //		(2) not an error, stopped, or unknown_final.
    //
    // Modify this return statement to change the definition of a "hit".
    // Probably should take into account the size of the perturbation?

    return (final.descriptor == error
	    || final.descriptor == stopped
	    || final.descriptor == unknown_final
	    || (final.descriptor == preservation
		&& inter.descriptor == non_resonance) ? 0 : 1);
}

// single_scatter_stats: Accumulate information on scattering outcomes.

void single_scatter_stats(scatter_profile & prof,
			  initial_state3 & init,
			  intermediate_state3 & inter,
			  final_state3 & final,
			  int rho_bin_index,
			  sigma_out & out,
			  stat_fp acc_stats,
			  int result)
{
    // Accumulate statistics:

    max_t = Starlab::max(max_t, final.time);
    max_err = Starlab::max(max_err, abs(final.error));

    // Accumulate results on data not available to higher levels:

    out.n_hits[inter.descriptor][final.descriptor][rho_bin_index]++;

    out.total_steps += final.n_steps;
    out.max_steps = Starlab::max(out.max_steps, final.n_steps);
    real logn = log10((real)Starlab::max(1, final.n_steps));
    int index = Starlab::min(N_STEP_BIN-1, (int)logn);
    out.step_counter[index]++;

    logn = log10((real)Starlab::max(1, inter.n_osc)) / log10(2.0);
    index = Starlab::min(N_OSC_BIN-1, (int)logn);
    out.osc_counter[index][final.descriptor]++;

    // Convenient to update these global counters here also:

    out.total_trials++;
    out.n_hit[rho_bin_index] += result;
    out.n_hit_tot += result;

    //------------------------------------------------------------------
    // Optional call to user-supplied function for more detailed
    // statistics gathering. The function must maintain its own internal
    // data in a way that can be extracted later.
    //------------------------------------------------------------------

    if (acc_stats)
	(*acc_stats)(prof, init, inter, final, rho_bin_index, out);
}

//=========================================================================

// get_sigma3: Run three-body scatterings for a given velocity at infinity
//             and return cross-sections for all possible outcomes.
//             Stop when the specified trial density is exceeded.
//
//	       The trials are performed by function "trial" and debugging
//             information is printed out by the function "print".
//
// In addition to determining cross-sections, get_sigma3 may also be used
// as a convenient means of driving a more specialized function acc_stats,
// which can gather additional statistics.  Get_sigma3 guarantees that the
// data sent to acc_stats are properly sampled over impact parameter space.

void get_sigma3(real max_trial_density,	// "density" of trials
		scatter_profile & prof,	// profile describing the run
		sigma_out & out,	// outcome structure
		int debug,		// output control
		real cpu_time_check,    // interval to print CPU time
		real dt_snap,		// snapshot output interval
		real snap_cube_size,    // size of snap cube
		int scatter_summary_flag,// output on each scattering
		stat_fp acc_stats,	// optional function to accumulate
					//     statistics on the runs
		print_fp print)		// function to print system
					//     status
{
    if (prof.rho_flag)			// Impact parameter has been specified.
	err_exit("get_sigma3: prof.rho_flag has been set!");

    // (This routine is intended for calculation of the cross-section
    // without restrictions on the impact parameter)

    if (prof.v_inf <= 0)
	err_exit("get_sigma3: sigma undefined for parabolic orbits");

    if (prof.r1 + prof.r2 >= 1)
	err_exit("get_sigma3: contact binary");

    // Initialize counters.

    int i_int, i_fin, rho_zone;

    for (i_int = 0; i_int < N_INTER; i_int++)
	for (i_fin = 0; i_fin < N_FINAL; i_fin++)
	    for (rho_zone = 0; rho_zone < N_RHO_ZONE_MAX; rho_zone++)
	    out.n_hits[i_int][i_fin][rho_zone] = 0;

    out.total_trials = 0;
    out.n_hit_tot = 0;
    for (rho_zone = 0; rho_zone < N_RHO_ZONE_MAX; rho_zone++)
	out.n_hit[rho_zone] = 0;
    
    out.max_steps = 0;
    out.total_steps = 0;
    for (int i = 0; i < N_STEP_BIN; i++) out.step_counter[i] = 0;

    for (int i = 0; i < N_OSC_BIN;  i++)
        for (int j = 0; j < N_FINAL; j++)
	    out.osc_counter[i][j] = 0;

    real cpu_init = cpu_time(); // The CPU time when we entered get_sigma3
    real cpu_save = cpu_init;

    int scatt_total = 0;	// Maintain these for use in parallel case
    real cpu_total = 0;

    out.i_max = 1;  // i value of safety bin, to be sampled, but supposed
                    // to contain no events, i.e. rho bins run from 0 to 
                    // i_max, with i_max empty

    // Take gravitational focusing into account in the allowed range
    // of rho squared.

    real v_sq = prof.v_inf * prof.v_inf;
    out.rho_sq_init = RHO_SQ_SCALE * (1 + 2 * (1 + prof.m3) / v_sq);

    // We normalize the density of trials by requiring a given number of
    // trials in an effective  pericenter target area, i.e. a typical area
    // of interest, measured in "pericenter" units (rho units with 
    // gravitational focusing scaling taken out) within which interesting
    // results are expected.
    //
    // For hard binaries (v_inf --> 0), we just take the geometrical area of 
    // the binary, since any encounter within this pericenter area is likely
    // to do something interesting.
    //
    // For soft binaries (v_inf --> infinity), we must take a smaller area,
    // by a factor 1/v_sq, since in the impulse approximation the same momentum
    // transfer at higher energy requires a closer encounter by a factor of
    // 1/v_inf.  Thus, the effective density is reduced by a factor 1/v_sq.

    real pericenter_target_factor = 1 / (1 + v_sq);

    out.trials_per_zone = (int) (max_trial_density * RHO_SQ_SCALE
				  / pericenter_target_factor + 0.5);
    out.central_trial_density = (out.trials_per_zone / RHO_SQ_SCALE)
					* pericenter_target_factor;

    // Some useful diagnostic feedback (before the fact):

    if (debug != 0) cerr << "==========\n";
    else cerr << endl;

    cerr << "get_sigma3:  v_inf = " << prof.v_inf
	 << ",  trials_per_zone = " << out.trials_per_zone << endl;

    real rho_sq_min, rho_sq_max;

    initialize_processors(out.trials_per_zone, debug);	// PVM or dummy call

    for (rho_zone = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
	 rho_zone <= out.i_max;
	 rho_zone++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

	// Perform out.trials_per_zone scatterings in each annulus,
	// and keep track of the results.

	int result = multiscatter3(prof, out,
				   rho_sq_min, rho_sq_max, rho_zone,
				   dt_snap, snap_cube_size,
				   cpu_time_check, cpu_init, cpu_save,
				   scatt_total, cpu_total, acc_stats,
				   debug, scatter_summary_flag);

	if (debug)
	    cerr << "rho_zone = " << rho_zone << ",  result = " << result;

	if (rho_zone == out.i_max && result > 0) out.i_max++;

	if (debug)
	    cerr << ",  out.i_max = " << out.i_max << endl;

    }

    cerr << "\ntotal scatterings = " << scatt_total
	 << ",  CPU time = " << cpu_total << endl;

    terminate_processors(debug);	// PVM or dummy call

    if (abs(debug) == 1) (*print)(prof.v_inf, out, debug);

    out.rho_max = sqrt(rho_sq_max);
    counts_to_sigma(out);
}

//=========================================================================

// Establish versions of get_sigma3 with specific defaults, using the
// print_status routine provided in this file.

// NOTE: This version fixes the use of the old "single_scatter" as
//       the basic scattering function.

// Add DEFAULT print_status and NULL acc_stats:

void get_sigma3(real max_trial_density,
		scatter_profile & prof,
		sigma_out & out,
		int debug,
		real cpu_time_check,
		real dt_snap,
		real snap_cube_size,
		int scatter_summary_flag)
{
    get_sigma3(max_trial_density, prof, out, debug, cpu_time_check,
	       dt_snap, snap_cube_size, scatter_summary_flag,
	       (stat_fp) NULL,
	       print_status);
}

// Add DEFAULT print_status to SPECIFIED acc_stats:

void get_sigma3(real max_trial_density,
		scatter_profile & prof,
		sigma_out & out,
		int debug,
		real cpu_time_check,
		real dt_snap,
		real snap_cube_size,
		int scatter_summary_flag,
		stat_fp acc_stats)
{
    get_sigma3(max_trial_density, prof, out, debug, cpu_time_check,
	       dt_snap, snap_cube_size, scatter_summary_flag,
	       acc_stats,
	       print_status);
}


// Add NULL acc_stats to SPECIFIED print_status:

void get_sigma3(real max_trial_density,
		scatter_profile & prof,
		sigma_out & out,
		int debug,
		real cpu_time_check,
		real dt_snap,
		real snap_cube_size,
		int scatter_summary_flag,
		print_fp print)
{
    get_sigma3(max_trial_density, prof, out, debug, cpu_time_check,
	       dt_snap, snap_cube_size, scatter_summary_flag,
	       (stat_fp) NULL,
	       print);
}
