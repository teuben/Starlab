
// flyby3.C: Perform a series of three-body scattering experiments,
//	     calculating statistics on changes in eccentricity and
//	     semi-major axis.

// Starlab application:  scatter3.

#include "scatter3.h"

static real e_init;

local real energy(sdyn3* b)
{
    real k = 0, u = 0, dissipation = 0;

    for_all_daughters(sdyn3, b, bi) {
	k += bi->get_mass() * bi->get_vel() * bi->get_vel();
	dissipation += bi->get_energy_dissipation();
	for (sdyn3 * bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister())
	  u -= bi->get_mass() * bj->get_mass()
	       / abs(bi->get_pos() - bj->get_pos());
    }

    return 0.5*k + u + dissipation;
}

// Print hook into low_n3_evolve:

local void print(sdyn3* b)
{
    kepler k;

    // Find the closest pair:

    sdyn3* b1 = NULL;
    sdyn3* b2 = NULL;
    real rmin2 = VERY_LARGE_NUMBER;

    for_all_daughters(sdyn3, b, bi)
	for (sdyn3* bj = bi->get_younger_sister(); bj != NULL;
	     bj = bj->get_younger_sister()) {
	    vec rij = bi->get_pos() - bj->get_pos();
	    real rij2 = rij * rij;
	    if (rij2 < rmin2) {
		b1 = bi;
		b2 = bj;
		rmin2 = rij2;
	    }
	}
    set_kepler_from_sdyn3(k, b1, b2);

    // Identify the third star:

    sdyn3* b3 = NULL;
    for_all_daughters(sdyn3, b, bk)
	if (b3 == NULL && bk != b1 && bk != b2) b3 = bk;
    vec cm = (b1->get_mass()*b1->get_pos() + b2->get_mass()*b2->get_pos())
	          / (b1->get_mass() + b2->get_mass());

    if (e_init == VERY_LARGE_NUMBER) e_init = energy(b);
    cout << "  t = "  << b->get_time() + b->get_time_offset()
	 << "  R = "  << abs(b3->get_pos() - cm)
	 << "  a = "  << k.get_semi_major_axis()
	 << "  e = "  << k.get_eccentricity()
	 << "  dE = " << energy(b) - e_init
	 << endl;
}

#define N_OUT		10	// feedback
#define N_BINS		24
#define BINS_PER_DECADE  2
#define MIN_STATS	 5

local int realcomp(const void * a, const void * b)	// For use by qsort
{
    if (*((real *) a) > *((real *) b))
	return 1;
    else if (*((real *) a) < *((real *) b))
	return -1;
    else
	return 0;
}

local void error_print(ostream& s,
		       int i,
		       intermediate_state3& inter,
		       final_state3& final)

{
    s << "  error at trial #" << i
      << ":  n_kepler = " << inter.n_kepler
      << "  n_osc = " << inter.n_osc
      << "  t = " << final.time
      << "  n_steps = " << final.n_steps
      << ":\n        " << state_string(inter.descriptor)
      << "  "  << state_string(final.descriptor)
      << endl;
}

// run_trials: Perform the experiments and get statistics.

local void run_trials(int n_trials, initial_state3 & init,
		      int out_flag, int planar_flag, int dump_flag,
		      int debug_flag, real dt_print)
{
    intermediate_state3 inter;
    final_state3 final;

    int total_ok = 0;
    real d_ecc_sum = 0, d_ecc_sq_sum = 0, d_ecc_min = 100, d_ecc_max = -1;
    real d_sma_sum = 0, d_sma_sq_sum = 0, d_sma_min = 100, d_sma_max = -1;

    real* d_ecc_list = new real[n_trials];
    real* d_sma_list = new real[n_trials];
    int d_ecc_histo[N_BINS];
    int d_sma_histo[N_BINS];

    real max_error = -VERY_LARGE_NUMBER, min_error = VERY_LARGE_NUMBER;
    real error_sum = 0, abs_error_sum = 0, error_sq_sum = 0;

    int i;
    for (i = 0; i < N_BINS; i++) d_ecc_histo[i] = d_sma_histo[i] = 0;

    init.r_init_min = init.r_init_max = abs(init.r_stop);

    if (debug_flag) cout << endl;

    for (i = 0; i < n_trials; i++) {

	randomize_angles(init.phase);
	if (planar_flag == 1)
	    init.phase.cos_theta = 1;	// Planar prograde
	else if (planar_flag == -1)
	    init.phase.cos_theta = -1;	// Planar retrograde

	e_init = VERY_LARGE_NUMBER;
	scatter3(init, inter, final,
		 VERY_LARGE_NUMBER,
		 VERY_LARGE_NUMBER,
		 VERY_LARGE_NUMBER,
		 0,
		 dt_print, print);

	// On return, the final state should be "stopped" and the
	// intermediate state should be "non_resonance".  Since kepler
	// extension is turned off when r_stop < 0 is specified, the 
	// only possibilities are (for outer separation sep and relative
	// velocity vr):
	//
	//	1. sep > r_stop with vr > 0 during direct integration
	//		==> n_osc = 1, n_kepler = 0,
	//		    final.time is when the criterion was met.
	//	2. apocenter occurred during direct integration
	//		==> n_osc = 1, n_kepler = 0, sep < r_stop, vr < 0,
	//		    final.time is when vr < 0 was first detected,
	//		    (within ~20 time units of apocenter).
	//
	// (Note that the second case guarantees d_sma > 0.)

/*
	cout << "inter:  n_osc = " << inter.n_osc
	     << "  n_kepler = " << inter.n_kepler
	     << "  " << state_string(inter.descriptor) << endl;
	cout << "final:  sep = " << final.outer_separation
	     << "  t = " << final.time << "  n_steps = " << final.n_steps
	     << "  " << state_string(final.descriptor) << endl;
*/

	// Check for errors:

	if (inter.n_kepler > 0
	    || inter.n_osc != 1
	    || final.time <= 0
	    || inter.descriptor != non_resonance
	    || final.descriptor != stopped) {
	    
	    if (debug_flag) 
		error_print(cout, i, inter, final);
	    else
		error_print(cerr, i, inter, final);

	} else {

	    // Gather statistics:

	    real d_ecc = final.ecc - init.ecc;
	    d_ecc_sum += d_ecc;
	    d_ecc_sq_sum += d_ecc*d_ecc;
	    d_ecc_min = min(d_ecc_min, d_ecc);
	    d_ecc_max = max(d_ecc_max, d_ecc);

	    // Eccentricity details:

	    *(d_ecc_list + total_ok) = d_ecc;

	    int ie = (int) (BINS_PER_DECADE * log10(max(1.e-30,abs(d_ecc))))
			+ N_BINS/2 - 1;
	    if (ie < 0)
		ie = 0;
	    else if (ie >= N_BINS/2)
		ie = N_BINS/2 - 1;

	    if (d_ecc > 0) 
		ie += N_BINS/2;
	    else
		ie = N_BINS/2 -1 - ie;

	    d_ecc_histo[ie]++;	    

	    real d_sma = final.sma - init.sma;
	    d_sma_sum += d_sma;
	    d_sma_sq_sum += d_sma*d_sma;
	    d_sma_min = min(d_sma_min, d_sma);
	    d_sma_max = max(d_sma_max, d_sma);

	    // Semi-major axis details:

	    *(d_sma_list + total_ok) = d_sma;

	    int ia = (int) (BINS_PER_DECADE * log10(max(1.e-30,abs(d_sma))))
			+ N_BINS/2 - 1;
	    if (ia < 0)
		ia = 0;
	    else if (ia >= N_BINS/2)
		ia = N_BINS/2 - 1;

	    if (d_sma > 0) 
		ia += N_BINS/2;
	    else
		ia = N_BINS/2 -1 - ia;

	    d_sma_histo[ia]++;	    

	    if (debug_flag) {
		cout << "  trial #" << i
		     << ":  d_ecc, d_sma = " << d_ecc << " " << d_sma
		     << "  bins " << ie << " " << ia << endl;
		cout << "  dr_minima: " << inter.r_min[0]
		     << "  " << inter.r_min[1]
		     << "  " << inter.r_min[2] << endl;
		cout << "  energy";
		if (e_init < VERY_LARGE_NUMBER) cout << " = " << e_init << " ";
		cout << " error = " << final.error << endl;
	    }

	    if (final.error > max_error) max_error = final.error;
	    if (final.error < min_error) min_error = final.error;

	    error_sum += final.error;
	    abs_error_sum += abs(final.error);
	    error_sq_sum += final.error * final.error;

	    total_ok++;
	    if (out_flag && total_ok % N_OUT == 0) cout << total_ok << endl;
	}
    }

    if (total_ok >= MIN_STATS) {

	real mean_error = error_sum / total_ok;
	real mean_abs_error = abs_error_sum / total_ok;
	cout << "\n  total hits = " << total_ok
	     << "  error range = " << min_error << " to " << max_error
	     << "\n                    mean  error  = " << mean_error
	     << "\n                    mean |error| = " << mean_abs_error
	     << "  sigma = "
	     << sqrt(error_sq_sum / total_ok - mean_error * mean_error)
	     << "\n\n";

	real mean_ecc = d_ecc_sum / total_ok;
	cout << "  <d_ecc> = " << mean_ecc << "  sigma = "
	     << sqrt(d_ecc_sq_sum / total_ok - mean_ecc * mean_ecc) << endl;
	cout << "  min d_ecc = " << d_ecc_min
	     << "  max = " << d_ecc_max << endl;

	real mean_sma = d_sma_sum / total_ok;
	cout << "\n  <d_sma> = " << mean_sma << "  sigma = "
	     <<  sqrt(d_sma_sq_sum / total_ok - mean_sma * mean_sma) << endl;
	cout << "  min d_sma = " << d_sma_min
	     << "  max = " << d_sma_max << endl;

	cout << "\nEccentricity details (d_ecc):\n";
	cout.precision(4);

	cout << "  10-percentiles:" << endl << " ";
	qsort((void*)d_ecc_list, (size_t)total_ok, sizeof(real), realcomp);
	real factor = 0.01*total_ok;
	for (i = 10; i <= 90; i += 10)
	    cout << " " << *(d_ecc_list + (int)(factor*i));
	cout << endl;

	cout << "  histogram (" << BINS_PER_DECADE << " bins per decade):\n";
	int j = 0;
	while (j < N_BINS) {	// cout and printf don't mesh on Solaris 2.3!
	    printf("\t");
	    for (i = 0; i < BINS_PER_DECADE && j < N_BINS; i++, j++)
		printf("  %5d", d_ecc_histo[j]);	
	    printf("\n");
	}

	cout << "\nSemi-major axis details (d_sma):\n";
	cout.precision(4);

	cout << "  10-percentiles:" << endl << " ";
	qsort((void*)d_sma_list, (size_t)total_ok, sizeof(real), realcomp);
	for (i = 10; i <= 90; i += 10)
	    cout << " " << *(d_sma_list + (int)(factor*i));
	cout << endl;

	cout << "  histogram (" << BINS_PER_DECADE << " bins per decade):\n";
	j = 0;
	while (j < N_BINS) {
	    printf("\t");
	    for (i = 0; i < BINS_PER_DECADE && j < N_BINS; i++, j++)
		printf("  %5d", d_sma_histo[j]);	
	    printf("\n");
	}

	if (dump_flag) {
	    cout << "\n  Raw data (total = " << total_ok << ", sorted):\n";

	    j = 0;
	    cout << "\n  Eccentricity:\t";
	    while (j < total_ok) {
		for (i = 0; i < 5 && j < total_ok; i++, j++)
		    printf("  %12.8f", *(d_ecc_list + j));
		cout << endl;
	    }

	    j = 0;
	    cout << "\n  Semi-major axis:\t";
	    while (j < total_ok) {
		for (i = 0; i < 5 && j < total_ok; i++, j++)
		    printf("  %12.8f", *(d_sma_list + j));
		cout << endl;
	    }
	}
    }
}

main(int argc, char **argv)
{
    int  n_trials = 100;
    int  seed = 0;
    int  out_flag = 0;
    int  dump_flag = 0;
    int  debug_flag = 0;
    int  planar_flag = 0;
    real dt_print = VERY_LARGE_NUMBER;

    extern char *poptarg;
    int c;
    const char *param_string = "A:dDe:m:M:n:opPr:R:s:t:v:";

    initial_state3 init;
    make_standard_init(init);

    init.rho = 5;
    init.r_stop = 100;
    init.v_inf = 0;

    // Defaults: n_trials = 100
    //		 ecc   = 0
    //		 rho   = 5
    //		 phase = random
    //		 v_inf = 0 (parabola, rho = periastron)
    //		 m1 = m2 = m3 = 0.5
    //           r_stop = 100
    //           no extra output

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'd': dump_flag = 1 - dump_flag;
		      break;
	    case 'D': debug_flag = 1 - debug_flag;
		      break;
	    case 'e': init.ecc = atof(poptarg);
		      break;
	    case 'm': init.m2 = atof(poptarg);
		      break;
	    case 'M': init.m3 = atof(poptarg);
		      break;
	    case 'n': n_trials = atoi(poptarg);
		      break;
	    case 'o': out_flag = 1 - out_flag;
		      break;
	    case 'p': planar_flag = 1;
		      break;
	    case 'P': planar_flag = -1;
		      break;
	    case 'r': init.rho = atof(poptarg);
		      break;
	    case 'R': init.r_stop = -atof(poptarg); // Kludge: stop at apocenter
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 't': dt_print = atof(poptarg);
		      debug_flag = 1;
		      break;
	    case 'v': init.v_inf = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	    }

    srandinter(seed);
    real peri = init.rho;

    // Force input rho to mean periastron in ALL cases...
    // Remember to scale the velocity by v_crit.

    if (init.v_inf > 0)
	init.rho = sqrt(init.rho * (init.rho + 2 * init.m3
				                 / ((1 - init.m2) * init.m2
				                       * pow(init.v_inf, 2))));

    // Basic output on the experiment:

    cout << "  Masses:"
         << "  m1 = " << 1 - init.m2
	 << "  m2 = " << init.m2
	 << "  m3 = " << init.m3
	 << endl;
    cout << "  v_inf = " << init.v_inf
	 << "  periastron = " << peri
	 << "  rho = " << init.rho
	 << "  ecc = " << init.ecc;
    if (planar_flag == 1)
	cout << "  planar prograde orbit";
    else if (planar_flag == -1)
	cout << "  planar retrograde orbit";
    cout << endl;

    cout << "  Start/stop at R > " << abs(init.r_stop)
	 << "  accuracy parameter = " << init.eta
	 << "  random seed = " << get_initial_seed()
	 << endl;

    cpu_init();
    run_trials(n_trials, init, out_flag, planar_flag, 
	       dump_flag, debug_flag, dt_print);
}
