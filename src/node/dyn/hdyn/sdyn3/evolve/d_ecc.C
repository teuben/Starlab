
// d_ecc.C: Perform a series of three-body scattering experiments,
//	  calculating statistics on the eccentricity change.

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
    cout << b->get_time() + b->get_time_offset() << " "
	 << abs(b3->get_pos() - cm) << " "
	 << k.get_eccentricity() << " "
	 << energy(b) - e_init << endl;
}

#define N_OUT		10	// feedback
#define N_BINS		24
#define BINS_PER_DECADE  4
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

// run_trials: Perform the experiments and get statistics.

local void run_trials(int n_trials, initial_state3 & init,
		      int out_flag, int planar_flag, int dump_flag,
		      int debug_flag, real dt_print)
{
    intermediate_state3 inter;
    final_state3 final;

    int total_ok = 0;
    real d_ecc_sum = 0, d_ecc_sq_sum = 0, d_ecc_min = 100, d_ecc_max = -1;

    real* de = new real[n_trials];
    int de_histo[N_BINS];

    int i;
    for (i = 0; i < N_BINS; i++) de_histo[i] = 0;

    init.r_init_min = init.r_init_max = init.r_stop;

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

	// Check for return passages:

	if (inter.n_kepler > 0 || inter.n_osc != 1
	     || inter.descriptor != non_resonance)
	    cerr << "  error at trial #" << i
		 << ":  n_kep = " << inter.n_kepler
		 << "  n_osc = " << inter.n_osc
		 << ":  " << state_string(inter.descriptor) << endl;
	else {

	    // Gather statistics:

	    real d_ecc = final.ecc - init.ecc;
	    d_ecc_sum += d_ecc;
	    d_ecc_sq_sum += d_ecc*d_ecc;
	    d_ecc_min = min(d_ecc_min, d_ecc);
	    d_ecc_max = max(d_ecc_max, d_ecc);

	    if (debug_flag) {
		cout << "  trial #" << i << ":  d_ecc = " << d_ecc << endl;
		cout << "  minima: " << inter.r_min[0]
		     << "  " << inter.r_min[1]
		     << "  " << inter.r_min[2] << endl;
		cout << "  energy = " << e_init
		     << "  error = " << final.error << endl;
	    }

	    *(de + total_ok) = d_ecc;

	    int ie = (int) (BINS_PER_DECADE * log10(d_ecc)) + N_BINS - 1;
	    if (ie < 0)
		ie = 0;
	    else if (ie >= N_BINS)
		ie = N_BINS - 1;

	    de_histo[ie]++;	    

	    total_ok++;
	    if (out_flag && total_ok % N_OUT == 0) cout << total_ok << endl;
	}
    }

    if (total_ok >= MIN_STATS) {
	real mean = d_ecc_sum / total_ok;
	cout << "  total hits = " << total_ok << "  <d_ecc> = " << mean
	     << "  sigma = " << sqrt(d_ecc_sq_sum / total_ok - mean * mean)
	     << endl;
	cout << "  min d_ecc = " << d_ecc_min
	     << "  max = " << d_ecc_max << endl;

	qsort((void*)de, (size_t)total_ok, sizeof(real), realcomp);

	cout.precision(4);
	real factor = 0.01*total_ok;
	cout << "  10-percentiles:" << endl << " ";
	for (i = 10; i <= 90; i += 10) cout << " " << *(de + (int)(factor*i));
	cout << endl;

	cout << "  histogram (" << BINS_PER_DECADE << " bins per decade):\n";
	int j = 0;
	while (j < N_BINS) {	// cout and printf don't mesh on Solaris 2.3!
	    printf("\t");
	    for (i = 0; i < BINS_PER_DECADE && j < N_BINS; i++, j++)
		printf("  %5d", de_histo[j]);	
	    printf("\n");
	}

	if (dump_flag) {
	    cout << "\n  Raw data (total = " << total_ok << ", sorted):\n";
	    j = 0;
	    cout << "\t";
	    while (j < total_ok) {
		for (i = 0; i < 5 && j < total_ok; i++, j++)
		    printf("  %12.8f", *(de + j));
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
    int  c;

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

    while ((c = pgetopt(argc, argv, "A:dDe:m:M:n:opPr:R:s:t:v:",
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
	    case 'R': init.r_stop = atof(poptarg);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 't': dt_print = atof(poptarg);
		      debug_flag = 1;
		      break;
	    case 'v': init.v_inf = atof(poptarg);
		      break;
            case '?': cerr << "usage: ecc [-A #] [-e #] [-m #] [-M #] [-n #] "
			   << "[-o] [-p] [-P] [-r #] [-R #] [-s #] [-t #] "
			   << "[-v #]"
			   << endl;
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

    cout << "  Start/stop at R > " << init.r_stop
	 << "  accuracy parameter = " << init.eta
	 << "  random seed = " << get_initial_seed()
	 << endl;

    init.r_init_min = init.r_init_max = init.r_stop;

    cpu_init();
    run_trials(n_trials, init, out_flag, planar_flag, 
	       dump_flag, debug_flag, dt_print);
}
