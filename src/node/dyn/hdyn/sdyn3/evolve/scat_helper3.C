
// scat_helper3.C: Helper functions to aid in setting up and manipulating 
//                 three-body scattering experiments.

// Starlab library function.

#include "scatter3.h"

// Translate enum states into strings:

char * state_string(intermediate_descriptor3 s)
{
    switch(s) {
        case non_resonance:          return "non_resonance";
        case hierarchical_resonance: return "hierarchical_resonance";
        case democratic_resonance:   return "democratic_resonance";
        case unknown_intermediate:   return "unknown_intermediate";
    }
}

char * state_string(final_descriptor3 s)
{
    switch(s) {
        case preservation:     	return "preservation";
        case exchange_1:  	return "exchange_1";
        case exchange_2:  	return "exchange_2";
        case ionization: 	return "ionization";
        case merger_binary_1: 	return "merger_binary_1";
        case merger_binary_2: 	return "merger_binary_2";
        case merger_binary_3: 	return "merger_binary_3";
        case merger_escape_1: 	return "merger_escape_1";
        case merger_escape_2: 	return "merger_escape_2";
        case merger_escape_3: 	return "merger_escape_3";
        case triple_merger: 	return "triple_merger";
        case error: 		return "error";
        case stopped: 		return "stopped";
        case unknown_final:	return "unknown_final";
    }
}

// Pretty-print a body array:

void print_bodies(ostream & s, body * system, int prec)
{
    s.precision(prec);

    s << "  body structure:" << endl;
    for (int k = 0; k < 3; k++) {

	// No need to print lines with zero index or negative mass...

	if (system[k].index > 0 && system[k].mass > 0) {
	    s << "    " << system[k].index << "  "   << system[k].mass;
	    for (int kp = 0; kp < 3; kp++) s << "  "   << system[k].pos[kp];
	    for (int kv = 0; kv < 3; kv++) s << "  "   << system[k].vel[kv];
	    s << endl;
	}
    }
}

// Pretty-print an initial state:

void print_initial(ostream & s, initial_state3 & i,
		   int bod_flag, int prec)
{
    s.precision(prec);

    s << "initial_state:" << endl;
    s << "  m1 = " << 1 - i.m2 << "  m2 = " << i.m2
      << "  m3 = " << i.m3
      << "  r1 = " << i.r1 << "  r2 = " << i.r2 << "  r3 = " << i.r3 << endl;
    s << "  a_binary = " << i.sma << "  e_binary = " << i.ecc << endl;
    s << "  v_inf = " << i.v_inf << "  rho = " << i.rho
      << "  r_init = " << i.r_init;
    if (i.r_stop < VERY_LARGE_NUMBER)
	s << "  r_stop = " << i.r_stop;
    else
	s << "  tidal_tol_factor = " << i.tidal_tol_factor;
    s << endl;
    s << "  phase = " << acos(i.phase.cos_theta) << "  "
      		      << i.phase.phi << "  " << i.phase.psi << "  "
      		      << i.phase.mean_anomaly << endl;

    if (bod_flag) print_bodies(s, i.system, prec);
}

// Pretty-print an intermediate state:

void print_intermediate(ostream & s, intermediate_state3 & i,
			int bod_flag, int prec)
{
    s.precision(prec);

    s << "intermediate_state:" << endl;
    s << "  " << state_string(i.descriptor)
      << "  n_osc = " << i.n_osc
      << "  n_kepler = " << i.n_kepler
      << endl;

    if (i.r_min_min * i.r_min_min < 0.9 * VERY_LARGE_NUMBER)
	{
	s << "  r_min = ";
	for (int k = 0; k < i.n_stars; k++)
	    s << "(" << i.index[k] << ") " << i.r_min[k] << "  ";
	s << "r_min_min = " << i.r_min_min << endl;
	}

    if (bod_flag) print_bodies(s, i.system, prec);
    
}

// Pretty-print a final state:

void print_final(ostream & s, final_state3 & f,
		 int bod_flag, int prec)
{
    s.precision(prec);

    s << "final_state:" << endl;
    s << "  " << state_string(f.descriptor);

    if (f.sma > 0) s << "  a_binary = " << f.sma << "  e_binary = " << f.ecc;
    s << endl;

    s << "  escaper = " << f.escaper ;

    if (f.outer_separation > 0)
      s << "  outer_separation = " << f.outer_separation 
	<< "  outer_virial_ratio = " << f.virial_ratio;
    s << endl;

    s << "  time = " << f.time << "  n_steps = " << f.n_steps
      << "  dE = " << f.error << endl;

    if (bod_flag) print_bodies(s, f.system, prec);

}

// Pretty-print scattering outcomes:

void print_scatter3_outcome(intermediate_state3& inter,
			    final_state3& final,
			    ostream& s)		// default = cerr
{
    s << state_string(inter.descriptor) << " "
      << state_string(final.descriptor) << endl;
}

void print_scatter3_summary(intermediate_state3& inter,
			    final_state3& final,
			    real cpu,
			    ostream& s)		// default = cerr
{
    int p = s.precision(STD_PRECISION);
    s << "  time = " << final.time
      << "  n_steps = " << final.n_steps
      << "  n_osc = " << inter.n_osc
      << "  n_kepler = " << inter.n_kepler
      << endl
      << "  error = " << final.error
      << "  total CPU time = " << cpu << " s"
      << endl;
    s.precision(p);
}

void print_scatter3_report(initial_state3& init,
			   intermediate_state3& inter,
			   final_state3& final,
			   real cpu,		// default = -1
			   int b_flag,		// default = 0
			   ostream& s)		// default = cerr
{
    int flag1 = b_flag, flag2 = b_flag;

    // No body output if b_flag = 0.
    // All body output if b_flag = 1 (cerr) or 3 (file).
    // No init body output if b_flag = 2.
    // No body output if b_flag = 2 and no merger occurred.

    // Note that the intermediate and final body arrays will be the
    // same unless a merger occurred, and neither will be particularly
    // interesting in the non-merger case.

    if (b_flag == 2) {
	flag1 = 0;
	if (final.descriptor < merger_binary_1) flag2 = b_flag = 0;
    }

    print_initial(s, init, flag1);
    print_intermediate(s, inter, flag2);
    print_final(s, final, b_flag);

    if (cpu > 0) s << "  CPU time = " << cpu << endl;
}

// Set up a random phase structure:

void randomize_angles(phase3 &p) // Establish random angles for scattering
{
    p.cos_theta = randinter(-1, 1);
    p.phi = randinter(0, TWO_PI);
    p.psi = randinter(0, TWO_PI);
    p.mean_anomaly = randinter(-PI, PI);
}

// System consists of "non-particles" until properly initialized

void initialize_bodies(body * system)
{
    for (int k = 0; k < 3; k++) {
        system[k].index = 0;
        system[k].mass = -1;
    }
}

// Set up a template initial state:

void make_standard_init(initial_state3 & init)
{
    init.m2 = 0.5;   	        // mass of secondary in target binary
    init.m3 = 0.5;   		// mass of projectile (m1 + m2 = 1)
    init.r1 = 0;   		// radius of star #1
    init.r2 = 0;   		// radius of star #2
    init.r3 = 0;   		// radius of star #3
    init.sma = 1;
    init.ecc = 0;    		// inner binary eccentricity

    init.v_inf = 1;  		// velocity at infinity, in units of v_crit
    init.rho = 0;    		// impact parameter

    // Not used (for bound systems):

    init.a_out = 0;
    init.e_out = 0;

    init.r_init_min = MIN_INITIAL_SEPARATION;
    init.r_init_max = MAX_INITIAL_SEPARATION;
    init.r_stop
          = VERY_LARGE_NUMBER;  // final separation
    init.tidal_tol_factor = DEFAULT_TIDAL_TOL_FACTOR;
    init.eta
          = DEFAULT_ETA;  	// accuracy parameter

    initialize_bodies(init.system);
    init.id = get_initial_seed() + get_n_rand();

    init.cpu_limit = VERY_LARGE_NUMBER;
    init.snap_limit = -1;
}
