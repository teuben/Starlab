
// sig_helper.C: Engine for determining cross-sections in
//               n-body scattering.

// Differences from the previous version:
//
//	more streamlined scan over impact parameter
//	decoupling of initialization, integration, and analysis functions
//	parallel operation, eventually

// Starlab library function.

#include "stdio.h"
#include "sigma.h"
#include "sigma_MPI.h"

#define DEFAULT_TIDAL_TOL_FACTOR 1e-6

// Pretty-print a scatter profile:
void print_profile(ostream & s, scatter_profile & p, int prec)
{
    s.precision(prec);

    s << "scatter_profile: " << endl
      << p.init_string << endl;
    s << "   Mt = " << p.mt << "  mp = " << p.mp << "  ap = " << p.ap 
      << "  peri = " << p.peri << "  rho = " << p.rho << "  rho^2 = " 
      << p.rho_sq_min << "  " << p.rho_sq_max << endl;
      
    if (p.v_inf >= 0) s << "  v_inf = " << p.v_inf;

    s << endl;

    if (p.rho_flag)
	s << "   rho2_min = " << p.rho_sq_min
	  << " rho2_max = " << p.rho_sq_max << endl;

    if (p.ecc_flag)
	s << "   binary_ecc = " << p.ecc << endl;

#if 0
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
#endif	
}

// Set up a template profile:

void make_standard_profile(scatter_profile & prof)
{
  prof.init_string = "-M 1 -v 2 -t -p";    // head on collision
  prof.mt = 1;
  prof.mp = 1;				// mass of projectile binary
  prof.ap = 1;                            // semi major axis of projectile 
  prof.peri = -1;  			// maximum pericenter distance
  prof.rho = -1;			        // selected impact parameter
  prof.rho_sq_min = -1;			// minimum impact parameter squared
  prof.rho_sq_max = -1;			// maximum impact parameter squared
  prof.v_inf = 1;                         // velocity at infinity
    
  prof.r1 = 0;				// radius of primary
  prof.r2 = 0;				// radius of secondary
  prof.r3 = 0;				// radius of third star
  prof.tidal_tol_factor = DEFAULT_TIDAL_TOL_FACTOR;
  prof.rho_flag = 0;			// choice of random or user-specified
  // 	impact parameter
  prof.ecc_flag = 0;			// choice random or user-specified
  //	eccentricity
  prof.ecc = 0;				// eccentricity if specified by user
  prof.eta = DEFAULT_ETA;			// overall/initial accuracy parameter
  
}

// Turn a scatter profile into a default initial state:

void sdyn_to_prof(sdyn *b, scatter_profile & prof)
{
  
    if (prof.mt > 1) {
	cerr << "prof_to_init: prof.mt = " << prof.mt << " > 1" << endl;
	exit(1);
    }
}


int n_coll(sigma_out out) {

    int total = 0;	

    // summ all encounters with a collision
    for_all_scatter_hist(scatter_hist, out.hi->get_first(), hi) {
      if(hi->get_scatter_discriptor() == two_body_collision ||
	 hi->get_scatter_discriptor() == multiple_body_collision)
	cerr << " COLLISION" << endl;
	total += hi->get_n_found();
    }

    return total;
}

#define output_header_1a \
  "   pres   exch ex_ion  s_ion  n_ion  t_ion coll_b coll_i  tripl n_tupl\n"
//       0      1      2      3      4      5      6      7      8      9

#define output_header_1b \
  "   2_coll   n_coll   unknown     error  stopped not_spec\n"
//         0        1         2         3        4        5

#define output_header_2 \
  "   pres   exch   ion  merg_b merg_e merg_3  error   stop\n"

#define output_header_3 \
  "   pres   exch ex_ion  s_ion  n_ion  t_ion  tripl n_tupl\n"
//       0      1      2      3      4      5      6      7


#define output_header_4 \
  "   coll_bin   coll_ion     s_coll     n_coll\n"
//           0          1          2          3


static char *f_lab[] = {"non_res", "res", "stop", "unknown"};
// (static here to keep Sun CC happy...)

// NOTE: All "print_sigma" functions start on the current line (i.e. they
//       do not skip a line at the start) and skip a line at the end.

char dummy_string[20];

void print_sigma_counts(sigma_out &out) {

    int n_scenarios = out.hi->get_last()->get_id_scenario();
    int n_total[number_of_scatter_discriptors];
    int i;
    for(i=0; i<number_of_scatter_discriptors; i++)
      n_total[i] = 0;

	// Combine some categories:

    for_all_scatter_hist(scatter_hist, out.hi->get_first(), hi) {
      for(i=preservation; i<number_of_scatter_discriptors; i++) 
	n_total[i] += hi->get_specific_counts((scatter_discriptor)i);
    }

    n_total[two_body_collision] = 0;
    n_total[multiple_body_collision] = 0;
    for_all_scatter_hist(scatter_hist, out.hi->get_first(), ha) {
      if(ha->get_n_coll() == 1)
	n_total[two_body_collision] += ha->get_n_found();
      else if(ha->get_n_coll() > 1)
	n_total[multiple_body_collision] += ha->get_n_found();
    }

    cerr << output_header_1a;
    int i_fin;
    for (i_fin = preservation; i_fin < two_body_collision; i_fin++) {

	    cerr << "  ";

	    sprintf(dummy_string, "%d", n_total[i_fin]);

	    // Use exactly 5 chars (+ 1 space) unless integer is very large.

	    for (i = 0; i < Starlab::max(0, 5 - (int)strlen(dummy_string)); i++)
		cerr << " ";
	    
	    cerr << dummy_string;
	}
    cerr << endl << endl;

    cerr << output_header_1b;
    for (i_fin = two_body_collision;
	 i_fin < number_of_scatter_discriptors; i_fin++) {

	    cerr << "  ";

	    sprintf(dummy_string, "%d", n_total[i_fin]);

	    // Use exactly 7 chars (+ 1 space) unless integer is very large.

	    for (int i = 0; i < Starlab::max(0, 7 - (int)strlen(dummy_string)); i++)
		cerr << " ";
	    
	    cerr << dummy_string;
	}
    cerr << endl << endl;
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


void print_sigma_array(sigma_out out, real scale_factor, int sqrt_flag) {

    cerr << output_header_1a;

    real s_total[number_of_scatter_discriptors];

    for(int r=0; r<2; r++) {

      for(int i=0; i<number_of_scatter_discriptors; i++)
	s_total[i] = 0;

      // Combine some categories:

      for_all_scatter_hist(scatter_hist, out.hi->get_first(), hi) {
	for(int i=preservation; i<number_of_scatter_discriptors; i++) 
	  if(r==0 && !hi->get_resonance() && !hi->get_stop())
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	  else if(r==1 && hi->get_resonance() && !hi->get_stop())
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	  else if(r==3 && hi->get_stop())
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
      }

      for (int i_fin = 0; i_fin < two_body_collision; i_fin++) {
	if (r==0 && i_fin == 0) 
	  fprintf(stderr, "       ");
	else {
	  real temp = (sqrt_flag == 0 ? s_total[i_fin]
		       : sqrt(s_total[i_fin]));
	  print_formatted(scale_factor * temp);
	}
      }
      cerr << "   " << f_lab[r] << endl;
    }
    cerr << endl;
}

void print_sigma_nonmergers(sigma_out out, real v2) {

    cerr << output_header_1a << endl; 

    real s_total[number_of_scatter_discriptors];
    real s_err_sq[number_of_scatter_discriptors];

    for (int i_int = 0; i_int < 2; i_int++) {

      for(int i=0; i<number_of_scatter_discriptors; i++) {
	s_total[i] = 0;
	s_err_sq[i] = 0;
      }

	// Cross-sections:

      for_all_scatter_hist(scatter_hist, out.hi->get_first(), hi) {

	for(int i=preservation; i<number_of_scatter_discriptors; i++) 
	  if(i_int==0 && !hi->get_resonance() && !hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
	  else if(i_int==1 && hi->get_resonance() && !hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
	  else if(i_int==2 && hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
      }

      int i_fin;
      for (i_fin = 0; i_fin <  two_body_collision; i_fin++) 
	if (i_int == 0 && i_fin == 0) 
	  fprintf(stderr, "       ");
	else 
	  print_formatted(v2 * s_total[i_fin]);
      cerr << "                 " << f_lab[i_int] << endl;

      // Errors:

      for (i_fin = 0; i_fin <  two_body_collision; i_fin++) 
	//	  if (i_fin!= collision_binary && i_fin!=collision_ionization) {
	if (i_int == 0 && i_fin == 0) 
	  fprintf(stderr, "       ");
      //	    else if(!(i_int==0 && i_fin==stable_higher_order))  
	else 
	  print_error(v2 * sqrt(s_err_sq[i_fin]));
      //	}
      cerr << "\n\n";
    }
}

void print_sigma_mergers(sigma_out out, real v2) {

    cerr << output_header_4 << endl; 

    real s_total[number_of_scatter_discriptors];
    real s_err_sq[number_of_scatter_discriptors];

    int i_fin, i_int, i;
    for (i_int = 0; i_int < 3; i_int++) {

      for(i=0; i<number_of_scatter_discriptors; i++) {
	s_total[i] = 0;
	s_err_sq[i] = 0;
      }


	// Cross-sections:

      for_all_scatter_hist(scatter_hist, out.hi->get_first(), hi) {
	for(i=preservation; i<number_of_scatter_discriptors; i++) 
	  if(i_int==0 && !hi->get_resonance() && !hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
	  else if(i_int==1 && hi->get_resonance() && !hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
	  else if(i_int==2 && hi->get_stop()) {
	    s_total[i] += hi->get_specific_sigma((scatter_discriptor)i);
	    s_err_sq[i] += hi->get_specific_sigma_err_sq((scatter_discriptor)i);
	  }
      }

	for (i_fin = 0; i_fin <  unknown; i_fin++) 
	  if (i_fin== collision_binary    || 
	      i_fin==collision_ionization ||
	      i_fin==two_body_collision   || 
	      i_fin==multiple_body_collision) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "               ");
	    else 
	      print_formatted(v2 * s_total[i_fin]);
	}
	cerr << "                 " << f_lab[i_int] << endl;

	// Errors:

	for (i_fin = 0; i_fin <  unknown; i_fin++) 
	  if (i_fin== collision_binary    || 
	      i_fin==collision_ionization ||
	      i_fin==two_body_collision   || 
	      i_fin==multiple_body_collision) {
	    if (i_int == 0 && i_fin == 0) 
		fprintf(stderr, "              ");
	    else 
	      print_error(v2 * sqrt(s_err_sq[i_fin]));
	}
	cerr << "\n\n";
    }
}

local void create_temp(char* temp, int length, char* string)
{
  //    for (int k = 0; k < length; k++) temp[k] = ' ';
  //    strcpy(temp, string);
  //    for (int k = 0; k < length; k++) if (temp[k] == '\0') temp[k] = ' ';
  //    temp[length-1] = '\0';
}

void print_sigma(sigma_out & out, real v2) {


  cerr << endl;
  int i;
  for(i=0; i<72; i++) cerr << "-";
  cerr << endl;
  for(i=0; i<72; i++) cerr << "-";
  cerr << endl;
  
    cerr << "Scenario count: " << endl;
    out.hi->put_scatter_hist(cerr);
    
    cerr << endl;
    for (i = 0; i < 72; i++) cerr << "-";
    cerr << endl;
    
    cerr.precision(6);
    cerr << "\n   rho_max = " << out.rho_max
	 << "  i_max = " << out.i_max << endl;
    cerr << "   central_trial_density = " << out.central_trial_density
	 << "  total_trials = " << out.total_trials << endl;

    cerr << "   average steps = " << out.total_steps / (real) out.total_trials
	 << "  maximum = " << out.max_steps << endl;
    cerr << "   timesteps    (bin factor = 10): ";
    for (i = 0; i < N_STEP_BIN; i++) cerr << " " << out.step_counter[i];
    cerr << endl;

    int coll_total = n_coll(out);

    cerr << "   oscillations (bin factor =  2): " ;
    //    for (i = 0; i < N_OSC_BIN; i++) cerr << " " << out.osc_counter[i];
    //    cerr << endl;

    //    int jmax = (int)triple_merger;
    //    if (coll_total == 0) jmax = (int)ionization;

    for (i = 0; i < N_OSC_BIN; i++) cerr << " " << out.osc_counter[i];
    cerr << endl;
    

    // Raw counts:
    
    cerr << endl << "raw counts:\n\n";
    print_sigma_counts(out);

    // Cross-sections:

    cerr << "v^2 sigma / (pi a^2):\n\n";
    print_sigma_array(out, v2, 0);

    // More detail and errors:

    // Non-mergers first (0-3, N_FINAL-2):

    print_sigma_nonmergers(out, v2);

    // Then mergers (4-10):

    //    if (coll_total > 0) 
    //      print_sigma_mergers(out, v2);

    cerr << endl;
    for (i = 0; i < 72; i++) cerr << "-";
    cerr << endl;

    cerr << flush;
}

// print_counts: print out enough information to determine cross-sections
//		 and errors, and to combine them with other runs.
//------------------------------------------------------------------------

// Encapsulate the following functions -- efficiency is not an issue here.

real zone_area(sigma_out& out, int i)	{        // (without the PI)

    real rho_sq_min = 0, rho_sq_max;

    rho_sq_max = out.rho_sq_init * pow((real)RHO_SQ_FACTOR, i);
    if (i > 0) rho_sq_min = rho_sq_max / RHO_SQ_FACTOR;

    return rho_sq_max - rho_sq_min;
}

real zone_density(sigma_out& out, int i)
{
  return out.trials_per_zone / zone_area(out, i);
}

real zone_weight(sigma_out& out, scatter_exp* hi, int i)
{
  return zone_area(out, i) / hi->get_nhits(i);
}

//------------------------------------------------------------------------

// counts_to_sigma: Convert raw counts into cross-sections and errors.
//		    This will update sigma from the current counts,
//		    making NO other changes in the sigma_out structure.

void counts_to_sigma(sigma_out & out)
{
    out.sigma_total = 0;
    out.sigma_total_err_sq = 0;

    int i_zone;
    out.sigma_total = 0;
    out.sigma_total_err_sq = 0;

    for_all_scatter_hist(scatter_hist, out.hi->get_first(), ha) {
      i_zone = 0;
      do {

	if(ha->get_nhits(i_zone)) {
	  real weight = zone_weight(out, ha, i_zone);
	  real delta_sigma =  weight * ha->get_nhits(i_zone);
	  real delta_err_sq = weight * weight
	    * ha->get_nhits(i_zone);

	  ha->inc_sigma(delta_sigma);
	  ha->inc_sigma_err_sq(delta_err_sq);

	  out.sigma_total += delta_sigma;
	  out.sigma_total_err_sq += delta_err_sq;
	  
	}
      }
      while(i_zone++<out.i_max); 
    }

    int i;
    for_all_scatter_hist(scatter_hist, out.hi->get_first(), ho) {
      for (i = 0; i < N_STEP_BIN; i++) 
	out.step_counter[i] += ho->get_step_counter(i);
      for (i = 0; i < N_OSC_BIN; i++) 
	out.osc_counter[i] += ho->get_osc_counter(i);
    }
} 

// Elementary statistics:

real max_t = 0;					       // global variables!!
real max_err = 0;

local void print_status(real v_inf, sigma_out & out, int debug) {

    cerr.precision(6);
    cerr.setf(ios::fixed, ios::floatfield);

    cerr << "counts_to_sigma(out);"<<endl;
    counts_to_sigma(out);
 
    cerr << "sigma[" << v_inf << "] = " << out.sigma_total
 	 << " +/- " << sqrt(out.sigma_total_err_sq)
	 << "  v^2 sigma = " << out.sigma_total*v_inf*v_inf
 	 << " +/- " << sqrt(out.sigma_total_err_sq)*v_inf*v_inf
	 << endl;
 
//    if (debug < 0)
	cerr << "  n_total = "	<< out.total_trials
	     << "  max_t = "	<< max_t
	     << "  max_err = "	<< max_err
	     << "  n_hit_tot = "<< out.n_hit_tot
	     << endl;
    }

// summarize_scattering_initial: Print a line specifying a
//				 scattering experiment.

void summarize_scattering_initial(scatter_profile & prof,
				  int n_rand,
				  real dt_snap,
				  real snap_cube_size)
{
    cerr.precision(6);
    cerr << "\nscatter"
	 << " -s " << get_initial_seed()
	 << " -N " << n_rand;

    PRL(prof.rho);
    cerr << prof.init_string << endl;

#if 0
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
#endif

    cerr.precision(6);
    if (dt_snap < VERY_LARGE_NUMBER) {
	cerr << " -D " << dt_snap
	     << " -C " << snap_cube_size;
    }

    cerr << endl;

//    print_initial(cerr, init);
}

// summarize_scattering_final: Print a line specifying a scattering result.

void summarize_scattering_final(scatter_exp exp,
				int level,
				real cpu)
{
  cerr << "outcome:  " << exp << &endl;
  //    print_scatter_outcome(inter, final, cerr);
  //    if (level > 1) print_scatter3_summary(inter, final, cpu, cerr);
  //    print_final(cerr, final);
}

//=========================================================================
//
// Functions to initialize, perform, and analyze a scattering experiment.
//
//=========================================================================

// single_scatter_init: Initialize a scattering experiment in the specified
//                      rho range with the given v_inf.

void single_scatter_init(sdyn* b, scatter_profile & prof,
			 real rho_sq_min, real rho_sq_max,
			 int & n_rand,
			 int scatter_summary, 
			 real dt_snap, real snap_cube_size)
{

    n_rand = get_n_rand();

    // Initialize the scattering experiment:

    sdyn_to_prof(b, prof);

    // Choose rho randomly from the allowed range:

    if (prof.rho_flag == 0)
	prof.rho = sqrt(randinter(rho_sq_min, rho_sq_max));
    else
	prof.rho = rho_sq_min;

    if (scatter_summary > 0)
	summarize_scattering_initial(prof, n_rand, dt_snap, snap_cube_size);
}


// single_scatter: Perform a scattering experiment with the specified init.


local void pp(sdyn* b, ostream & s, int level = 0) {

    s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    if(b != b->get_root()) {
      b->pretty_print_node(s);
      s << " \t"<< b->get_mass() << " \t"
		<< b->get_pos() << " (r= " << abs(b->get_pos()) << ")   " 
		<< b->get_vel() << " (v= " << abs(b->get_vel()) << ")" << endl;
      //	<< "r= " << abs(b->get_pos()) << "    " 
      //	<< "v= " << abs(b->get_vel()) << endl;
    }

    for (sdyn * daughter = b->get_oldest_daughter();
	 daughter != NULL;
	 daughter = daughter->get_younger_sister())
	pp(daughter, s, level + 1);	
}

int single_scatter(sdyn* b, scatter_input input, 
		   scatter_exp &experiment) {

  /*
    real eta = input.eta;
    real delta_t = input.delta_t;
    real dt_out = input.dt_out;
    real cpu_time_check = input.cpu_time_check;
    real dt_snap = input.dt_snap;
    real ttf = input.tidal_tol_factor;
    real snap_cube_size = input.snap_cube_size;
    int debug = input.debug;
    
  return single_scatter(b, experiment,
  eta, delta_t, dt_out,
			cpu_time_check,
			dt_snap,
			ttf, 
			snap_cube_size,
			debug);
			}

			int single_scatter(sdyn *b,
		   scatter_exp &experiment,
		   real eta, real delta_t, real dt_out,
		   real cpu_time_check,
		   real dt_snap,
		   real ttf, 
		   real snap_cube_size,
		   int debug) {
  */
  // Perform a scattering experiment with initialized nbody b
    // result = false if the outcome is a "non-resonant preservation".
    // (terminology from scatter3
  
  scatter(b, input, experiment);
  //    scatter(b, eta, delta_t, dt_out, cpu_time_check, dt_snap, 
  //	    ttf, snap_cube_size, debug, experiment);

    // Return a "hit" if the result was
    //
    //		(1) not a flyby (preservation non_resonance)
    //		(2) not an error, stopped, or unknown_final.
    //
    // Modify this return statement to change the definition of a "hit".
    // Probably should take into account the size of the perturbation?

    // return true if preservation and non_resonant
    if((experiment.get_scatter_discriptor()==unknown  ||
	experiment.get_scatter_discriptor()==error    ||
	experiment.get_scatter_discriptor()==stopped) ||
       (experiment.get_scatter_discriptor()==preservation &&
	!experiment.get_resonance()))
      return 0;
    else {
      experiment.inc_n_hits(experiment.get_nzone());
    }

    experiment.inc_n_hit(experiment.get_nzone());
    return 1;
}

// single_scatter_stats: Accumulate information on scattering outcomes.

void single_scatter_stats(scatter_exp* exp,
			  sigma_out & out)
			  //		 stat_fp acc_stats,
{
    // Accumulate statistics:

    max_t = Starlab::max(max_t, exp->get_time());
    max_err = Starlab::max(max_err, abs(exp->get_energy_error()));

    // Accumulate results on data not available to higher levels:

    //    exp->set_nzone(out.n_zone);
    //    exp->inc_n_hits(out.n_zone);

    out.total_steps += exp->get_n_steps();
    out.max_steps = Starlab::max(out.max_steps, exp->get_n_steps());
    real logn = log10((real)Starlab::max(1, exp->get_n_steps()));
    int index = Starlab::min(N_STEP_BIN-1, (int)logn);
    exp->inc_step_counter(index);

    logn = log10((real)Starlab::max(1, exp->get_form_changes())) / log10(2.0);
    index = Starlab::min(N_OSC_BIN-1, (int)logn);
    exp->inc_osc_counter(index);

    // Convenient to update these global counters here also:

    out.total_trials++;
    //    exp->inc_n_hit(out.n_zone);
    int result = 1;
    out.n_hit_tot += result;

    //------------------------------------------------------------------
    // Optional call to user-supplied function for more detailed
    // statistics gathering. The function must maintain its own internal
    // data in a way that can be extracted later.
    //------------------------------------------------------------------

    //    if (acc_stats)
      //	(*acc_stats)(prof, init, inter, final, rho_bin_index, out);

    out.hi->add_scatter_hist(*exp, out.n_zone);
    //    out.hi->put_scatter_hist(cerr);

}

//=========================================================================

// get_sigma:  Run n-body scatterings for a given velocity at infinity
//             and return cross-sections for all possible outcomes.
//             Stop when the specified trial density is exceeded.
//
//	       The trials are performed by function "trial" and debugging
//             information is printed out by the function "print".
//
// In addition to determining cross-sections, get_sigma may also be used
// as a convenient means of driving a more specialized function acc_stats,
// which can gather additional statistics.  Get_sigma guarantees that the
// data sent to acc_stats are properly sampled over impact parameter space.

void get_sigma(sigma_input &init, MPI_Datatype inputtype,
	       scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  int i; 
  int k = 0;
  real init_dens = init.max_trial_density;

  int intermediate_sigma = 1;
  if (intermediate_sigma==1) {
    int min_dens = 4;
    while (init_dens > min_dens) {
      k++;
      init_dens /= 4;
    }
  }

  sigma_out out;    

  int random_seed = srandinter(init.seed, init.n_rand);

  // initialize scatter_hist
  sdyn *b = mkscat(init.init_string, init);
  make_tree(b, !DYNAMICS, STABILITY, K_MAX, init.debug);
  cerr << "Initial configuration: " << endl;
  b->set_name("root");
  vec center = b->get_pos();
  print_structure_recursive(b, 0., center, true, true, 4);
  out.hi = initialize_scatter_hist(b);
  delete b;
  b = NULL;   // safety

  init.n_rand_inc = get_n_rand() - init.n_rand; 
  cerr << "Random seed: input seed " << init.seed
       << "  actual seed = " << random_seed << endl
       << "             n_rand = " << init.n_rand
       << "  n_rand_inc = " << init.n_rand_inc << endl;

  int scatt_total = 0;	// Maintain these for use in parallel case
  real cpu_total = 0;

  for (real max_dens = init_dens; k >= 0; k--, max_dens *= 4) {

    out.total_trials = 0;
    out.n_hit_tot = 0;
    out.max_steps = 0;
    out.total_steps = 0;
    for (i = 0; i < N_STEP_BIN; i++) out.step_counter[i] = 0;
    for (i = 0; i < N_OSC_BIN;  i++) out.osc_counter[i] = 0;

    real cpu_init = cpu_time(); // The CPU time when we entered get_sigma3
    real cpu_save = cpu_init;

    out.i_max = 1;  // i value of safety bin, to be sampled, but supposed
                    // to contain no events, i.e. rho bins run from 0 to 
                    // i_max, with i_max empty

    // Take gravitational focusing into account in the allowed range
    // of rho squared.

    real v_sq = init.v_inf * init.v_inf;

    real m_total = init.pmass + init.tmass;
    if (v_sq > 0)
      out.rho_sq_init = RHO_SQ_SCALE * (1 + 2*m_total/v_sq);
    else
      out.rho_sq_init = RHO_SQ_SCALE;

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

    real tmp = 0.5 + init.max_trial_density * RHO_SQ_SCALE 
             / pericenter_target_factor;
    out.trials_per_zone = (int)tmp;
    out.central_trial_density = (out.trials_per_zone / RHO_SQ_SCALE)
					* pericenter_target_factor;

    // Some useful diagnostic feedback (before the fact):

    if (init.debug != 0) cerr << "==========\n";
    else cerr << endl;

    cerr << "get_sigma:  v_inf = " << init.v_inf
	 << ",  trials_per_zone = " << out.trials_per_zone << endl;

    real rho_sq_min, rho_sq_max;


    for (out.n_zone = 0, rho_sq_min = 0, rho_sq_max = out.rho_sq_init;
	 out.n_zone <= out.i_max;
	 out.n_zone++, rho_sq_min = rho_sq_max, rho_sq_max *= RHO_SQ_FACTOR) {

        out.rho_max = sqrt(rho_sq_max);
	// added (line below) for MPI operation.
	experiment.set_nzone(out.n_zone);
        init.n_experiments = out.trials_per_zone;
        init.rho_sq_min = rho_sq_min;
        init.rho_sq_max = rho_sq_max;
	init.peri = -1;
        int n_hits = multiscatter(out, init, inputtype,
				  experiment, scatter_exp_type,
				  cpu_save, scatt_total, cpu_total);

	if (init.debug)
	  cerr << "rho_zone = " << out.rho_max
	       << ",  result = " << n_hits << endl;

	if (out.n_zone == out.i_max && n_hits > 0) out.i_max++;

	if (init.debug)
	    cerr << ",  out.i_max = " << out.i_max << endl;

	//print_sigma(out, init.v_inf*init.v_inf);
    }

    cerr << "\ntotal scatterings = " << scatt_total
	 << ",  CPU time = " << cpu_total << endl;

    //    terminate_processors(debug);	// PVM or dummy call

    //    if (abs(debug) == 1) (*print)(prof.v_inf, out, debug);

    out.rho_max = sqrt(rho_sq_max);

    //out.hi->put_scatter_hist(cerr);
    counts_to_sigma(out);
    print_sigma(out, init.v_inf * init.v_inf);

  }


  terminate_all_processors();

}


//=========================================================================
