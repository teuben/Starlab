
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  sigma.h: definitions for determining cross sections
 *.............................................................................
 *    version 1:  March 1994   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of state structures
 *       ....
 *.............................................................................
 */

#ifndef  STARLAB_SIGMAN_H
#  define  STARLAB_SIGMAN_H

#include  "scatter.h"

#define USE_specific_function
#undef USE_specific_function

#define ECC_TOL 2         // Factor by which initial binary periastron must
		          // exceed the sum of the component radii. 

#define RHO_SQ_SCALE 0.5  // Inner disk of radius 1 will be uniformly sampled
#define RHO_SQ_FACTOR 2   // Each new zone doubles the total area

#define N_STEP_BIN 10     // Number of time-step bins
#define N_OSC_BIN 15      // Number of oscillation bins

// Structure describing the result of a series of scattering experiments:

class sigma_input : public scatter_input {
  public:

  real  rho_sq_min;
  real  rho_sq_max;
  real  peri;
  real  pmass;
  real  tmass;
  real  max_trial_density;
  real  v_inf;

  //  int print_counts = FALSE;
  //  int intermediate_sigma = TRUE;

  sigma_input() : scatter_input() {
    v_inf = 1;
    rho_sq_min = rho_sq_max = -1;
    peri = -1;
    pmass = 0;
    tmass = 1;
    max_trial_density = 1;
  }

  friend ostream& operator<<(ostream& s, scatter_input&);
#ifdef USE_MPI
  MPI_Datatype sigma_input::initialize_data_structures_MPI();
#else
  MPI_Datatype sigma_input::initialize_data_structures_MPI() {
    MPI_Datatype dummy;
    return dummy;
}
#endif
};

typedef struct {

    // history
    scatter_hist *hi;

    // Binning:
    real rho_sq_init;           // squared radius of innermost rho zone
    real central_trial_density;	// current density of trials in innermost zone
    int n_zone;                


    // Safety zone quantities:

    int  n_hit_tot;		// total number of hits
    int  n_hit[N_RHO_ZONE_MAX]; // number of hits per rho^2 zone
    int  i_max;                 // i value of the (empty) top zone

    // Diagnostics:

    int  total_trials;
    int  total_steps;
    int  max_steps;

    int  step_counter[N_STEP_BIN];
    int  osc_counter[N_OSC_BIN];

    // "Real" quantities of interest:

    real rho_max;		// maximum impact parameter (from i_max)
    int  trials_per_zone;

    int  n_hits[number_of_scatter_discriptors][N_RHO_ZONE_MAX]; 
    real sigma[number_of_scatter_discriptors];         
    real sigma_err_sq[number_of_scatter_discriptors]; 

    real sigma_total;
    real sigma_total_err_sq;

} sigma_out;

//-----------------------------------------------------------------------------

#if 0

// Standard framework for determining cross-sections:

typedef  int  (*trial_fp)(initial_state3 &,
			  intermediate_state3 &,
			  final_state3 &,
			  real, real, real);
typedef  void (*stat_fp)(scatter_profile &,
			 initial_state3 &,
			 intermediate_state3 &,
			 final_state3 &,
			 int,
			 sigma_out &);
typedef  void (*print_fp)(real, sigma_out &, int);

#endif

void get_sigma(sigma_input &input, MPI_Datatype inputtype,
	       scatter_exp &experiment, MPI_Datatype scatter_exp_type);
#if 0
void get_sigma(real max_dens, scatter_profile &, sigma_out &,
	       real eta, real delta_t, real dt_out, 
	       int debug=0, real cpu_time_check = VERY_LARGE_NUMBER,
	       real dt_snap=VERY_LARGE_NUMBER, real snap_cube_size =0,
	       int scatter_summary_level=0);
#endif
#if 0

void get_sigma3(real, scatter_profile &, sigma_out &,
		int, real, real, real, int,
		trial_fp);
void get_sigma3(real, scatter_profile &, sigma_out &,
		int, real, real, real, int,
		stat_fp);
void get_sigma3(real, scatter_profile &, sigma_out &,
		int, real, real, real, int,
		print_fp);
void get_sigma3(real, scatter_profile &, sigma_out &,
		int, real, real, real, int,
		trial_fp, stat_fp);
void get_sigma3(real, scatter_profile &, sigma_out &,
		int, real, real, real, int,
		trial_fp, stat_fp, print_fp);
#endif

// Helpers:

real zone_area(sigma_out&, int);
real zone_density(sigma_out&, int);
real zone_weight(sigma_out& out, scatter_exp* hi, int i);
                                        // (weight = 1/density, in fact)

void print_profile(ostream&, scatter_profile&, int prec = 6);
void make_standard_profile(scatter_profile &); // Set up a template profile
//void prof_to_init(scatter_profile &, initial_state3 &);

int  n_coll(sigma_out);

void counts_to_sigma(sigma_out &);
void print_all_sigma_counts(sigma_out &, ostream& s = cerr);
void print_sigma(sigma_out &, real);
void print_sigma_counts(sigma_out &);
int specific_counts(scatter_exp *hi, scatter_discriptor discription);


void print_sigma_array(sigma_out, real, int sqrt_flag = 0);
void print_sigma_nonmergers(sigma_out, real);
void print_sigma_mergers(sigma_out out, real);
//void print_sigma_err_array(real[][N_FINAL], real);

// For use in sigma3, rate3, etc:

void summarize_scattering_initial(scatter_profile & prof,
				  int n_rand,
				  real dt_snap,
				  real snap_cube_size);

void summarize_scattering_final(scatter_exp exp,
				int level,
				real cpu);


void single_scatter_init(sdyn *b, scatter_profile & prof,
			 real rho_sq_min, real rho_sq_max,
			 int & n_rand,
			 int scatter_summary,
			 real dt_snap, real snap_cube_size);

int single_scatter(sdyn* b, scatter_input input, 
		   scatter_exp &experiment);

int single_scatter(sdyn *b,
		   scatter_exp &experiment,
		   real eta, real delta_t, real dt_out,
		   real cpu_time_check,
		   real dt_snap,
		   real ttf,
		   real snap_cube_size,
		   int debug);


void single_scatter_stats(scatter_exp* exp,
			  sigma_out & out);

#if 0
void single_scatter_stats(scatter_profile & prof,
			  scatter_exp* exp,
			  int rho_bin_index,
			  sigma_out & out,
			  //			  stat_fp acc_stats,
			  int result);
#endif

int  multiscatter(sigma_out &out, 
		  sigma_input &input, 
		  MPI_Datatype inputtype,
		  scatter_exp &experiment, 
		  MPI_Datatype scatter_exp_type,
		  real &cpu_save, int& scatt_total, real& cpu_total);

int multiscatter(scatter_profile & prof, sigma_out & out,
		 //		 real rho_sq_min, real rho_sq_max, 
		 real eta, real delta_t, real dt_out,
		 real dt_snap, real ttf, real snap_cube_size,
		 real cpu_time_check, real cpu_init, real &cpu_save,
		 int& scatt_total, real& cpu_total, 
		 // stat_fp acc_stats,
		 int debug, int scatter_summary_flag);

sdyn* mkscat(int, char**, sigma_input&);
sdyn* mkscat(char*,  sigma_input &input);

// see in specific_function.C
void print_scatter_specific_information(sdyn *b,
					sigma_input input, 
					scatter_exp experiment);

#endif
