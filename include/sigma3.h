       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  sigma3.h: definitions for determining cross sections
 *.............................................................................
 *    version 1:  March 1994   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of state structures
 *       ....
 *.............................................................................
 */

#ifndef  STARLAB_SIGMA3_H
#  define  STARLAB_SIGMA3_H

#include  "scatter3.h"

#define ECC_TOL 2         // Factor by which initial binary periastron must
		          // exceed the sum of the component radii. 

#define RHO_SQ_SCALE 0.5  // Inner disk of radius 1 will be uniformly sampled
#define RHO_SQ_FACTOR 2   // Each new zone doubles the total area
#define N_RHO_ZONE_MAX 20 // Maximum number of rho zones

#define N_STEP_BIN 10     // Number of time-step bins
#define N_OSC_BIN 15      // Number of oscillation bins

// Structure describing the high-level profile of a scattering experiment:

typedef struct {
    real m2;				// mass of binary secondary
    real m3;				// mass of third star
    real r1;				// radius of primary
    real r2;				// radius of secondary
    real r3;				// radius of third star
    real v_inf;				// velocity at infinity
    real tidal_tol_factor;		// tidal perturbation at start/stop
    int  rho_flag;			// choice of random or user-specified
	                                // 	impact parameter
    real rho_min;			// minimum impact parameter
    real rho_max;			// maximum impact parameter
    int  ecc_flag;			// choice random or user-specified
					//	eccentricity
    real ecc;				// eccentricity if specified by user
    phase3_flag phase_flag; 	        // randomization options for angles
    phase3 phase;			// angles if specified externally
    real eta;				// overall/initial accuracy parameter

    intermediate_descriptor3		// specification of the final state
    intermediate_target;
    final_descriptor3 final_target1;
    final_descriptor3 final_target2;

} scatter_profile;

// Structure describing the result of a series of scattering experiments:

typedef struct {

    // Binning:

    real rho_sq_init;           // squared radius of innermost rho zone
    real central_trial_density;	// current density of trials in innermost zone

    // Safety zone quantities:

    int  n_hit_tot;		// total number of hits
    int  n_hit[N_RHO_ZONE_MAX]; // number of hits per rho^2 zone
    int  i_max;                 // i value of the (empty) top zone

    // Diagnostics:

    int  total_trials;
    int  total_steps;
    int  max_steps;
    int  step_counter[N_STEP_BIN];

    int  osc_counter[N_OSC_BIN][N_FINAL];

    // "Real" quantities of interest:

    real rho_max;		// maximum impact parameter (from i_max)
    int  trials_per_zone;

    int  n_hits[N_INTER][N_FINAL][N_RHO_ZONE_MAX]; // breakdown of hits by type
    real sigma[N_INTER][N_FINAL];         // breakdown of cross-section by type
    real sigma_err_sq[N_INTER][N_FINAL]; // breakdown of squared errors by type

    real sigma_total;
    real sigma_total_err_sq;

} sigma_out;

//-----------------------------------------------------------------------------

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

void get_sigma3(real, scatter_profile &, sigma_out &,
		int debug = 0,
		real cpu_time_check = VERY_LARGE_NUMBER,
		real dt_snap =  VERY_LARGE_NUMBER,
		real snap_cube_size = 0,
		int scatter_summary_flag = 0);

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

// Helpers:

real zone_area(sigma_out&, int);
real zone_density(sigma_out&, int);
real zone_weight(sigma_out&, int);	// (weight = 1/density, in fact)

void print_profile(ostream&, scatter_profile&, int prec = 6);
void make_standard_profile(scatter_profile &); // Set up a template profile
void prof_to_init(scatter_profile &, initial_state3 &);

int  n_coll(int[][N_FINAL][N_RHO_ZONE_MAX], int);

void counts_to_sigma(sigma_out &);
void print_all_sigma3_counts(sigma_out &, ostream& s = cerr);
void print_sigma3(sigma_out &, real);
void print_sigma3_counts(int[][N_FINAL][N_RHO_ZONE_MAX], int);
void print_sigma3_array(real[][N_FINAL], real, int sqrt_flag = 0);
void print_sigma3_err_array(real[][N_FINAL], real);
void print_sigma3_nonmergers(real[][N_FINAL], real[][N_FINAL], real);
void print_sigma3_mergers(real[][N_FINAL], real[][N_FINAL], real);

// For use in sigma3, rate3, etc:

void summarize_scattering_initial(initial_state3 & init,
				  int n_rand,
				  real dt_snap,
				  real snap_cube_size);

void summarize_scattering_final(intermediate_state3 & inter,
				final_state3 & final,
				int level,
				real cpu);

void single_scatter_init(scatter_profile & prof,
			 real rho_sq_min, real rho_sq_max,
			 initial_state3 & init, int & n_rand,
			 int scatter_summary,
			 real dt_snap, real snap_cube_size);

void set_decc_crit(real decc);

int single_scatter(initial_state3 & init,
		   intermediate_state3 & inter,
		   final_state3 & final,
		   real cpu_time_check,
		   real dt_snap,
		   real snap_cube_size);

void single_scatter_stats(scatter_profile & prof,
			  initial_state3 & init,
			  intermediate_state3 & inter,
			  final_state3 & final,
			  int rho_bin_index,
			  sigma_out & out,
			  stat_fp acc_stats,
			  int result);

int multiscatter3(scatter_profile & prof, sigma_out & out,
		  real rho_sq_min, real rho_sq_max, int rho_zone,
		  real dt_snap, real snap_cube_size,
		  real cpu_time_check, real cpu_init, real &cpu_save,
		  int& scatt_total, real &cpu_total, stat_fp acc_stats,
		  int debug, int scatter_summary);

#endif
