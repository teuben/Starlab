
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  scatter.h: definitions for scattering experiments
 *.............................................................................
 *    version 1:  Dec 1993   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of state structures
 *  2) declaration of an integrator
 *       ....
 *.............................................................................
 */

#ifndef  STARLAB_SCATTER_H
#  define  STARLAB_SCATTER_H

#ifdef USE_MPI
#include  "mpi++.h"
#else
#include  "../src/node/dyn/hdyn/sdyn/evolve/localmpi++.h"
#endif

#include  "../src/node/dyn/hdyn/sdyn/evolve/scatter_hist.h"
#include  "sdyn.h"

#define DATA_COMING_TAG    0
#define DATA_REQUEST_TAG   1
#define WRITE_REQUEST_TAG  2
#define WRITE_READY_TAG    3
#define SEND_HOLD_TAG      4
#define STOP_HOLD_TAG      5

#define N_FINAL 255 // Obscure historical number inherited from scatter3
                    // gives now the number of scenarios, which can be large!
#define N_INTER 255 // see N_FINAL

#define DYNAMICS	1
#define STABILITY	1
#define K_MAX		2

//-----------------------------------------------------------------------------

#if 0
enum intermediate_descriptor {
    non_resonance=0, hierarchical_resonance, democratic_resonance,
    unknown_intermediate, number_of_intermediate_descriptors
};

enum final_descriptor {
  preservation=0, non_preservation, error, stopped, unknown_final,
  number_of_final_descriptors
};
#endif

// Integration parameters:

#define CHECK_INTERVAL  	20.0
#define DEFAULT_ETA		0.05

// Scattering parameters:

#define LARGE_SEPARATION 10.0      // To be determined (e.g. large mass ratios)
#define LARGE_SEPARATION_FACTOR 10.0 // For all mass ratios
#define TIDAL_TOL_FACTOR 1e-6
#define ENERGY_SAFETY_FACTOR 0.01

#define MIN_INITIAL_SEPARATION 	  10.0
#define MAX_INITIAL_SEPARATION 	  VERY_LARGE_NUMBER

#define ENERGY_TOLERANCE 1e-4        // Maximum absolute energy error allowed
#define MERGER_ENERGY_TOLERANCE 1e-3 // Relax tolerance for *relative* error
				     // in case of merger


class scatter_input {
public:
  char init_string[255];

  real  delta_t;
  real  eta;
  real  tidal_tol_factor;
  real  dt_out;
  real  dt_snap;
  real  snap_cube_size;
  real  cpu_time_check;

  int  n_experiments;
  int  seed;
  int  n_rand;
  int  n_rand_inc;
  int  pipe;
  int  debug;
  int  verbose;

  scatter_input() {

    delta_t = dt_out = dt_snap = VERY_LARGE_NUMBER;
    eta = 0.02;    
    tidal_tol_factor = 1.e-6;
    snap_cube_size = 10;
    cpu_time_check = 3600;
    n_experiments = 1;
    seed = n_rand = n_rand_inc = pipe = debug = 0;
    verbose = 0;
  }

  friend ostream& operator<<(ostream& s, scatter_input&);
};

// Structure describing the high-level profile of a scattering experiment:

typedef struct {
    char* init_string;                  // Initialization for mkscat
    int n_string;                       
    real mt;				// mass of target binary
    real mp;				// mass of projectile binary
    real ap;                            // semi major axis of projectile 
    real peri;  			// maximum pericenter distance
    real rho;			        // selected impact parameter
    real rho_sq_min;			// minimum impact parameter squared
    real rho_sq_max;			// maximum impact parameter squared
    real v_inf;                         // velocity at infinity

    real r1;				// radius of primary
    real r2;				// radius of secondary
    real r3;				// radius of third star
    real tidal_tol_factor;		// tidal perturbation at start/stop
    int  rho_flag;			// choice of random or user-specified
	                                // 	impact parameter
    int  ecc_flag;			// choice random or user-specified
					//	eccentricity
    real ecc;				// eccentricity if specified by user
  //phase3_flag phase_flag; 	        // randomization options for angles
  //phase3 phase;			// angles if specified externally
    real eta;				// overall/initial accuracy parameter

  //    intermediate_descriptor3	// specification of the final state
  //    intermediate_target;
  //    final_descriptor3 final_target1;
  //    final_descriptor3 final_target2;

} scatter_profile;

// Simple structure representing the dynamical state of a particle:

typedef struct {
    int index;		// > 0 for a real particle
    real mass;		// >= 0 for a real particle
    vec pos;
    vec vel;
} body;

// Helpers:

void print_bodies(ostream&, body*, int prec = 6);

void initialize_bodies(body *);
void set_kepler_from_sdyn(kepler&, sdyn *, sdyn *);

sdyn* mkscat(int, char**, scatter_profile&);
sdyn* mkscat(char*, scatter_profile&);
sdyn* mkscat(int, char**);
sdyn* mkscat(char*);



// Structures defining the initial, intermediate, and final states of a
// single scattering.  "Initial" structure now modified to allow dual
// use in both scattering and bound configurations.

typedef struct {

    int  np_star;       // number of stars in projectile
    int  nt_star;       // number of stars in target
    real v_inf;		// com projectile velocity at infinity (unit = v_crit)
    real rho;		// com projectile impact parameter

  //    nbodystate * system;
} initial_state;

void scatter(sdyn* b, real eta,
	     real delta_t, real dt_out, real cpu_time_check,
	     real dt_snap, real ttf, real snap_cube_size,
	     int debug,
	     scatter_exp &experiment);

void scatter(sdyn* b, scatter_input input, 
	     scatter_exp &experiment);

int extend_or_end_scatter(sdyn * b, real ttf, bool debug);
void ppn(sdyn* b, ostream & s, int level = 0);
void parse_string(char* s, int& argc, char* argv[]);

#endif


