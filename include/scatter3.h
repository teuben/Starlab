/*
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~
*/

/*
 *  scatter3.h: definitions for scattering experiments
 *.............................................................................
 *    version 1:  Dec 1993   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of state structures
 *  2) declaration of an integrator
 *       ....
 *.............................................................................
 *
 * Note use of C-style comments, as this may conceivably be included in
 * C or C++ programs.
 */

#ifndef  STARLAB_SCATTER3_H
#  define  STARLAB_SCATTER3_H

#ifndef C_ONLY
#  include  "sdyn3.h"
#else
#  include  "c_stdinc.h"
#endif

/*---------------------------------------------------------------------------*/

/* Possible states of a three-body run: */

#define N_INTER 4

enum intermediate_descriptor3 {
    non_resonance, hierarchical_resonance, democratic_resonance,
    unknown_intermediate
};

#define N_FINAL 14

enum final_descriptor3 {
    preservation, exchange_1, exchange_2, ionization,
    merger_binary_1, merger_binary_2, merger_binary_3,
    merger_escape_1, merger_escape_2, merger_escape_3,
    triple_merger, error, stopped, unknown_final
};

/*
// The inner binary always lies in the (x-y) plane, with its long
// axis along the x-axis. The normal to the plane of the outer orbit
// makes an angle theta with respect to the z axis; its projection
// onto the (x-y) plane makes an angle phi with the x-axis--that is,
// theta and phi are the usual spherical polar coordinates describing
// the normal vector.

// All orbital motions are, by definition, counterclockwise around the
// normal vector. Thus, planar prograde orbits have theta = 0, planar
// retrograde orbits have theta = pi.

// The angle psi measures rotation in the plane of the outer orbit.
// If the outer orbit lies in the (x-y) plane, psi measures the
// the orientation of the periastron, measured counterclockwise from
// to the x-axis. For an outer orbit in the (x-z) plane, psi measures
// angle from x towards -z; for the (y-z) plane, psi is from y towards z.
*/

/* Phase angles defining the three-body scattering configuration: */

typedef struct {
    real cos_theta;
    real phi;
    real psi;
    real mean_anomaly;
} phase3;

/* Flags to distinguish between random and specific angles: */

typedef struct {
    int  cos_theta;
    int  phi;
    int  psi;
    int  mean_anomaly;
} phase3_flag;

/*---------------------------------------------------------------------------*/

/* Integration parameters: */

#define CHECK_INTERVAL  20
#define DEFAULT_ETA	0.05

/* Scattering parameters: */

#define LARGE_SEPARATION 10        /* To be refined (e.g. large mass ratios) */
#define LARGE_SEPARATION_FACTOR 10 /* For all mass ratios */
#define ENERGY_SAFETY_FACTOR 0.01

#define MIN_INITIAL_SEPARATION 	  10
#define MAX_INITIAL_SEPARATION 	  VERY_LARGE_NUMBER

#define ENERGY_TOLERANCE 1e-4        /* Maximum absolute energy error allowed */
#define MERGER_ENERGY_TOLERANCE 1e-3 /* Relax tolerance for *relative* error */
				     /* in case of merger */

/* Simple structure representing the dynamical state of a particle.
 * Place the index at the end for alignment reasons explained in
 * .../scatter3/interfaces/f_interface.h.
 */

typedef struct {
    real mass;		/* >= 0 for a real particle */
    real pos[3];	/* Can't use vectors here because g++ 2.5.8 */
    real vel[3];	/* complains that they are improperly initialized... */
    int index;		/* > 0 for a real particle */
} body;

/*
// Structures defining the initial, intermediate, and final states of a
// single scattering.  "Initial" structure now modified to allow dual
// use in both scattering and bound configurations.
*/

typedef struct {
    real m2;		/* mass of secondary in target binary (m1 + m2 = 1) */
    real m3;		/* mass of incoming projectile			    */
    real r1;		/* radius of primary				    */
    real r2;		/* radius of secondary				    */
    real r3;		/* radius of third star				    */
    real sma;		/* inner binary semi-major axis			    */
    real ecc;		/* inner binary eccentricity			    */

    /* Outer orbit (unbound): */

    real v_inf;		/* projectile velocity at infinity (unit = v_crit)  */
    real rho;		/* projectile impact parameter			    */

    /* Outer orbit (bound): */

    real a_out;		/* outer binary semi-major axis			    */
    real e_out;		/* inner binary eccentricity			    */

    real r_init_min;	/* min initial distance from projectile to target   */
    real r_init_max;	/* max initial distance from projectile to target   */
    real r_init;	/* actual initial separation			    */
    real r_stop;	/* max permitted distance from escaper to binary    */
    real tidal_tol_factor; /* tidal perturbation at start/stop		    */
    phase3 phase;	/* phase angles specifying the initial state	    */
    real eta;		/* accuracy parameter				    */
    body system[3];	/* array giving details of initial state	    */
    int id;		/* identifier					    */

    real cpu_limit;	/* handy to keep some limits here...		    */
    int snap_limit;

} initial_state3;

typedef struct {
    int  n_osc;			/* number of "oscillations" in min(r)	    */
    int  n_kepler;		/* number of analytic continuations	    */
    int  n_stars;		/* final number of stars		    */
    int  index[3];		/* final labels				    */
    real r_min[3];		/* minimum interparticle separations	    */
    real r_min_min;		/* absolute minimum separation		    */
    enum intermediate_descriptor3 descriptor;
    body system[3];		/* details of intermediate state	    */
    int id;			/* identifier (should agree with init.id)   */
    int n_snap;			/* snapshot counter			    */
} intermediate_state3;

typedef struct {
    real sma;			/* final binary semi-major axis		    */
    real ecc;			/* final binary eccentricity		    */
    real outer_separation;      /* final separation between third star	    */
                                /* and binary (if any)			    */
    enum final_descriptor3 descriptor;    
    int  escaper;		/* ID of escaping star (0 if none exists)   */
    real error;			/* rel. energy error (unit=binary energy)   */
    xreal time;			/* termination time			    */
    int  n_steps;		/* number of integration steps		    */
    real virial_ratio;		/* final ratio K.E./|P.E.| of outer orbit   */
    body system[3];		/* array giving details of final state	    */
				/* (non-particle <==> index = 0, mass < 0)  */
    int id;			/* identifier (should agree with init.id)   */
} final_state3;

/*---------------------------------------------------------------------------*/

#ifndef C_ONLY

/* Function to perform a single scattering experiment: */

void  scatter3(initial_state3&, intermediate_state3&, final_state3&,
	       real cpu_time_check = 3600,
	       real dt_out = VERY_LARGE_NUMBER,
	       real dt_snap = VERY_LARGE_NUMBER,
	       real snap_cube_size = 10,
	       real dt_print = VERY_LARGE_NUMBER,
	       sdyn3_print_fp p = NULL);

/* Helpers: */

char * state_string(intermediate_descriptor3); /* If state is X, print "X" */
char * state_string(final_descriptor3);

void print_bodies(ostream&, body*, int prec = 6);
void print_initial(ostream&, initial_state3&,
		   int bod_flag = 0, int prec = 6);
void print_intermediate(ostream&, intermediate_state3&,
			int bod_flag = 0, int prec = 6);
void print_final(ostream&, final_state3&,
		 int bod_flag = 0, int prec = 6);

void print_scatter3_outcome(intermediate_state3&,
			    final_state3&,
			    ostream& s = cerr);
void print_scatter3_summary(intermediate_state3&,
			    final_state3&,
			    real,
			    ostream& s = cerr);
void print_scatter3_report(initial_state3&,
			   intermediate_state3&,
			   final_state3&,
			   real cpu = -1,
			   int bod = 0,
			   ostream& s = cerr);

void randomize_angles(phase3 &); /* Establish random angles for scattering */

void initialize_bodies(body *);
void make_standard_init(initial_state3 &);   /* Set up template initial state */

// In scat_init.C:

void set_orientation(kepler &k, phase3 &p);
sdyn3 * set_up_dynamics(real m2, real m3, kepler & k1, kepler & k3);
void sdyn3_to_system(sdyn3 * root, body * system);

// In extend_or_end.C:

void merge_collisions(sdyn3 * b);

int extend_or_end_scatter(sdyn3 * b,
			  initial_state3& init,
			  intermediate_state3 * inter = NULL,
			  final_state3 * final = NULL);

#endif
#endif
