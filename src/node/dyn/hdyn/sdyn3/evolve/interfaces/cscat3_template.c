
/*
 * cscat3_template.c:  C interface to three-body scattering experiments.
 *		       The "c_" routines are C-callable links to Starlab
 *		       C++ functions.
 */

#define   C_ONLY		/* Necessary to avoid C++ parts of header */
#include "scatter3.h"
#include "c_interface.h"

local void print_system(body* system)
{
    int i, k;

    printf("    system:\n");
    for (i = 0; i < 3; i++) {
	printf("        %d (%d):  mass:  %f\n",
	       i, system[i].index, system[i].mass);
	printf("                pos: ");
	for (k = 0; k < 3; k++) printf("  %f", system[i].pos[k]);
	printf("\n");
	printf("                vel: ");
	for (k = 0; k < 3; k++) printf("  %f", system[i].vel[k]);
	printf("\n");
    }
}

void main(int argc, char **argv)
{
    int  seed 	= 0;		/* Seed for random number generator	   */
				/* (0 means system chooses number)	   */
    int  n_rand = 0;		/* Number of times to invoke the generator */
				/* before starting "for real"		   */
    real dt_out
	  = VERY_LARGE_NUMBER;	/* Output time interval			   */
    real dt_snap
	  = VERY_LARGE_NUMBER;	/* Snapshot time interval		   */
    real snap_cube_size = 10;	/* Snapshots only if bodies lie in cube	   */

    real cpu_time_check = 3600;	/* Interval for checking CPU time (sec)	   */

    initial_state3 init;	/* These structures define the initial,	   */
    intermediate_state3 inter;	/* intermediate, and final states of the   */
    final_state3 final;		/* system.  See scatter3.h for details.	   */

    int k;

    /* Initialize the system. */

    c_srandinter(seed, n_rand);
    c_make_standard_init(&init);
    c_initialize_angles(&init, 0, 0, 0.0);

    /* Modify, e.g. with */  init.v_inf = 0.1;
    /*               or  */  init.phase.mean_anomaly = c_randinter(0.0, 2*PI);

    /* Perform the scattering experiment. */

    c_cpu_init();
    c_scatter3(&init, &inter, &final, cpu_time_check,
	       dt_out, dt_snap, snap_cube_size);

    /* Print all available information on the initial state. */

    printf("Initial state:\n");
    printf("    m1         = %f\n", 1-init.m2);	 /* convention: m1+m2 = 1 */
    printf("    m2         = %f\n", init.m2);
    printf("    m3         = %f\n", init.m3);
    printf("    r1         = %f\n", init.r1);
    printf("    r2         = %f\n", init.r2);
    printf("    r3         = %f\n", init.r3);
    printf("    sma        = %f\n", init.sma);
    printf("    ecc        = %f\n", init.ecc);
    printf("    v_inf      = %f\n", init.v_inf);
    printf("    rho        = %f\n", init.rho);
    printf("    r_init_min = %f\n", init.r_init_min);
    printf("    r_init_max = %e\n", init.r_init_max);
    printf("    r_init     = %f\n", init.r_init);
    printf("    r_stop     = %e\n", init.r_stop);
    printf("    tidal_tol  = %e\n", init.tidal_tol_factor);
    printf("    eta        = %f\n", init.eta);
    printf("    phase:\n");
    printf("        cos_theta    = %f\n", init.phase.cos_theta);
    printf("        phi          = %f\n", init.phase.phi);
    printf("        psi          = %f\n", init.phase.psi);
    printf("        mean_anomaly = %f\n", init.phase.mean_anomaly);
    print_system(init.system);

    /* Print all available information on the intermediate state. */

    printf("\nIntermediate state:\n");
    printf("    n_osc      = %d\n", inter.n_osc);
    printf("    n_kepler   = %d\n", inter.n_kepler);
    printf("    n_stars    = %d\n", inter.n_stars);
    printf("    index      =");
    for (k = 0; k < 3; k++ ) printf(" %d", inter.index[k]);
    printf("\n");
    printf("    r_min      =");
    for (k = 0; k < 3; k++ ) printf(" %e", inter.r_min[k]);
    printf("\n");
    printf("    r_min_min  = %e\n", inter.r_min_min);
    printf("    descriptor = "), fflush(stdout);
    c_print_intermediate_descriptor(inter.descriptor);
    print_system(inter.system);

    /* Print all available information on the final state. */

    printf("\nFinal state:\n");
    printf("    sma              = %f\n", final.sma);
    printf("    ecc              = %f\n", final.ecc);
    printf("    outer_separation = %f\n", final.outer_separation);
    printf("    escaper          = %d\n", final.escaper);
    printf("    error            = %e\n", final.error);
    printf("    time             = %f\n", final.time);
    printf("    n_steps          = %d\n", final.n_steps);
    printf("    virial_ratio     = %f\n", final.virial_ratio);
    printf("    descriptor       = "), fflush(stdout);
    c_print_final_descriptor(final.descriptor);
    print_system(final.system);

    printf("\nCPU time = %f\n", c_cpu_time());
}
