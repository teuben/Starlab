
// Interface to let C routines use (C++) Starlab scattering routines.

#include "scatter3.h"

extern "C" {

    /* --------- Initialize, return CPU time: --------- */

    void c_cpu_init()
	{cpu_init();}

    real c_cpu_time()
	{return cpu_time();}


    /* --------- Initialize, return, print random numbers: --------- */

    int c_srandinter(int seed, int n_rand)
	{return srandinter(seed, n_rand);}

    real c_randinter(real a, real b)
	{return randinter(a, b);}

    void c_print_initial_random_parameters() {
	cerr << "Random seed = " << get_initial_seed()
	     << "  n_rand = " << get_n_rand() << flush;
    }


    /* --------- Simple command-line parsing: --------- */

    int c_pgetopt(int argc, char** argv, char* param_string) /* CL parsing */
	{return pgetopt(argc, argv, param_string);}


    /* --------- Initialize standard initial state for scattering: --------- */

    void c_make_standard_init(initial_state3* init)
	{make_standard_init(*init);}

    void c_initialize_angles(initial_state3* init,
			int planar_flag, bool psi_flag, real psi) {
	randomize_angles(init->phase);
	if (planar_flag == 1) {
	    init->phase.cos_theta = 1;			// Planar prograde
	    if (psi_flag) init->phase.psi = psi * M_PI / 180.0;
	} else if (planar_flag == -1) {
	    init->phase.cos_theta = -1;			// Planar retrograde
	    if (psi_flag) init->phase.psi = psi * M_PI / 180.0;
	}
    }


    /* --------- Print information on scattering: --------- */

    void c_print_scatter3_info(initial_state3* init,
			       intermediate_state3* inter,
			       final_state3* final,
			       bool Q_flag, bool q_flag, int b_flag,
			       real cpu) {
	cerr << ":  ";
	print_scatter3_outcome(*inter, *final, cerr);
	if (Q_flag) print_scatter3_summary(*inter, *final, cpu, cerr);
	if (!q_flag) print_scatter3_report(*init, *inter, *final,
					   cpu, b_flag, cerr);
    }

    void c_print_intermediate_descriptor(intermediate_descriptor3 descriptor)
	{cout << state_string(descriptor) << endl;}

    void c_print_final_descriptor(final_descriptor3 descriptor)
	{cout << state_string(descriptor) << endl;}


    /* --------- Perform a scattering experiment: --------- */

    void c_scatter3(initial_state3* init,
		    intermediate_state3* inter,
		    final_state3* final,
		    real cpu_time_check,
		    real dt_out,
		    real dt_snap,
		    real snap_cube_size) {
	scatter3(*init, *inter, *final, cpu_time_check,
		 dt_out, dt_snap, snap_cube_size);
    }
}

	
