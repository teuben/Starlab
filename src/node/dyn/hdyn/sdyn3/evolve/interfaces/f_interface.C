
// Interface to let FORTRAN routines use (C++) Starlab scattering routines.

#include "scatter3.h"
#include "fc_interface.h"

// The FORTRAN common blocks are now precisely C++ structures,
// so no additional conversion should be required, apart from
// Sun's idiotic trailing underscores on global names...

#ifdef FORTRAN_TRAILING_UNDERSCORE
#  define INIT	init_
#  define INTER	inter_
#  define FINAL	final_
#else
#  define INIT	init
#  define INTER	inter
#  define FINAL	final
#endif

extern initial_state3 INIT;
extern intermediate_state3 INTER;
extern final_state3 FINAL;

// WARNING: This is likely to be machine-dependent, and hence unreliable,
// because we cannot be certain how any given compiler will arrange data
// in memory.

// We are ASSUMING here that the FORTRAN common blocks defined in
// f_scatter3.h are aligned by f77 in the same way as the C/C++ init,
// inter, and final structures are aligned by (g)cc.  If problems arise,
// they usually do so in the "system" portions of the structures.

extern "C" {

    /* --------- Initialize, return CPU time: as for C --------- */

    void f_cpu_init()
	{cpu_init();}

    real f_cpu_time()
	{return cpu_time();}


    /* --------- Initialize, return, print random numbers: --------- */

    int f_srandinter(int* seed, int* n_rand)
	{return srandinter(*seed, *n_rand);}

    real f_randinter(real* a, real* b)
	{return randinter(*a, *b);}

    void f_print_initial_random_parameters()
	{c_print_initial_random_parameters();}


    /* --------- Initialize standard initial state for scattering: --------- */

    void f_make_standard_init()
	{make_standard_init(INIT);}

    void f_initialize_angles(int* planar_flag, int* psi_flag, real* psi)
	{c_initialize_angles(&INIT, *planar_flag, (bool) *psi_flag, *psi);}


    /* --------- Print information on scattering: --------- */

    void f_print_scatter3_info(int* Q_flag, int* q_flag, int* b_flag,
			       real* cpu)
	{c_print_scatter3_info(&INIT, &INTER, &FINAL,
			       (bool) *Q_flag, (bool) *q_flag, (bool) *b_flag,
			       *cpu);}

    void f_print_intermediate_descriptor(int* descriptor)
	{cout << "    descriptor =   "
	      << state_string((intermediate_descriptor3) *descriptor) << endl;}

    void f_print_final_descriptor(int* descriptor)
	{cout << "    descriptor       =   "
	      << state_string((final_descriptor3) *descriptor) << endl;}


    /* --------- Perform a scattering experiment: --------- */

    void f_scatter3(real* cpu_time_check,
		    real* dt_out,
		    real* dt_snap,
		    real* snap_cube_size)
	{scatter3(INIT, INTER, FINAL, *cpu_time_check,
		  *dt_out, *dt_snap, *snap_cube_size);}

    // Note: For unknown reasons, the system reports uncleared IEEE
    // errors (Inexact;  Division by Zero; Invalid Operand) when
    // this routine is used on Suns.  The results appear to be fine.

    // This is almost certainly related to data alignment problems
    // in the "system" part of the structures and common blocks,
    // but I'm not sure why that should be the case...
}
