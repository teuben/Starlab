#define USE_GRAPE
=======

       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// sys_stats:  Print out various diagnostic statistics on the input
////             system.  These include:
////
////                 system time, number, mass, mass distribution
////                 relaxation time
////                 system energy
////                 core parameters
////                 lagrangian radii
////                     for quartiles [default]
////                     for ten-percentiles
////                     for "special" choice of Lagrangian masses,
////                         currently 0.005, 0.01, 0.02, 0.05,
////                                   0.1, 0.25, 0.5, 0.75, 0.9
////                 mass distribution by lagrangian zone
////                 anisotropy by lagrangian zone
////                 binary parameters
////
////             In addition, Lagrangian radii are written to the root
////             dyn story.  If N^2_ops	is selected, then core parameters
////             and particle densities are also written to the root
////             and particle dyn stories.
////
////             If N^2_ops is selected, Lagrangian radii are computed
////             relative to the density center, as are binary radial
////             coordinates.  Otherwise, the modified center of mass
////             (center of mass with outliers excluded)  is used.
////                 
//// Options:    -b    specify level of binary statistics [2]
////                       0:		  none
////                       1 (or no arg): short binary output
////                       2:		  full binary output
////             -B    include binary evolution [get from snapshot]
////             -e    recalculate the total energy (requires -n) [yes]
////             -l    specify percentile choice [2]
////                       0:             quartiles (1-3)
////                       1:             10-percentiles (10-90)
////                       2 (or no arg): nonlinear Lagrangian masses
////                                      (0.5, 1, 2, 5, 10, 25, 50, 75, 90%)
////             -n    perform/don't perform actions requiring O(N^2)
////                   operations (e.g. computation of energy, core radius,
////                   density, bound pairs) [yes]
////             -o    pipe system to cout [no]
////
//// Notes:  The hdyn sys_stats *does* take tidal fields into account.
////
////         If GRAPE is available, we will compute the energies even if
////         the "-n" flag is set false.

#include "hdyn.h"

#ifdef TOOLBOX

#include "../evolve/kira_grape_include.C"	// GRAPE-specific functions:
						// Invoke by setting USE_GRAPE
						// at compilation time.

local void check_set_tidal(hdyn *b)
{
    b->set_tidal_field(0);

    bool verbose = false;			// turn on debugging here
    if (verbose) cerr << endl;

    real initial_mass = get_initial_mass(b, verbose);
    real initial_r_virial = get_initial_virial_radius(b, verbose);
    real initial_r_jacobi = get_initial_jacobi_radius(b, initial_r_virial,
						      verbose);
    if (initial_r_jacobi > 0) {
	int tidal_field_type = 0;
	set_tidal_params(b, verbose,
			 initial_r_jacobi,
			 initial_mass,
			 tidal_field_type);
    }

    if (verbose) cerr << endl;
}

local void sys_stats(hdyn* b,
		     real energy_cutoff,
		     bool verbose,
		     bool binaries,
		     bool long_binary_output,
		     int  which_lagr,
		     bool print_time,
		     bool compute_energy,
		     bool allow_n_sq_ops)
{
    // Make densities using GRAPE (if available).

    vector cod_pos, cod_vel;
    // Reinstated compute_densities by SPZ on April 2001
    compute_densities(b, cod_pos, cod_vel);

    // Suppress densities for now, and write density time to prevent
    // sys_stats from doing the calculation.

    putrq(b->get_dyn_story(), "density_time", b->get_real_system_time());
    for_all_nodes(hdyn, b, bb)
	putrq(bb->get_dyn_story(), "density_time", b->get_real_system_time());

    // Invoke the dyn version with appropriate hdyn extensions, as in kira.

    sys_stats(b, energy_cutoff,
	      verbose,
	      binaries,
	      long_binary_output,
	      which_lagr,
	      print_time,
	      compute_energy,
	      allow_n_sq_ops,
	      get_energies_with_external, // uses calculate_internal_energies
					  // from kira_grape_include.C
	      print_dstar_params,
	      print_dstar_stats);
}

// Main code follows dyn/util/sys_stats.C version.

main(int argc, char **argv)
{
    check_help();

    bool binaries, long_binary_output, B_flag, verbose, out, n_sq, calc_e;
    int which_lagr;

    // (Defaults are set in dyn/util/sys_stats.C)

    if (!parse_sys_stats_main(argc, argv,
			      which_lagr,
			      binaries, long_binary_output, B_flag,
			      calc_e, n_sq, out, verbose)) {
	get_help();
	exit(1);
    }

    if (!n_sq) {

#ifndef USE_GRAPE
	calc_e = false;
#endif

    }

    // Loop over input until no more data remain.

    hdyn *b;
    int i = 0;

    while (b = get_hdyn(cin)) {

	// Set up tidal and stellar structures.

	check_set_tidal(b);

	check_addstar(b);
	if (B_flag || check_kira_flag(b, "kira_evolve_binaries"))
	    b->set_use_dstar(true);

	if (i++ > 0) cerr << endl;

	sys_stats(b, 0.5, verbose, binaries, long_binary_output,
		  which_lagr, true, calc_e, n_sq);

	if (out) put_node(cout, *b);

	rmtree(b);	// causes core dump if B_flag is enabled...
    }
}

#endif
