
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// The help information below is copied directly from the dyn version.
// The only difference between this program and the dyn version is that
// possible GRAPE and dstar extensions are used here.

//// sys_stats:  This is the "hdyn" (kira output) version of sys_stats.
////             It contains all of the "dyn" sys_stats functionality (and
////             will work with standard dyn data files), but includes
////             additional functions specific to the hdyn class.
////
////             Print out various diagnostic statistics on the input system.
////             These include:
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
////                   density, bound pairs) [no]
////             -o    pipe system to cout [no]

#include "hdyn.h"

#ifdef TOOLBOX

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

    bool conf = false;
    // Loop over input until no more data remain.

    hdyn *b;
    int i = 0;

    set_hdyn_check_timestep(false);
    while (b = get_hdyn()) {

	if (!conf) {

	    unsigned int config = kira_config(b);	// default settings
	    kira_print_config(config);

	    if (!n_sq) {

	        if (!b->has_grape()) {
		    // calc_e = false;
		}
	    }

	    conf = true;
	}

	check_addstar(b);
	check_set_external(b, true);	// true ==> verbose output
	cerr << endl;

	if (B_flag || check_kira_flag(b, "kira_evolve_binaries"))
	    b->set_use_dstar(true);

	if (i++ > 0) cerr << endl;

	// Simply invoke the dyn function with appropriate hdyn
	// extensions, as in kira_log.C.

	sys_stats(b,
		  0.5,			// energy cutoff
		  verbose,
		  binaries,
		  long_binary_output,
		  which_lagr,
		  true,			// print_time
		  calc_e,
		  n_sq,
		  kira_calculate_energies,
		  print_dstar_params,
		  print_dstar_stats);

	if (out) {
	    b->log_history(argc, argv);
	    put_node(b);
	}

	rmtree(b);	// causes core dump if B_flag is enabled...
    }
}

#endif
