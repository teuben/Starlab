#include "dyn.h"

local void hsys_stats(dyn* b,
		      real energy_cutoff = 0.5,
		      bool verbose = true,
		      bool binaries = true,
		      int which_lagr = 1,
		      bool print_time = true,
		      bool compute_energy = false,
		      bool allow_n_sq_ops = false)
{
    sys_stats(b, energy_cutoff, verbose, binaries,
	      which_lagr, print_time, compute_energy, allow_n_sq_ops, false);
//	      print_dstar_params,
//	      print_dstar_stats,
//	      get_energies_with_tidal);	// uses calculate_internal_energies
					// from kira_harp_include.C
}

main(int argc, char **argv)
{
    dyn *b;

    while (b = get_dyn(cin)) {

	put_story(cerr, *b->get_dyn_story());

	check_addstar(b);
//	b->set_use_dstar(true);

	hsys_stats(b);

	rmtree(b);	// causes core dump if use_dstar enabled...
    }
}
