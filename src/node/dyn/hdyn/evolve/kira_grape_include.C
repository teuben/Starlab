
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// kira_grape_include.C:  All conditional code that depends on GRAPE
//			  availability.  Not compiled into any library,
//			  but instead included directly into kira and
//			  the hdyn version of sys_stats.
//
// USE_GRAPE is true only during the building of kira_grape4/6
// and sys_stats_grape, and must be set in the appropriate Makefile.
//
// Functions defined:
//
//	bool kira_use_grape
//
//	void calculate_internal_energies
//	void kira_calculate_top_level_acc_and_jerk
//	void compute_densities
//
// The last three functions in effect act as switches to the appropriate
// GRAPEx routines:
//
//	grape_calculate_energies
//	grape_calculate_acc_and_jerk
//	grape_calculate_densities

bool kira_use_grape()
{
#if defined(USE_GRAPE)
    return true;
#else
    return false;
#endif
}

void kira_calculate_energies(dyn* b, real eps2, 
			     real &potential, real &kinetic, real &total,
			     bool cm)
{
    // Coerce hdyn::calculate_internal_energies into the same calling
    // sequence as dyn::calculate_energies, for use by sys_stats...
    // Discard eps2 (--> 0).

    // Function calculate_internal_energies() is not a library function.
    // It is compiled directly into kira via kira_grape_include.h,
    // allowing selection of the GRAPE/non-GRAPE option at compile time.

    calculate_internal_energies((hdyn*)b, potential, kinetic, total, cm);
}

void kira_top_level_energies(dyn *b, real eps2,
			     real& potential_energy,
			     real& kinetic_energy)
{
    // Another lookalike, this time to perform the operation of
    // dyn::get_top_level_energies() using the GRAPE if possible.

    real energy;
    kira_calculate_energies(b, eps2,
			    potential_energy, kinetic_energy, energy,
			    true);
}

void calculate_internal_energies(hdyn* b,
				 real& epot, real& ekin, real& etot,
				 bool cm,		// default = false
				 bool use_grape)	// default = true
{
    // Compute the total internal energy; also compute the "pot" class data.

    // Notes from Steve (8/99):
    //
    //	- grape_calculate_energies recomputes hdyn::pot, but does
    //	  *not* include the tidal terms.
    //
    //  - new code uses the hdyn version of calculate_energies() (see
    //    ../util/hdyn_tt.C), which sets the hdyn::pot member data,
    //    and also omits the tidal terms.
    //
    // The cm flag specifies that we should use the center-of-mass
    // approximation, i.e. compute the top-level energies only.
    // This is what we want for scale.  Implemented for GRAPE by
    // Steve, 7/01.

#if defined(USE_GRAPE)

    if (use_grape)					// tautology?

	grape_calculate_energies(b, epot, ekin, etot, cm);

    else

	calculate_energies(b, b->get_eps2(), epot, ekin, etot, cm);

#else

    calculate_energies(b, b->get_eps2(), epot, ekin, etot, cm);

#endif

}

void kira_calculate_top_level_acc_and_jerk(hdyn ** next_nodes,
					   int n_next,
					   xreal time,
					   bool & restart_grape)
{
#if defined(USE_GRAPE)

    grape_calculate_acc_and_jerk(next_nodes, n_next,
				 time, restart_grape);
    restart_grape = false;

#else

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	if (bi->is_top_level_node())
	    bi->top_level_node_real_force_calculation();
    }

#endif
}

void compute_densities(hdyn* b, vector& cod_pos, vector& cod_vel)
{
#if defined(USE_GRAPE)

    cerr << "Computing densities using GRAPE..." << endl;
    real cpu0 = cpu_time();

    // The second argument determines the squared radius at which
    // particles are deemed to have zero densities.  This allows
    // discrimination against low-density particles, and also limits
    // costly repeat GRAPE calls.

    grape_calculate_densities(b, 0.1);		// (densities are saved in
						// particle dyn stories)

    real cpu1 = cpu_time();

    compute_mean_cod(b, cod_pos, cod_vel);
    real cpu2 = cpu_time();

    cerr << "CPU times:  density " << cpu1 - cpu0
	 << "  cod " << cpu2 - cpu1
	 << endl;

#endif
}
