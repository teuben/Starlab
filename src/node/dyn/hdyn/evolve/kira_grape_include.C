
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
//			  but instead included directly into kira.C.
//			  scale_grape.C, and the hdyn version of
//			  sys_stats.C when the executable is built.
//
// Convenient to separate these functions from the rest of kira and
// other code, both for identification and ease of reuse.
//
// USE_GRAPE should only be true during the building of kira_grape4/6,
// scale_grape, and sys_stats_grape, and must be set in the appropriate
// Makefile.
//
// Functions defined:
//
//	bool kira_use_grape
//
//	void kira_calculate_internal_energies
//	void kira_calculate_energies
//	void kira_top_level_energies
//	void kira_calculate_top_level_acc_and_jerk
//	void kira_compute_densities
//	void kira_synchronize_tree
//
// The last six functions in effect act as switches to the appropriate
// GRAPEx routines.  If USE_GRAPE is defined at build time, they cause
// the GRAPE libraries to be loaded.
//
//	grape_calculate_energies
//	grape_calculate_acc_and_jerk
//	grape_calculate_densities

// Note that these functions are declared in hdyn.h, so they are generally
// accessible, but *only* if loaded by an executable including this file...

bool kira_use_grape()
{
#if defined(USE_GRAPE)
    return true;
#else
    return false;
#endif
}

void kira_calculate_internal_energies(hdyn* b,
				      real& epot, real& ekin, real& etot,
				      bool cm,		// default = false
				      bool use_grape)	// default = true
{
    // Compute the total internal energy; also compute the "pot" class
    // datum.

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

    // This function is called by hdyn::merge_nodes() and kira routine
    // calculate_energies_with_external().  Note that use by a member
    // function apparently does *not* mean that all tools using hdyns
    // need to include this file...

#if defined(USE_GRAPE)

    if (use_grape)					// tautology?

	grape_calculate_energies(b, epot, ekin, etot, cm);

    else

	calculate_energies(b, b->get_eps2(), epot, ekin, etot, cm);

#else

    calculate_energies(b, b->get_eps2(), epot, ekin, etot, cm);

#endif

}

void kira_calculate_energies(dyn* b, real eps2, 
			     real &potential, real &kinetic, real &total,
			     bool cm)
{
    // Provide an hdyn function with a calling sequence that can be
    // substituted for the dyn function dyn::calculate_energies when
    // called by sys_stats...  Discard eps2 (--> 0).

    // Called by hdyn/evolve/kira_log.C (kira function print_statistics())
    // and hdyn/util/sys_stats.C (standalone tool), in each case as an
    // argument to the dyn version of sys_stats.

    kira_calculate_internal_energies((hdyn*)b, potential, kinetic, total, cm);
}

void kira_top_level_energies(dyn *b, real eps2,
			     real& potential_energy,
			     real& kinetic_energy)
{
    // Another lookalike, this time to perform the operation of
    // dyn::get_top_level_energies() using the GRAPE if possible.
    // Used only in the standalone tool scale_grape.C as an argument
    // to the dyn function scale.

    real energy;
    kira_calculate_energies(b, eps2,
			    potential_energy, kinetic_energy, energy,
			    true);
}

void kira_calculate_top_level_acc_and_jerk(hdyn ** next_nodes,
					   int n_next,
					   xreal time,
					   bool & restart_grape)
{
    // Switch between GRAPE and non-GRAPE determination of the forces
    // on the particles in the specified list.  Called only by 
    // calculate_acc_and_jerk_for_list() in kira_ev.C.

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

void kira_compute_densities(hdyn* b, vector& cod_pos, vector& cod_vel)
{
    // Density computation (currently limited to GRAPE systems).
    // Called only by log_output() in kira_log.C.

#if defined(USE_GRAPE)

    cerr << "Computing densities using GRAPE..." << endl;
    real cpu0 = cpu_time();

    // The second argument determines the squared radius at which
    // particles are deemed to have zero densities.  This allows
    // discrimination against low-density particles, and also limits
    // costly repeat GRAPE calls.

    grape_calculate_densities(b, 0.1);		// (densities are saved in
						//  particle dyn stories)
    real cpu1 = cpu_time();
    compute_mean_cod(b, cod_pos, cod_vel);
    real cpu2 = cpu_time();

    cerr << "CPU times:  density " << cpu1 - cpu0
	 << "  cod " << cpu2 - cpu1
	 << endl;

#else

    // Skip as too expensive if no GRAPE is available...

    cerr << "Skipping density calculation..." << endl;

#endif

}

void kira_synchronize_tree(hdyn *b)
{
    // GRAPE replacement for synchronize_tree().  Synchronize all
    // top-level nodes.  Called from integrate_list() in kira.C and
    // hdyn::merge_nodes().  Somewhat more elaborate than other
    // functions in this file,/ as the entire algorithm is contained
    // here.

#if defined(USE_GRAPE)

    // Code is similar to that in integrate_list(), but only top-level
    // nodes are considered and we don't check for errors in function
    // correct_and_update.  Possibly should merge this with (part of)
    // integrate_list() and drop synchronize_tree() completely.
    //							 (Steve, 1/02)

    // Make a list of top-level nodes in need of synchronization.

    xreal sys_t = b->get_system_time();

    cerr << endl
	 << "synchronizing tree using GRAPE at time " << sys_t
	 << endl << flush;

    int n_next = 0;
    for_all_daughters(hdyn, b, bi)
	if (bi->get_time() < sys_t) n_next++;

    hdyn **next_nodes = new hdynptr[n_next];
    n_next = 0;
    for_all_daughters(hdyn, b, bi)
	if (bi->get_time() < sys_t) next_nodes[n_next++] = bi;

    // Integrate all particles on the list.  Start by computing forces.
    // (Assume exact = false and ignore_internal = true.)

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	predict_loworder_all(bi, sys_t);
	bi->clear_interaction();
	bi->top_level_node_prologue_for_force_calculation(false);
    }

    bool restart_grape = false;
    grape_calculate_acc_and_jerk(next_nodes, n_next, sys_t, restart_grape);

    int ntop = get_n_top_level();
    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	bi->inc_direct_force(ntop-1);
	bi->top_level_node_epilogue_force_calculation();
	bi->inc_steps();
    }

    // Complete calculation of accs and jerks by correcting for C.M.
    // interactions and applying external fields.

    kira_options *ko = b->get_kira_options();

    for (int i = 0; i < n_next; i++)
	next_nodes[i]->set_on_integration_list();

    if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list) {
	bool reset = false;
	correct_acc_and_jerk(b, reset);			// old version
    } else
	correct_acc_and_jerk(next_nodes, n_next);	// new version

    for (int i = 0; i < n_next; i++)
	next_nodes[i]->clear_on_integration_list();

    if (b->get_external_field() > 0) {

        // Add external forces.

        for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    real pot;
	    vector acc, jerk;
	    get_external_acc(bi, bi->get_pred_pos(), bi->get_pred_vel(),
			     pot, acc, jerk);
	    bi->inc_pot(pot);
	    bi->inc_acc(acc);
	    bi->inc_jerk(jerk);
	}
    }

    // Apply corrector and redetermine timesteps.

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	bi->correct_and_update();
	bi->init_pred();
	bi->store_old_force();
    }

    delete [] next_nodes;

    cerr << endl
	 << "end of synchronization"
	 << endl << flush;

#else

    cerr << endl
	 << "synchronizing tree without GRAPE at time " << b->get_system_time()
	 << endl;

    synchronize_tree(b);

#endif
}
