
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions associated with stellar evolution.
//
// Externally visible function:
//
//	bool evolve_stars

#include "hdyn.h"
#include "star/dstar_to_kira.h"

#define MINIMUM_MASS_LOSS 1.e-16  // Do not allow for mass loss below
                                  // the numerical precision

#define MINIMUM_REPORT_MASS_LOSS  1.e-6

local void dissociate_binary(hdyn* bi)
{
    if (bi->get_kepler() == NULL) {
	cerr << "dissociate_binary: no kepler!\n";
	return;
    }

    // Discard kepler if sudden mass loss took place.

    // The orbital phase of the binary is effectively randomized
    // before the sudden mass loss is applied, since the binary
    // dynamics is evolved forward to the current system time.

    if (bi->get_kira_diag()->report_binary_mass_loss) {
	int p = cerr.precision(INT_PRECISION);
	cerr << "\nsudden mass loss from binary, ";
	cerr << bi->format_label() << endl;

	PRC(bi->get_time()); PRL(bi->get_system_time());
	PRL(bi->get_kepler()->get_time());

	// pp3(bi->get_top_level_node(), cerr);

	if (has_dstar(bi))
	    ((double_star*)(bi->get_parent()
			      ->get_starbase()))
			      ->dump(cerr);
	// bi->get_kepler()->print_all(cerr);
	cerr.precision(p);
    }

    // The kepler structure has not yet been modified due
    // to rapid binary evolution (e.g. supernova), so the
    // component positions and velocities retain their values
    // prior to the start of the evolution step in which the
    // supernova occurred.  Any mass transfer occurring during
    // the step, but prior to the supernova, has already been
    // accounted for.  Component masses will be adjusted
    // individually later (as two single stars).

    // A combination of fast and several slow mass-loss episodes
    // is currently handled as a single slow episode followed by
    // a fast episode.

    bi->update_dyn_from_kepler();	// updates kepler to current time

    // Note that this update will in general produce a tidal error since it
    // changes the phase of an unperturbed binary.  However, this error will
    // be absorbed in "kira_counters->de_total" after reinitialization.

    // bi->set_first_timestep();
    // Not necessary because entire system is scheduled
    // to be reinitialized.

    delete bi->get_kepler();

    bi->set_kepler(NULL);
    bi->get_binary_sister()->set_kepler(NULL);

    bi->set_unperturbed_timestep(-VERY_LARGE_NUMBER);
    bi->get_binary_sister()->set_unperturbed_timestep(-VERY_LARGE_NUMBER);

//    if (bi->get_kira_diag()->report_binary_mass_loss) {
        int p = cerr.precision(HIGH_PRECISION);
	cerr << endl << "dissociate_binary: deleted kepler for "
	     << bi->format_label() << ":" << endl;
	PRI(4); PRC(bi->get_time()); PRL(bi->get_system_time());
        cerr.precision(p);
//    }

    if (has_dstar(bi)) {
	bool update_dynamics[2] = {false, false};
	create_or_delete_binary(bi->get_parent(),
				update_dynamics);	// Defaults probably
							// safe here because
							// merger won't occur.
	
    }
		
    // NOTE: In case of sudden mass loss, probably have
    // to recompute accs and jerks, at least of neighbors.
}



real cpu;

local void print_start_evolution(char* s, real t)
{
    cerr << "\n===============\n"
	 << s << " evolution start at time " << t
	 << endl;
    cpu = cpu_time();
}

local void print_end_evolution(char* s, bool correct_dynamics)
{
    real delta_cpu = cpu_time() - cpu;
    PRC(delta_cpu);
    PRL(correct_dynamics);
    cerr << s << " evolution end\n===============\n";
}

bool evolve_stars(hdyn* b,
		  int full_dump)	// default = 0
{
    bool correct_dynamics = false;

    // The following if statement was moved to evolve_system() 4/99.
    // All stars (sstar and dstar) are now updated once we enter
    // this function.

    // if (fmod(b->get_system_time(), dt_sstar) == 0.0 &&

    if (b->get_use_sstar()) {
	
	if (b->get_kira_diag()->report_stellar_evolution)
	    print_start_evolution("Stellar", b->get_system_time());

	correct_dynamics |= stellar_evolution(b);

	if (b->get_kira_diag()->report_stellar_evolution)
	    print_end_evolution("Stellar", correct_dynamics);

	b->get_kira_counters()->step_stellar++;
    }

    if (b->get_use_dstar()) {

	if (b->get_kira_diag()->report_stellar_evolution)
	    print_start_evolution("Binary", b->get_system_time());

	correct_dynamics |= binary_evolution(b, full_dump);

	if (b->get_kira_diag()->report_stellar_evolution)
	    print_end_evolution("Binary", correct_dynamics);
    }

    // Evolution is over.  See if we need to correct the system for
    // mass loss, mergers, etc.  By construction, dynamics is synchronized
    // (apart from unperturbed binaries) before entering this function.
    // Binary evolution has just updated unperturbed binaries to their
    // current times ( < system_time in most cases).

    if (correct_dynamics) {

	// Correct_dynamics is true only if a significant change in
	// the mass or semi-major axis of some star or binary has
	// occurred.  However, once triggered, *all* pending binary
	// mass loss is corrected for.

	b->get_kira_counters()->step_correct++;

	real time = b->get_system_time();

	if (b->get_kira_diag()->report_stellar_evolution) {
	    cerr << "Correcting dynamics at system time "
		 << time << endl;
	}

	// Shouldn't be necessary to synchronize, since all nodes
	// but unperturbed binaries will just have been advanced
	// (see evolve_system in kira.C).
	//
	// Note that (kira_)synchronize_tree() does NOT touch unperturbed
	// binaries (handles only top-level nodes.)

	predict_loworder_all(b, b->get_system_time());	    // unnecessary??
	synchronize_tree(b);

	if (b->get_kira_diag()->report_stellar_evolution)
	    cerr << "After synchronize_tree" << endl;

	if (b->get_kira_diag()->report_stellar_evolution)
	    cerr << "After predict_loworder_all" << endl;

	real mass0 = total_mass(b);

	real epot0, ekin0, etot0;
	calculate_energies_with_external(b, epot0, ekin0, etot0);

	if (b->get_kira_diag()->report_stellar_evolution) {
	    cerr << "After calculate_energies_with_external"<<endl;
	    PRC(epot0); PRC(ekin0); PRL(etot0);
	}

	// Change stellar masses after evolution.

	if (b->get_kira_diag()->report_stellar_evolution)
	    cerr << "\n----------\nCorrecting dynamical masses..." << endl;

	real de_kick = 0;

	for_all_leaves(hdyn, b, bi) {
	    real dm = 0;

	    // Order of business:	slow mass loss
	    //				fast mass loss and dissociation
	    //				single stars

	    if (bi->get_kepler()) {
		if (has_dstar(bi)) {	// elder binary component only...

		    hdyn * parent = bi->get_parent();

		    // Update all binaries for slow mass loss.

		    dm = get_total_mass(parent) - parent->get_mass();

		    // Note: dm is the CHANGE in mass of the binary system.
		    // It may be positive or negative.

		    // Sudden_mass_loss is the mass LOST by the system
		    // in an impulsive fashion.  It is always positive,
		    // by convention.

		    // Determine mass changes due to fast and slow processes
		    // (positive numbers mean mass increase!).

		    real dm1_fast = -sudden_mass_loss(bi);
		    real dm2_fast = -sudden_mass_loss(bi->get_binary_sister());
		    real dm_fast  = dm1_fast + dm2_fast;
		    real dm_slow  = dm - dm_fast;

		    if (dm != 0 || dm_slow != 0 || dm_fast != 0) {

			if (b->get_kira_diag()->report_stellar_evolution &&
			    (abs(dm)      >= MINIMUM_REPORT_MASS_LOSS ||
			     abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS ||
			     abs(dm_fast) >= MINIMUM_REPORT_MASS_LOSS)) {

			    cerr << "Binary evolution mass loss from "
				 << bi->format_label();
			    cerr << ", parent = "
				 << bi->get_parent()->format_label()
				 << ", at time " << bi->get_time() << endl;
			    PRC(dm); PRC(dm_slow); PRC(dm_fast);
			    cerr << " (" << dm1_fast <<" "<< dm2_fast<<")"
				 <<endl;
			}
		    }
			
		    // Note: sudden_mass_loss() is reset to 0 after use.

		    if (b->get_kira_diag()->report_stellar_evolution &&
			abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS &&
			b->get_kira_diag()->report_binary_mass_loss) {

			b->get_kira_counters()->step_dmslow++;

			int p = cerr.precision(INT_PRECISION);
			cerr << "slow mass loss from binary "
			     << bi->format_label() << "  "; PRL(dm_slow);

			PRC(bi->get_time()); PRL(bi->get_system_time());
			PRL(bi->get_kepler()->get_time());

			// pp3(bi->get_top_level_node(), cerr);

			((double_star*)(parent->get_starbase()))->dump(cerr);

			cerr.precision(p);
		    }

		    // The kepler and dyn structures currently reflect the
		    // binary configuration at the start of the current
		    // evolution step.  Update them now to reflect changes
		    // due to binary evolution.

		    // Only correct for slow mass loss here.  Fast mass loss
		    // is accounted for after the binary is dissociated.

		    update_kepler_from_binary_evolution(parent, dm_fast);

		    if (b->get_kira_diag()->report_stellar_evolution &&
			abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS &&
			b->get_kira_diag()->report_binary_mass_loss)
			cerr << "After kepler update..." << flush;

		    update_dyn_from_binary_evolution(parent, bi, dm1_fast,
						                 dm2_fast);

		    if (b->get_kira_diag()->report_stellar_evolution &&
			abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS &&
			b->get_kira_diag()->report_binary_mass_loss)
			cerr << "dyn updated..." << flush;

		    // Neither update_kepler_from_binary_evolution nor
		    // update_dyn_from_binary_evolution modify the parent
		    // mass.  Correct it here, changing component offsets
		    // accordingly.

		    correct_leaf_for_change_of_mass(parent, dm_slow);

		    // Correcting the parent but not the components in this
		    // way effectively isotropizes the mass loss.

		    if (b->get_kira_diag()->report_stellar_evolution &&
			abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS &&
			b->get_kira_diag()->report_binary_mass_loss) {
			cerr << "leaf corrected" << endl;
			// pp3(bi->get_top_level_node(), cerr);
		    }

		    if (bi->get_kepler()) {

			// Kepler structure still exists.  Set the new
			// unperturbed timestep (must extend at least as
			// far as system_time and preferably past apocenter).

			real usteps = bi->get_unperturbed_steps();

#if 0
			cerr << "in evolve_stars for "
			     << bi->format_label() << "  ";
			PRL(usteps);
#endif

			if (usteps > 0) {
			    usteps *= bi->get_timestep();
			    bi->set_unperturbed_timestep(usteps);
			    bi->get_binary_sister()
			      ->set_unperturbed_timestep(usteps);
			}

		    } else {

#if 0
			cerr << "evolve_stars: no kepler for "
			     << bi->format_label() << endl;
#endif
		    }
			    
		    if (b->get_kira_diag()->report_stellar_evolution &&
			abs(dm_slow) >= MINIMUM_REPORT_MASS_LOSS &&
			b->get_kira_diag()->report_binary_mass_loss) {

			cerr << "After unperturbed timestep: "<<endl;

			PRC(bi->get_time());
			PRL(bi->get_unperturbed_timestep());
			PRL(bi->get_time() + bi->get_unperturbed_timestep());
			PRL(bi->get_system_time());
		    }

		    // Start to handle any fast evolution by dissociating
		    // the binary.  (The rest is done by treating the
		    // components as single stars.)

		    if (dm_fast != 0) {
			if (b->get_kira_diag()->report_stellar_evolution) {
			    cerr << "evolve_stars: SN in binary: ";
			    PRL(dm_fast);
			}

			// Note that we always dissociate a binary following
			// fast mass loss (unperturbed motion may subsequently
			// restart).  Phase is effectively randomized, so
			// we shouldn't trust acc, jerk, or dt...  Other
			// than merger, this is the only way a kepler can
			// be deleted in this function.

			b->get_kira_counters()->step_dmfast++;
			dissociate_binary(bi);

			// Do NOT use the time step associated with bi, as
			// it refers to the phase at which the unperturbed
			// motion began and may be inappropriate to the
			// phase now.  Recompute using kepler_step().

			real kep_dt = kepler_step(bi);
			real dt = bi->get_timestep();

			while (dt > kep_dt) dt /= 2;

			bi->set_timestep(dt/2);		// (2 = safety factor)
		    }

		} else if (bi->is_leaf() && !bi->get_use_dstar()) {

		    // Pretty suspect to have stellar evolution with
		    // no binary evolution, although it is possible
		    // that the dstar has already been deleted...

		    if (b->get_kira_diag()->report_stellar_evolution)
			cerr << "evolve_stars: SN in non-dstar binary "
			     << bi->format_label() << endl;

		    dissociate_binary(bi);

		    // Recompute time step (see notes above).

		    real kep_dt = kepler_step(bi);
		    real dt = bi->get_timestep();

		    while (dt > kep_dt) dt /= 2;

		    bi->set_timestep(dt/2);
		}

	    }	// end if (bb->get_kepler()) { ...


	    // Deal with mass loss (single stars or binary CM).

	    // If a merger occurred, everything except a possible change
	    // in mass of the parent node should have been taken care of
	    // in binary_evolution.  The old center of mass is now a leaf;
	    // its stellar mass should reflect the effects of mass loss,
	    // just as with other single stars.

	    if (bi->is_valid()) {		     // (just in case...)

		vector dv = anomalous_velocity(bi);  // |dv| > 0 <==> supernova
		real dv_sq = square(dv);

		dm = get_total_mass(bi) - bi->get_mass();

		hdyn *bt = NULL;
		if (full_dump && dv_sq > 0) {

		    // Print out the state of bi before the supernova.

		    bt = bi->get_top_level_node();
		    put_node(cout, *bt, false, 1);
		}

		// For single stars, dm is the total mass loss.  We draw no
		// distinction between fast and slow evolution.

		// For ex-binary components, dm is the fast mass loss only.
		// (Slow mass loss was completely accounted for in the
		// preceding loop.)  Check for a kepler to avoid applying
		// rounding error corrections to binary component masses.
		// (Note that we now delete the kepler structure in the
		// case of fast mass loss, so the second check should
		// always return true.)

		if (abs(dm) > MINIMUM_MASS_LOSS && bi->get_kepler() == NULL) {
		
		    correct_leaf_for_change_of_mass(bi, dm);

		    if (b->get_kira_diag()->report_stellar_evolution &&
			b->get_kira_diag()->report_stellar_mass_loss &&
			abs(dm) >= MINIMUM_REPORT_MASS_LOSS) {

			cerr << "Single star " << bi->format_label()
			     << " lost mass:  dm = " << dm << endl;
			if (bi->is_low_level_node()) {
			    cerr << "star is leaf" << endl;
			    // pp3(bi->get_top_level_node(), cerr);
			}
		    }

		    if (bi->get_mass() < 0) {
			cerr << "\nEvolved star:  negative mass of star "
			     << bi->format_label()
			     << ", dm = " << dm << endl;
			// pp3(bi, cerr);
			hdyn * btop = bi->get_top_level_node();
			// pp3(btop, cerr);
			((star*)btop->get_starbase())->dump(cerr);
		    }
		}
	
		// Add kick velocity, if any...
	
		if (dv_sq > 0) {

		    b->get_kira_counters()->total_kick++;

		    if (full_dump) {

			// Print out the state of bi after the supernova.

			put_node(cout, *bt, false, 1);
		    }

		    correct_leaf_for_change_of_vector(bi, dv, &hdyn::get_vel,
							      &hdyn::inc_vel);

		    vector vnew = hdyn_something_relative_to_root(bi,
							      &hdyn::get_vel);
		    real dde_kick = 0.5 * bi->get_mass()
					* (square(vnew) - square(vnew-dv));

		    de_kick += dde_kick;

		    cerr << endl;
		    PRC(bi->format_label()), PRL(dde_kick);
		    PRL(dv);
		    pp3(bi->get_top_level_node(), cerr);
		}

	    }	// end if (bi->is_valid()) { ...

	}	// end for_all_leaves(hdyn, b, bi) { ...


	// *Don't* reset center of mass; recompute accs and jerks
	// on top-level nodes (latter not needed?).

	// b->to_com();

	// calculate_acc_and_jerk_on_all_top_level_nodes(b);
	// predict_loworder_all(b, b->get_system_time());

	real epot, ekin, etot;
	calculate_energies_with_external(b, epot, ekin, etot);

	real de_total = etot - etot0;

	if (b->get_kira_diag()->report_stellar_evolution) {
	    cerr << "End of correction: "; PRL(de_total);
	    PRC(epot); PRC(ekin); PRL(etot);
	}

	// Update statistics on mass and energy change components:

	b->get_kira_counters()->dm_massloss += mass0 - total_mass(b);
	
	b->get_kira_counters()->de_total += de_total;
	b->get_kira_counters()->de_massloss += de_total - de_kick;
	b->get_kira_counters()->de_kick += de_kick;

	predict_loworder_all(b, b->get_system_time());

	if (b->get_kira_diag()->report_stellar_evolution)
	    cerr << endl << "initialize_system_phase2 called from evolve_stars"
		 << endl;

	// ALWAYS reinitialize the system following stellar evolution...

       	initialize_system_phase2(b, 5);		// default set_dt

       	b->set_mass(total_mass(b));

	// print_recalculated_energies(b);

    }	// end if (correct_dynamics) { ...

    // test_kepler(b);

    // pp3("(21,100021)");
    // pp3("(23,100023)");


//      if (b11 && b11->get_kepler()) {
//  	cerr << "*************************************" << endl;
//  	PRL(0);
//  	b11->get_kepler()->print_all(cerr);
//  	pp3(b11->get_parent());
//  	cerr << "*************************************" << endl;
//      }


    if (b->get_kira_diag()->report_stellar_evolution) {
	cerr << "\nEnd of evolve_stars: "; PRL(correct_dynamics);
	PRL(b->get_kira_counters()->de_kick);
	cerr << "===============\n\n";
    }

    return correct_dynamics;
}

