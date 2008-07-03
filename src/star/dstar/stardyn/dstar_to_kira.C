
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// dstar_to_kira.C 
//
// Set of dstar functions referenced by dstar_kira (and associated
// hdyn functions).
//
// Externally visible functions:
//
//	bool create_or_delete_binary
//	bool binary_is_merged
//	bool binary_evolution
//	void update_kepler_from_binary_evolution
//
// An interesting mixture of dyn.h and hdyn.h references!
//
// SPZ+SMcM July 1998
//
// Globally visible functions:
//
//	binary_is_merged		bool function
//	create_or_delete_binary		switch:  calls adddouble
//						 or delete_double
//					called from:
//					    binary_evolution (this file)
//					    integrate_unperturbed_motion
//							(hdyn_unpert.C)
//					    dissociate_binary
//							(kira_stellar.C)
//	binary_evolution		evolve a binary node
//					called from evolve_stars
//							(kira_stellar.C)
//	update_kepler_from_binary_evolution	what it says...
//
// Local functions:
//
//	need_dstar
//	binary_wants_to_be_updated
//	evolve_an_unperturbed_binary	called by delete_double and
//					binary_evolution
//	initial_binary_period
//	delete_double
//	delete_double			delete double_star part of a node
//					called by create_or_delete_binary
//
// Functions delete_double() and binary_evolution() use the function
// hdyn::merge_nodes().  This may be OK because the only way to get to
// the reference in this file is via a call that originated in kira.

#include "dstar_to_kira.h"

#define REPORT_ADD_DOUBLE        false
#define REPORT_DEL_DOUBLE        false
#define REPORT_EVOLVE_DOUBLE     false

local bool need_dstar(hdyn* bi)
{
    if (bi->is_low_level_node() && (bi->get_elder_sister() == NULL)
	&& (bi->get_kepler() != NULL) && bi->get_fully_unperturbed()) {
	return true;
    } else {
	return false;
    }
}

// safety function, but not used in practice.
bool other_urgent_binary_matters(hdyn *bi)  		// bi is binary CM
{
    bool update_binary = false;
     
    starbase *cm = bi->get_starbase();
    starbase *ps = ((star*)cm)->get_primary();
    starbase *ss = ((star*)cm)->get_secondary();

    real pericenter = cm->get_semi()*(1-cm->get_eccentricity());

    real rp = ps->get_effective_radius();
    real rs = ss->get_effective_radius();

    if(pericenter <= rp * cnsts.parameters(tidal_circularization_radius) ||
       pericenter <= rs * cnsts.parameters(tidal_circularization_radius)) {

      update_binary = true;
    }
   
    return update_binary;
}


local bool binary_wants_to_be_updated(hdyn *b)		// b is binary CM
{
    bool update_binary = false;
     
    hdyn* od = b->get_oldest_daughter();
    real stellar_time = b->get_starbase()
			 ->conv_t_dyn_to_star(od->get_time());
    				      // 3/99 was get_system_time()
    real current_age  = b->get_starbase()->get_current_time();
    real dt_max       = b->get_starbase()->get_evolve_timestep();

    //if (stellar_time >= current_age + dt_max)
    if (stellar_time >= current_age 
                      + cnsts.star_to_dyn(binary_update_time_fraction)*dt_max)
        update_binary = true;

//    (SPZ: 21 Mar 2001) safety switched off.
//    if(other_urgent_binary_matters(b))
//      update_binary = true;

    return update_binary;
}

local void evolve_an_unperturbed_binary(hdyn *bi,	// bi is binary CM
					bool& check_binary_mass,
					bool& check_binary_sma,
					bool& check_binary_merger,
					bool force_evolve=false)
{ 
    check_binary_mass = check_binary_sma = check_binary_merger = false;

    real old_dyn_mass, old_dyn_sma;
    real new_dyn_mass_from_star, new_dyn_sma_from_star;

    hdyn* od = bi->get_oldest_daughter();
    real stellar_evolution_time = bi->get_starbase()
                                    ->conv_t_dyn_to_star(od->get_time());
    						 // 3/99 was get_system_time()

    kepler * p_kep = od->get_kepler();
    old_dyn_mass = bi->get_mass();
    old_dyn_sma  = p_kep->get_semi_major_axis();

    if (REPORT_EVOLVE_DOUBLE) {
	cerr << "Before evolution:"
	     << " sma= " << bi->get_starbase()->get_semi()
	     << " ecc= " << bi->get_starbase()->get_eccentricity() 
	     << " Time = "<<stellar_evolution_time <<"      ";
	put_state(make_state(dynamic_cast(double_star*, bi->get_starbase())));
	cerr << endl;
	((double_star*) (bi->get_starbase()))->dump(cerr);
	cerr << endl;
    }

    if (!force_evolve) {
	if (binary_wants_to_be_updated(bi)) 
	    bi->get_starbase()->evolve_element(stellar_evolution_time);
	else
	    ((double_star*)(bi->get_starbase()))->try_zero_timestep();
    }
    else {
      // **** Look also in double_star::evolve_the_binary() ***
      double_star *the_binary = dynamic_cast(double_star*, bi->get_starbase());
      the_binary->evolve_the_binary(stellar_evolution_time);
      // SPZ 8 Febr 2003.
      // this line: 
      // bi->get_starbase()->evolve_element(stellar_evolution_time);
      // was there instead of the current line:
      // bi->get_starbase()->evolve_the_binary(stellar_evolution_time);
      // which is commented out! 
      // The reasoning is as follows:
      // 1) for isolated binary evolution we would prefer to have the
      //    binary updated only when either of the stars require an upate.
      //    There is basically no external influence on the binary.
      // 2) However, when the binary can be perturbed by external means 
      //    (via the parameter force_update) the binary SHOULD be 
      //    updated irrespective of any other parameters or circumstances.
    }
      
    new_dyn_mass_from_star = get_total_mass(bi);
    new_dyn_sma_from_star = bi->get_starbase()->
			    conv_r_star_to_dyn(bi->get_starbase()->get_semi());

    // Relative position and velocity should be changed also
    // for the upper parent nodes!!!

    if (abs(new_dyn_mass_from_star-old_dyn_mass)/old_dyn_mass
	>=cnsts.star_to_dyn(stellar_mass_update_limit)) {  // MASS_UPDATE_LIMIT
	check_binary_mass=true;
    }
    if (abs(new_dyn_sma_from_star-old_dyn_sma)/old_dyn_sma
	>=cnsts.star_to_dyn(semi_major_axis_update_limit)) { //SMA_UPDATE_LIMIT
	check_binary_sma=true;
    }
    if (binary_is_merged(od)) {
	cerr << endl << "evolve_an_unperturbed_binary:  binary "
	     << bi->format_label()
	     << " is apparently merged" << endl;
	check_binary_merger = true;
    }
      
    star * primary = ((star*) (bi->get_starbase()))->get_primary();
    star * secondary = ((star*) (bi->get_starbase()))->get_secondary();
    real m1 = primary->get_total_mass();
    real m2 = secondary->get_total_mass();
    real semi = bi->get_starbase()->get_semi();
    real ecc = bi->get_starbase()->get_eccentricity();
    const char * t1 = type_string(primary->get_element_type());
    const char * t2 = type_string(secondary->get_element_type());

    if (REPORT_EVOLVE_DOUBLE) 
	if (check_binary_sma || check_binary_mass || check_binary_merger) {
	    cerr << " check: "
		 <<check_binary_mass<<" "
		 <<check_binary_sma<<" "
		 <<check_binary_merger<<endl;
	    bi->pretty_print_node(cerr); PRC(bi->get_time());
	    PRL(od->get_time());
	    PRC(bi->get_system_time()); PRL(stellar_evolution_time);
	    PRC(t1); PRC(m1); PRC(t2); PRC(m2); PRC(semi); PRL(ecc);
	    cerr << "After evolution:"
		 << " sma= " << semi
		 << " ecc= " << ecc 
		 << " Time = "<<stellar_evolution_time <<"      ";
	    put_state(make_state(dynamic_cast(double_star*,
					      bi->get_starbase())));
	    cerr << endl;
	    ((double_star*) (bi->get_starbase()))->dump(cerr);
	    cerr << endl;
	}

    //bi->get_starbase()->get_seba_counters()->step_dstar++;
}

local real initial_binary_period(double_star * ds)
{
    double_init *di = ds->get_initial_conditions();

    real pdays = 2*PI*sqrt(pow(di->semi*cnsts.parameters(solar_radius), 3.)
			   / (cnsts.physics(G)*di->mass_prim
			      *(1+di->q)*cnsts.parameters(solar_mass)))
			   /  cnsts.physics(days);
			   return pdays;
}

// Helper version of delete_double -- overloaded!

local void delete_double(hdyn *b)			// b is binary CM
{
    if (REPORT_DEL_DOUBLE) {
	cerr << "Unsafe double_star deletion: "<<endl;
	b->pretty_print_node(cerr); cerr << endl;
	((star*)(b->get_starbase()))->dump(cerr);
    }
    double_star *the_binary = dynamic_cast(double_star*, b->get_starbase());

    // First make sure that the binary is evolved to exact time.

    // Individual time of internal binary motion should be used here,
    // NOT the system time!!!

    hdyn* od = b->get_oldest_daughter();
    real stellar_time = b->get_starbase()
			 ->conv_t_dyn_to_star(od->get_time());
    				      // 3/99 was get_system_time()
    the_binary->evolve_the_binary(stellar_time);
    the_binary->set_node(NULL);

    story* old_story = the_binary->get_star_story();
    the_binary->set_star_story(NULL);

    if (old_story)
	delete old_story;
      
    // (Keep the story if we want to add a story chapter about the
    // binary evolution.)

    // Now safely delete the binary.

    delete the_binary;
      
    // Create a new starbase to the dyn.  See note in other delete_double.

    // b->set_starbase(new starbase(b));

    b->set_starbase(new starbase());
    b->get_starbase()->set_node(b);
}

//-----------------------------------------------------------------------------
//
// delete_double  -- the real "delete double()"
//
// Delete the double_star attached to the parent of the leaves.
// Does not check whether or not a double_star indeed exists.
//
// May set two "update_dynamics" flags:
//		0:	binary mass or semimajor axis changed during
//			the final binary evolution step
//		1:	binary orbit has changed significantly since
//			the dynamical binary became unperturbed
//
// Called only by create_or_delete_binary.
//
//-----------------------------------------------------------------------------

local void delete_double(hdyn *b,		    // b is binary CM
			 bool * update_dynamics,
			 bool allow_evolution_before_termination,
			 int full_dump = 0)
{
    if (b->get_oldest_daughter()->get_kepler() == NULL) {
	cerr << "delete_double: deleting double for "
	     << b->format_label() << " without a kepler pointer"
	     << endl;
	delete_double(b);
	return;
    }

    // cerr << endl << "in delete_double" << endl;
  
    bool kepler_created = false;

    hdyn* od = b->get_oldest_daughter();
    if (!od->get_kepler()) {

	// This cannot occur -- already tested above...

	cerr << "No Kepler in delete_double; create one."<<endl;
	new_child_kepler(b, od->get_time(), MAX_CIRC_ECC); 
	kepler_created = true;     // 3/99 was get_system_time()
    }
    else {
	od->get_kepler()->
	    set_circular_binary_limit(MAX_CIRC_ECC);	// Probably unnecessary
	dyn_to_child_kepler(b, od->get_time()); // 3/99 was get_system_time()
    }

    if (REPORT_DEL_DOUBLE) {
	cerr << "kepler created in delete_double ";
	cerr << "at address " << od->get_kepler() 
	     << endl;
    }
  
    kepler* johannes = od->get_kepler();
  
    if (REPORT_DEL_DOUBLE) {
	cerr << "delete double star "; b->pretty_print_node(cerr);
	cerr<<endl;
	((double_star*)(b->get_starbase()))->dump(cerr);
    }
    
    double_star *the_binary = dynamic_cast(double_star*, b->get_starbase());

    if (REPORT_DEL_DOUBLE) 
	put_state(make_state(the_binary));

    // First make sure that the binary is evolved to exact time.

    real stellar_time = b->get_starbase()
			 ->conv_t_dyn_to_star(od->get_time());
    				      // 3/99 was get_system_time()

    bool change_mass = false, change_sma = false, check_merger = false;
    bool force_evolve = true;

    // Optionally evolve the binary (to its current dynamical time)
    // before deleting it.  If no evolution occurs, all dynamics
    // flags will be false on return.

    if (allow_evolution_before_termination) {

	evolve_an_unperturbed_binary(b,
				     change_mass,
				     change_sma,
				     check_merger,
				     force_evolve);
    }

    update_dynamics[0] |= (change_mass || change_sma);

    if (REPORT_DEL_DOUBLE) {
	PRL(update_dynamics[0]);
	the_binary->dump(cerr);
	put_state(make_state(the_binary));
    }

    if (check_merger) {

	hdyn* bcoll = od->get_binary_sister();

	cerr <<"delete_double:  merger within "
	     << b->format_label()
	     << " triggered by ";
	cerr << bcoll->format_label();
	int p = cerr.precision(INT_PRECISION);
	cerr << " at time (tsys)" << b->get_system_time()
	     << "(binary: "<< od->get_time() <<")"<< endl; 
	cerr.precision(p);

	cerr << "merging nodes(1)..." << endl << flush;

	od->merge_nodes(bcoll, full_dump);
	update_dynamics[0] = true;

	cerr << "...done" << endl << flush;
      
	// Components are not deleted in merge_nodes.  Do this now.
	// Also, if unperturbed binaries are resolved in perturber lists,
	// must check and correct the lists now.

	if (RESOLVE_UNPERTURBED_PERTURBERS) {
	    hdynptr del[2];
	    del[0] = od;
	    del[1] = bcoll;
	    correct_perturber_lists(b->get_root(), del, 2, b);
	}

	if (kepler_created) delete johannes;

	delete od;
	delete bcoll;

	delete the_binary;

	return; 
    }
    else {

	// Unperturbed motion is about to be terminated.  If the final
	// orbit is radically different from the initial one, we will
	// have to recompute the time step on returning to kira...

	if (the_binary->get_period()
	    / initial_binary_period(the_binary) < 0.9) {
	    real dm_fast = sudden_mass_loss(b); 
	    update_kepler_from_binary_evolution(b, dm_fast);
	    update_dynamics[1] = true;
	    if (REPORT_DEL_DOUBLE) 
		PRL(update_dynamics[1]);
	}

	// Could keep the story, possibly to be added to the root or
	// component stories, or we might add a story chapter about
	// the binary evolution.  However, for now, any story text in
	// the parent node is simply discarded.

	story* old_story = the_binary->get_star_story();
	the_binary->set_star_story(NULL);

	if (old_story)
	    delete old_story;
      
//       cerr << "before deleting the_binary" << endl;
//       PRL(b);
//       PRL(b->get_starbase());
//       PRL(b->get_starbase()->get_star_story());

	the_binary->set_node(NULL);
	delete the_binary;
      
	// Create a new starbase to the dyn.

	// Note from Steve to Simon, 8/27/98.  Deleting the_binary here
	// has the effect of corrupting the star_story pointer associated
	// with b (see commented lines before and after the deletion).
	// However, "new starbase(b)" uses any existing (i.e. non-null)
	// star_story pointer it finds, so it propogates the corruption into
	// the newly created starbase.  This can have disastrous effects
	// much later when we come to print or delete b!
	//
	// (It took me a day to find this bug!)

//       cerr << "after deleting the_binary" << endl;
//       PRL(b);
//       PRL(b->get_starbase());
//       PRL(b->get_starbase()->get_star_story());

//	 b->set_starbase(new starbase(b));

	b->set_starbase(new starbase());
	b->get_starbase()->set_node(b);

	// This creates a new starbase with an empty star story.
	// Actually, adddouble sets star_story = NULL, so NULL should
	// be OK here too...

	// If we want to preserve the old_story, don't delete above,
	// and do this:
	// 
	//	b->get_starbase()->set_star_story(old_story);
    }
      
    if (kepler_created) {
	if (od->get_kepler()) {
	    od->set_kepler(NULL);
	    od->get_binary_sister()->set_kepler(NULL);
	    delete johannes;
	}
    }
}



//------------------------------------------------------------------------
//
// Global functions:
//
//------------------------------------------------------------------------

bool create_or_delete_binary(hdyn *bi,		      // pointer to parent node
			     bool * update_dynamics,
			     int full_dump /*= 0*/ )  // default = false
{
// Two "update_dynamics" flags may be modified (by delete_double):
//
//		0:	binary mass or semimajor axis changed during
//			the final binary evolution step
//		1:	binary orbit has changed significantly since
//			the dynamical binary became unperturbed
//
// They are *not* initialized by this function.

    if (REPORT_ADD_DOUBLE || REPORT_DEL_DOUBLE )
	cerr << "create_or_delete_binary" << endl;
  
    if (has_dstar(bi->get_oldest_daughter())) {
	if (REPORT_DEL_DOUBLE )
	    cerr << "create_or_delete_binary: delete existing binary from "
		 << bi->format_label()<<endl;

	// Setting this flag false means that the binary will not
	// be allowed to evolve before its dstar is deleted.  In that
	// case, update_dynamics will be set false on return.

	bool allow_evolution_before_termination = false;

	delete_double(bi, update_dynamics,
		      allow_evolution_before_termination,
		      full_dump);

	//bi->get_starbase()->get_seba_counters()->del_dstar++;

	return true;
    }
    else {
	if (REPORT_ADD_DOUBLE)
	    cerr << "create_or_delete_binary: create new binary to "
		 << bi->format_label()<<endl;

	adddouble(bi, bi->get_oldest_daughter()->get_time());
				     // 3/99 was get_system_time()

	//bi->get_starbase()->get_seba_counters()->add_dstar++;

	return true;
    }

    return false;
}

bool binary_is_merged(dyn* bi)
{
    starbase * p_star = bi->get_parent()->get_starbase();
    if (bi->get_kepler()
	&& p_star->get_element_type()==Double) {      
	if (p_star->get_bin_type()==Merged)
	    return true;
    }
    return false;
}

bool binary_evolution(hdyn *b,		// root node
		      int full_dump)	// default = 0

// Return value is true iff a system reinitialization will be needed
// after returning to kira.

{
    bool update_dynamics[2] = {false, false};
    		
    // int p = cerr.precision(HIGH_PRECISION);
    // cerr << "entering binary_evolution at time "
    //      << b->get_system_time() << endl;
    // cerr.precision(p);

    for_all_leaves(hdyn, b, bi) {

    	if (need_dstar(bi) && (! has_dstar(bi))) {

	    hdyn *  parent = bi->get_parent();

	    create_or_delete_binary(bi->get_parent(),
				    update_dynamics, full_dump);

	    if (REPORT_ADD_DOUBLE) {
	        cerr << "new double star created at time "
		     << bi->get_time()<<endl;
	        pp3(parent, cerr);
	    }
	}

	if (has_dstar(bi)) {			    // bi is binary component

	    if (need_dstar(bi)) {

		if (binary_wants_to_be_updated(bi->get_parent())) {

		    if (REPORT_EVOLVE_DOUBLE) 
			cerr << "\nEvolving binary " << bi->format_label()
			     << endl;

		    // The procedure below is effectively similar to that
		    // followed in delete_double.

		    bool change_mass, change_sma, check_merger;
		    bool force_evolve = false;
		
		    evolve_an_unperturbed_binary(bi->get_parent(),
						 change_mass,
						 change_sma,
						 check_merger,
						 force_evolve);
		  
		    update_dynamics[0] |= (change_mass || change_sma);

		    if (REPORT_DEL_DOUBLE)
			PRL(update_dynamics[0]);

		    if (check_merger) {

			hdyn* par = bi->get_parent();
			hdyn* bcoll = bi->get_binary_sister();

			cerr << "binary_evolution:  merger within "
			     << par->format_label()
			     << " triggered by ";
			cerr << bcoll->format_label();
			int p = cerr.precision(INT_PRECISION);
			cerr << " at time " << b->get_time() << endl;
			cerr.precision(p);

// 			cerr << endl << "computing energies before merger:"
// 			     << endl << flush;
// 			real epot0, ekin0, etot0;
// 			calculate_energies_with_external(b, epot0,
// 							 ekin0, etot0);
// 			cerr << "OK" << endl << endl << flush;
			cerr << "merging nodes(2)..." << endl;

			bi->merge_nodes(bcoll, full_dump);

			cerr << "...done" << endl << flush;
// 			cerr << endl << "computing energies after merger (#1):"
// 			     << endl << flush;
// 			calculate_energies_with_external(b, epot0,
// 							 ekin0, etot0);
// 			cerr << "OK" << endl << endl << flush;

			// Components are not deleted by merge_nodes.  Do this
			// here, and adjust perturber lists if necessary.
		    
			if (RESOLVE_UNPERTURBED_PERTURBERS) {
//			if (true) {
			    hdynptr del[2];
			    del[0] = bi;
			    del[1] = bcoll;
			    cerr << "correcting perturber lists"
			 	 << endl << flush;
			    correct_perturber_lists(par->get_root(),
						    del, 2, par);
			}

			delete bi;
			delete bcoll;

// 			cerr << endl << "computing energies after merger (#2):"
// 			     << endl << flush;
// 			calculate_energies_with_external(b, epot0,
// 							 ekin0, etot0);
// 			cerr << "OK" << endl << endl << flush;

			// Looks as though everything is now taken care of
			// in case of a merger.  The binary is replaced by
			// its center of mass node and all stellar information
			// is updated.  Dynamics will be corrected on return.
			// Should be safe to continue the tree traversal at
			// the parent node...

			// *Don't* return immediately.  Set update_dynamics
			// true because binary evolution has changed the tree.

			// return true;

			update_dynamics[0] = true;
			bi = par;
		    }
		}

	    } else {
		
		if (REPORT_DEL_DOUBLE) {
		    cerr << "delete double star  at time "
			 << bi->get_time() << endl;
		    pp3(bi->get_parent(), cerr);
		    PRL(bi->is_low_level_node());PRL(bi->get_elder_sister());
		    PRL(bi->get_kepler());PRL(bi->get_fully_unperturbed());
		    PRL(need_dstar(bi));PRL(has_dstar(bi));
		}

		// This is not very safe.
		// double_star needs to be updated before kepler is deleted.
		// Do it now in delete_double, not very pretty.
		// (SPZ:8/02/1998)
		// delete_double(bi->get_parent());

		create_or_delete_binary(bi->get_parent(),
					update_dynamics, full_dump);
	    }
	}
    }

    return update_dynamics[0];
}

void update_kepler_from_binary_evolution(hdyn* the_binary_node, 
					 real dm_fast)
{
    if (REPORT_EVOLVE_DOUBLE)
	((double_star*)(the_binary_node->get_starbase()))->dump(cerr);

    kepler *johannes = the_binary_node->get_oldest_daughter()->get_kepler();
    real sma = the_binary_node->get_starbase()->conv_r_star_to_dyn(
                the_binary_node->get_starbase()->get_semi());
    real ecc = the_binary_node->get_starbase()->get_eccentricity();

    real m_total = get_total_mass(the_binary_node) - dm_fast;

    real mean_anomaly = johannes->get_mean_anomaly();

    // real peri = johannes->get_periastron();
    real peri = sma * (1-ecc);
    hdyn * b = the_binary_node->get_oldest_daughter();
    johannes->transform_to_time(b->get_time());

    real time = johannes->get_time();
    if (REPORT_EVOLVE_DOUBLE) {
	PRC(time); PRC(sma); PRL(m_total); PRC(mean_anomaly); PRL(peri);
    }

    make_standard_kepler(*johannes, time, m_total, -0.5 * m_total / sma, ecc,
			 peri, mean_anomaly, 0);

    // Note from Steve, 10/9/98: make_standard_kepler creates a new
    // kepler with the specified properties.  However, the default
    // action is to aligh the orbit with the coordinate axes, placing
    // the binary in the x-y plane and losing orientation information.
    // The trailing "0" here retains the old orientation vectors.

    // Feb 26, 1998.
    // I believe pred_advance must be called here to
    // set the periastron time

    johannes->pred_advance_to_periastron();

    if (johannes->get_pred_time() < b->get_time()) {
	cerr << "update_kepler_from_binary, time went backwards"<<endl;
	PRC(the_binary_node->format_label());  
	PRC(johannes->get_pred_time()); PRL(b->get_time());
    }
}
