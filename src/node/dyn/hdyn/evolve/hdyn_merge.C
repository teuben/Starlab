
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ 
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ 
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ 
 //                                                       //            _\|/_
//=======================================================//              /|\ 

//
//  hdyn_merge.C: functions related to physical mergers in kira
//		  (formerly part of hdyn_tree.C).
//.............................................................................
//    version 1:  Jun 2000	Steve McMillan
//.............................................................................
//
//  Externally visible functions:
//
//      void   hdyn::merge_logs_after_collision
//	hdyn*  hdyn::merge_nodes
//	bool   hdyn::check_merge_node
//
//.............................................................................

#include "hdyn.h"
#include <star/dstar_to_kira.h>
#include <star/single_star.h>

#define CALCULATE_POST_COLLISION_ON_GRAPE

// Local function (incomplete).

local hdyn* apply_tidal_dissipation(hdyn* bi, hdyn* bj, kepler k)
{
    // Change the relative orbit of bi and bj to take tidal dissipation
    // into account.  Kepler structure k describes the relative motion
    // of bi and bj on entry.  The system has just passed periastron where
    // tidal dissipation is to be applied.

    print_encounter_elements(bi, bj, "Tidal encounter");


    //-----------------------------------------------------------------


    cerr << endl << "**** tidal dissipation suppressed ****"
	 << endl << endl;
    return NULL;


    //-----------------------------------------------------------------


    real t_peri = k.return_to_periastron();
    if (t_peri > bi->get_time())
	warning("apply_tidal_dissipation: stars not past periastron");

    real m_tot = bi->get_mass() + bj->get_mass();
    vector r_rel = k.get_rel_pos();
    vector v_rel = k.get_rel_vel();

    real de_diss = 0;
    for_all_daughters(hdyn, bi->get_parent(), bb) {
	if (bb->get_radius() > 0) {
	    real dde_diss = 0;
	    real rr = bb->get_radius()/k.get_separation();
	    real eta_pt = sqrt(bb->get_mass() / m_tot) * pow(rr, -1.5);
	    stellar_type type = bb->get_starbase()->get_element_type();
	    dde_diss = tf2_energy_diss(eta_pt, type)
			 + pow(rr, 2) * tf3_energy_diss(eta_pt, type);

	    cerr << endl;
	    PRC(type); PRC(bb->get_radius()), PRL(k.get_separation());
	    PRC(rr); PRC(eta_pt); PRL(dde_diss);

	    de_diss += pow(bb->get_binary_sister()->get_mass(), 2)
			* pow(rr, 6) * dde_diss / bb->get_radius();
	}
    }

    real m_red = bi->get_mass()*bj->get_mass() / m_tot;
    real kinetic = 0.5*m_red*square(v_rel);
    real potential = -m_red*m_tot/abs(r_rel);

    cerr << endl << "at periastron:  U = " << potential
	 << "  K = " << kinetic
	 << "  total = " << potential + kinetic << endl;
    cerr << "de_diss = " << de_diss << endl;

    // Reduce energy by de_diss, keep angular momentum constant,
    // and compute the new orbit.

    real new_E = k.get_energy() - de_diss;
    real J = k.get_angular_momentum();

    real new_sma = -0.5*m_tot/new_E;	// Note: sma may have either sign
    					// 	 -- not kepler convention

    real new_ecc_2 = 1 - J*J/(m_tot*new_sma);

    if (new_ecc_2 < 0) {
	cerr << "immediate merger with J^2 > Ma" << endl;
	return bj;
    }

    cerr << "old orbit parameters:  sma = " << k.get_semi_major_axis()
	 << "  ecc = " << k.get_eccentricity()
	 << "  periastron = " << k.get_periastron() << endl;

    real new_ecc = sqrt(new_ecc_2);
    k.set_semi_major_axis(abs(new_sma));
    k.set_eccentricity(new_ecc);
    k.initialize_from_shape_and_phase();

    cerr << "expected orbit parameters:  sma = " << abs(new_sma)
	 << "  ecc = " << new_ecc << endl;
    cerr << "new orbit parameters:  sma = " << k.get_semi_major_axis()
	 << "  ecc = " << k.get_eccentricity()
	 << "  periastron = " << k.get_periastron() << endl;

    real old_phi = bi->get_mass()*bi->get_pot()
			+ bj->get_mass()*bj->get_pot()
			    - 2*potential;

    k.transform_to_time(bi->get_time());

    r_rel = k.get_rel_pos();
    v_rel = k.get_rel_vel();

    vector r_com = (bi->get_mass()*bi->get_pos()
		     + bj->get_mass()*bj->get_pos()) / m_tot;
    vector v_com = (bi->get_mass()*bi->get_vel()
		     + bj->get_mass()*bj->get_vel()) / m_tot;

    // Adjust positions and velocities of the interacting stars.

    bi->set_pos(r_com - bj->get_mass()*r_rel/m_tot);
    bj->set_pos(r_com + bi->get_mass()*r_rel/m_tot);
    bi->set_vel(v_com - bj->get_mass()*v_rel/m_tot);
    bj->set_vel(v_com + bi->get_mass()*v_rel/m_tot);

    // Recompute accs and jerks and update energy bookkeeping.

    hdynptr list[2] = {bi, bj};
    bool _false = false;
    calculate_acc_and_jerk_for_list(bi->get_root(), list, 2,
				    bi->get_time(),
				    _false, _false, _false, _false);

    bi->get_kira_counters()->de_total -= de_diss;
    bi->get_kira_counters()->de_tidal_diss -= de_diss;

    potential = -m_red*m_tot/abs(k.get_separation());
    real new_phi = bi->get_mass()*bi->get_pot()
			+ bj->get_mass()*bj->get_pot()
			    - 2*potential;

    PRC(old_phi); PRC(new_phi); PRL(potential);
    PRL(new_phi-old_phi);

    return NULL;
}

// Flag close encounter distances between stars.  Output occurs on any
// unbound encounter, and on the first encounter in a bound system.

local void print_perioapo_clustron(hdyn* bi) {

  hdyn *bj = bi->get_root();

  real d2cd_2 = -1;
  real d2cd_1 = -1;
  real d2cd   = -1;

    // Variables:	bi is star, bj is cluster com
    //			pa_time is time of previous coll
    //			d2cd is squared distance to cluster com
    //			d2cd_1 is previous d2cd
    //			d2cd_2 is previous d2cd_1
    //			pcp_time is time of last periclustron
    //			pca_time is time of last apoclustron
    //			pcp_cntr counts bound periclustron
    //			pca_cntr counts bound apoclustron
    //
    // All but bi and bj are stored in the bi log story.

#if 0
    cerr << "\nIn print_close_encounter with bi =  " << bi->format_label()
	 << "  at time " << bi->get_system_time() << endl;
    cerr << "     "; print_coll(bi,2);
    cerr << "     "; print_nn(bi,2); 
    cerr << "     (bj = " << bj->format_label();
    cerr << ":  coll = "; print_coll(bj);
    cerr << ",  nn = "; print_nn(bj); cerr << ")" << endl;

#endif
  
    // Retrieve coll information from the log story.

    d2cd_2 = getrq(bi->get_log_story(), "d2cd_1");
    d2cd_1 = getrq(bi->get_log_story(), "d2cd");
    d2cd   = square(bi->get_pos() - bj->get_pos());

//    PRC(bi->format_label());PRC(d2cd_2);PRC(d2cd_1);PRL(d2cd);

    if (d2cd_2 > 0) {

      if (d2cd_2 > d2cd_1 && d2cd >= d2cd_1) {    // just passed pericluctron

	int pcp_cntr = 0;
	real E = get_total_energy(bi, bj);

	if (E > 0) {

	  // Always print unbound encounter elements.

	  print_encounter_elements(bi, bj, "Perioclustron", false);

	}
	else {
	  if (find_qmatch(bi->get_log_story(), "pcp_cntr")) 
	    pcp_cntr = getiq(bi->get_log_story(), "pcp_cntr");

	  pcp_cntr++;
	  
	  if (pcp_cntr == 1) 
	    print_encounter_elements(bi, bj, "First periclustron", false);
	  else
	    print_encounter_elements(bi, bj, "Periclustron", false);
	}

	// Save data on last pericenter (bound or unbound).
	// Note that an unbound pericenter resets pcc_cntr to zero.

	putiq(bi->get_log_story(), "pcp_cntr", pcp_cntr);
	putrq(bi->get_log_story(), "pcp_time", bi->get_time());  
	
      }
      else if (d2cd_2 < d2cd_1 && d2cd <= d2cd_1) { //just passed apocluctron

	int pca_cntr = 0;
	real E = get_total_energy(bi, bj);

	if (E > 0) {

	  // Always print unbound encounter elements.
	  
	  print_encounter_elements(bi, bj, "Apoclustron", false);

	}
	else {
	  if (find_qmatch(bi->get_log_story(), "pca_cntr")) 
	    pca_cntr = getiq(bi->get_log_story(), "pca_cntr");

	  pca_cntr++;

	  if (pca_cntr == 1) 
	    print_encounter_elements(bi, bj, "First apoclustron", false);
	  else
	    print_encounter_elements(bi, bj, "Apoclustron", false);
	}

	// Save data on last pericenter (bound or unbound).
	// Note that an unbound pericenter resets pca_cntr to zero.

	putiq(bi->get_log_story(), "pca_cntr", pca_cntr);
	putrq(bi->get_log_story(), "pca_time", bi->get_time());  
      }
      else {
	  
	// Unlikely multiple encounter(?).  Has been known to occur
	// when one incoming star overtakes another.

	//	  cerr << endl << "print_close_encounter:  "
	//	       << "bi = " << bi->format_label()
	//	       << " at time " << bi->get_system_time() << endl;
	//	  cerr << "     current coll = " << bj->format_label();
	//	  cerr << endl;
	//	  cerr << "  (periclustron counter = "
	//	       << getiq(bi->get_log_story(), "pcp_cntr") << ")\n";
      }

    }

    putrq(bi->get_log_story(), "d2cd_1", d2cd_1);  
    putrq(bi->get_log_story(), "d2cd", d2cd);  
    //    putrq(bi->get_log_story(), "pcd_time", bi->get_time());  

#if 0
    cerr << "\nAt end of print_close_encounter at time "
	 << bi->get_system_time() << endl;
    cerr << "bi = " << bi->format_label();
#endif

}

hdyn* hdyn::check_periapo_node() {

  // for now only for single stars.
  //  if (is_leaf() && (is_top_level_node() || parent->is_top_level_node())

    print_perioapo_clustron(this);

}


hdyn* hdyn::check_merge_node()
{
    // Note: merger criterion is based on coll pointer.

    // In the case of a collision we simply return a pointer to the
    // colliding star and take no other action.  Otherwise, we do
    // bookkeeping, apply tidal effects, etc. in place.

    // In constant-density case, collision criterion also covers
    // tidal destruction by a massive black hole so long as the
    // "radius" of the hole (based on constant density) is much
    // greater than the radius of a particle.
  
    if (use_dstar) {
	if (is_low_level_node()) {
	    if (binary_is_merged(this)) {
		cerr << "check_merge_node: merging a binary with logs:"
		     << endl;
		print_log_story(cerr);
		cerr << "and\n";
		get_binary_sister()->print_log_story(cerr);
		return get_binary_sister();
	    }
	}
    }

    // Note multiple functionality: check for and implement stellar
    // mergers; implement tidal dissipation at periastron.

    // Check for encounters moved to kira.C by Steve, 3/24/00.

    if (stellar_capture_criterion_sq <= 0) return NULL;

    hdyn *coll = get_coll();

    if (is_leaf() && coll != NULL && coll->is_leaf() && (kep == NULL)) {

	// Check for tidal capture and physical stellar collision.

	if (!coll->is_valid())
	    warning("check_merge_node: invalid coll pointer");

	real sum_of_radii_sq = pow(get_sum_of_radii(this, coll), 2);


	// (NB: d_coll_sq = distance_squared - sum_of_radii_sq)
      
	if (get_d_coll_sq() <= (stellar_capture_criterion_sq-1)
						* sum_of_radii_sq) {

	    // Note that the pair must approach within
	    // stellar_capture_criterion in order for
	    // stellar_merger_criterion to be checked.

	    // cerr << endl << "Capture between " << format_label();
	    // cerr << " and " << coll->format_label()
	    //      << " at time = " << get_time() << endl;

	    // pp3_minimal(get_top_level_node(), cerr);
	    // if (get_top_level_node()->get_nn())
	    //     pp3_minimal(get_top_level_node()->get_nn(), cerr);
	
	    // For a 2-body orbit, use kepler pericenter rather than
	    // sampled orbit points to determine whether or not a
	    // collision will occur.  Probably desirable to have
	    // stellar_capture_criterion a few times greater than
	    // stellar_merger_criterion to ensure that we get close
	    // enough to the collision for a more precise determination
	    // of its elements.

	    // Don't use the kepler orbit if the two stars aren't sisters!
	    // Don't apply tidal dissipation if the two stars aren't sisters!

	    // Original:
	    //
	    // if (d_coll_sq > (stellar_merger_criterion_sq-1)
	    //					* sum_of_radii_sq) {

	    // BUG corrected by Steve 12/98:
	    //
	    // if (k.get_periastron() > (stellar_merger_criterion_sq-1)
	    //					* sum_of_radii_sq) {

	    real sep_sq = get_d_coll_sq();
	    kepler k;

	    if (get_parent() == coll->get_parent()) {
		coll->synchronize_node();
		initialize_kepler_from_dyn_pair(k, this, coll, false);
		sep_sq = min(sep_sq, k.get_periastron() * k.get_periastron()
						- sum_of_radii_sq);
	    }

	    // Check for collision.

	    if (sep_sq <= (stellar_merger_criterion_sq-1)
	    					* sum_of_radii_sq)
		return coll;

	    // Check for tidal dissipation.

	    if (get_parent() != coll->get_parent()) return NULL;

	    if (find_qmatch(get_log_story(), "tidal_dissipation")) {
		rmq(get_log_story(), "tidal_dissipation");
		return apply_tidal_dissipation(this, coll, k);
	    }

	    real t_peri = k.get_time_of_periastron_passage();
	    if (time < t_peri && time + timestep >= t_peri)
		putiq(get_log_story(), "tidal_dissipation", 1);

	}
    }

    return NULL;
}


// correct_nn_pointers:  Check nn and coll pointers of all nodes in the
//			 system, replacing them by cm if they lie on
//			 the specified list.
//
// These checks should probably be applied every time correct_perturber_lists
// is invoked...

local void correct_nn_pointers(hdyn * b, hdyn ** list, int n, hdyn * cm)
{
    for_all_nodes(hdyn, b, bb) {
	for (int i = 0; i < n; i++) {
	    if (bb->get_nn() == list[i]) bb->set_nn(cm);
	    if (bb->get_coll() == list[i]) bb->set_coll(cm);
	}
    }
}

void hdyn::merge_logs_after_collision(hdyn *bi, hdyn* bj) {

  // log only the moment of the last merger
  putrq(get_log_story(), "last_merger_time", bi->get_time());
  
  if(find_qmatch(bi->get_log_story(), "black_hole"))
    putiq(get_log_story(), "black_hole",
	  getiq(bi->get_log_story(), "black_hole"));
  else if(find_qmatch(bj->get_log_story(), "black_hole"))
    putiq(get_log_story(), "black_hole",
	  getiq(bj->get_log_story(), "black_hole"));
}


#define CONSTANT_DENSITY	1	// Should be a parameter...
#define MASS_LOSS		0.0
#define KICK_VELOCITY		0.0

hdyn* hdyn::merge_nodes(hdyn * bcoll)
{

    // This function actually does the work of merging nodes.
    // It returns a pointer to the newly merged CM.

    if (diag->tree) {
	if (diag->tree_level > 0) {
	    cerr << endl << "merge_nodes: merging leaves ";
	    pretty_print_node(cerr);
	    cerr << " and ";
	    bcoll->pretty_print_node(cerr);

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << " at time " << system_time << endl;
	    cerr.precision(p);

	    if (diag->tree_level > 1) {
		cerr << endl;
		pp3(this, cerr);
		pp3(bcoll, cerr);
		if (diag->tree_level > 2) {
		    cerr << endl;
		    put_node(cerr, *this,  options->print_xreal);
		    put_node(cerr, *bcoll, options->print_xreal);
		}
		cerr << endl;
	    }
	}
    }

    bool is_pert = (get_kepler() == NULL);
    real sum_of_radii_sq = pow(get_sum_of_radii(this, bcoll), 2);
//    real sum_of_radii_sq = pow(radius+bcoll->radius, 2);

    cerr << "merge_nodes: calling synchronize_tree..." << flush;
    PRL(cpu_time());
    synchronize_tree(get_root());
    cerr << "back" << endl << flush;
    PRL(cpu_time());

    if (is_top_level_node() || get_binary_sister() != bcoll) {

        // Nodes to be merged ('this' and bcoll) are not binary
	// sisters.  Force them to become sisters prior to actually
	// merging them, temporarily placing them in the same
	// subtree if necessary.

	cerr << "creating CM node" << endl;

	hdyn* ancestor = common_ancestor(this, bcoll);
	
	if (ancestor->is_root()) {

	    // Nodes are not in the same subtree.  Combine their
	    // subtrees into a single clump to allow use of function
	    // move_node below.  Undo the combination before 
	    // continuing, if necessary.
	    
	    bool decombine = !is_top_level_node()
				&& !bcoll->is_top_level_node();

	    cerr << "merge_nodes: call combine_top_level_nodes" << endl;

	    // PRL(this);
	    // PRL(get_top_level_node());
	    // PRL(bcoll);
	    // PRL(bcoll->get_top_level_node());

	    combine_top_level_nodes(get_top_level_node(),
				    bcoll->get_top_level_node());

	    // pp2(get_top_level_node(), cerr);

	    if (get_binary_sister() != bcoll)
	        move_node(bcoll, this);		// Move bcoll to become
		  				// sister of this.
	    
	    if (decombine) {

		// cerr << "call split_top_level_node 3" << endl;

		split_top_level_node(get_top_level_node());

	    }
	    
	} else {

	    // pp2(get_top_level_node(), cerr);

	    move_node(bcoll, this);		// Move bcoll to become
						// sister of this.
	}		
    }
    
    // Nodes 'this' and bcoll are binary sisters (in all cases).
    // Merge them into their center-of-mass node and remove them.
    
    // pp2(get_top_level_node(), cerr);
    
    // Compute energies prior to merger, for bookkeeping purposes:

    cerr << "calculating energies..." << endl << flush;
    
    real epot0, ekin0, etot0;

    //calculate_energies(get_root(), eps2, epot0, ekin0, etot0);//dyn function
    // replaced by GRAPE friendly function (SPZ, March 2001)
    calculate_internal_energies(get_root(), epot0, ekin0, etot0);
    PRC(epot0); PRC(ekin0); PRL(etot0);

    hdyn* cm = get_parent();
    
    real epot_int = -mass * bcoll->mass / abs(pos - bcoll->pos);
    real ekin_int = 0.5 * (mass * square(vel)
			   + bcoll->mass * square(bcoll->vel));
    real etot_int = epot_int + ekin_int;
    
    PRC(epot_int); PRC(ekin_int); PRL(etot_int);
    PRL(cpu_time());
    //pp3(cm, cerr);
    //pp3(cm->get_root(), cerr);

    label_merger_node(cm);

    // Default mass-radius relation:

    cm->set_radius(radius + bcoll->radius);
    PRL(cm->format_label());

    // pp3(cm->get_root(), cerr);
    
    real dm;
    vector dv;

    if (!use_sstar) {

	// Test code for systems without stellar evolution.
	// Set mass loss and velocity kick.
	
	dm = -MASS_LOSS*cm->mass;
	
	real costh = randinter(-1, 1);
	real sinth = sqrt(1 - costh*costh);
	if (randinter(0, 1) < 0.5) sinth = -sinth;
	real phi = TWO_PI*randinter(0, 1);
	dv = KICK_VELOCITY * sqrt(cm->mass/cm->radius)
			   * vector(sinth*cos(phi),
				    sinth*sin(phi),
				    costh);

	print_encounter_elements(this, bcoll, "Collision");

	if (CONSTANT_DENSITY) {

	    // Override default M-R relation by requiring that the
	    // (average) density remain constant:
	    //
	    //	    M_cm / R_cm^3 = 0.5 * (M1 / R1^3 + M2 / R2^3)

	    real rho1 = 0, rho2 = 0;
	    if (radius > 0) rho1 = mass/pow(radius, 3);
	    if (bcoll->radius > 0) rho2 = bcoll->mass/pow(bcoll->radius, 3);

	    if (rho1+rho2 > 0)
		cm->set_radius(pow(2*(1-MASS_LOSS)*cm->mass
				    / (rho1 + rho2), 1.0/3));
	    else
		cm->set_radius(0);

	} else
	    cm->set_radius((1-MASS_LOSS)*cm->get_radius());

	cerr << "non-stellar merger: new radius = "
	     << cm->get_radius() << endl;

    } else {

	// Real code uses stellar evolution code for data to
	// initialize collision product.
	
	hdyn * primary = this;
	hdyn * secondary = bcoll;
	
	cerr << "merge_nodes: merging two STARS"<<endl;
	// pp2(cm->get_root(), cerr);

	((star*)(primary->sbase))->dump(cerr, false);
	((star*)(secondary->sbase))->dump(cerr, false);
	
	if (!merge_with_primary(dynamic_cast(star*, primary->sbase),
			        dynamic_cast(star*, secondary->sbase))) {
	    primary = bcoll;
	    secondary = this;
	}

	print_encounter_elements(primary, secondary, "Collision");
	
	((star*)(primary->sbase))
	                ->merge_elements(((star*)(secondary->sbase)));
	
	cerr << "Merger product: "<< endl;
	cerr << format_label() << " (";
	//	        put_state(make_star_state(sbase), cerr);
	        put_state(make_star_state(primary), cerr);
	cerr << "; M=" << get_total_mass(primary) << " [Msun], "
	     << " R=" << get_effective_radius(primary) << " [Rsun]). " << endl;
	((star*)(primary->sbase))->dump(cerr, false);

	real old_mass = cm->mass;

	// For now, keep the disintegrated star and make a log of it.
	// Special case if the merger product ceases to exist.
	//
	// if (((star*)(primary->sbase))->get_element_type()==Disintegrated) {
	//   cerr << "Disintegrated star, set CM mass to zero" << endl;
	//   cm->mass = 0;
	// }
	
	// Must (1) merge the logs of primary and secondary, and      <---
	//      (2) add a line to the log describing the merger.      <---
	
	// Attach the new merger star to cm.
	// Old information of star in CM (if any) is lost here......
	
	cm->sbase = primary->sbase;
	cm->sbase->set_node(cm);
	primary->sbase = NULL;
	
	real new_mass = get_total_mass(cm);
	dm = new_mass - old_mass;
	dv = anomalous_velocity(cm);
    }

    cerr << "new node name = " << cm->format_label() << endl;
    PRL(cpu_time());

    correct_leaf_for_change_of_mass(cm, dm);
    if (cm->mass < 0) {
	cerr << "check and merge, negative mass ";
	PRL(cm->mass);
    }

    cm->set_coll(NULL);			// Set its coll to NULL.
    
    // Add kick velocity, if any...
    
    if (square(dv) > 0)
	correct_leaf_for_change_of_vector(cm, dv, &hdyn::get_vel,
					  &hdyn::inc_vel);
    
    cm->set_oldest_daughter(NULL);

    // Do not try to free memory by simply 
    //		delete bcoll;
    //		delete this;
    // causes a not so very strange error.
    
    real epot, ekin, etot;
#ifdef CALCULATE_POST_COLLISION_ON_GRAPE
    // replaced by GRAPE friendly function (SPZ, March 2001)
    calculate_internal_energies(get_root(), epot, ekin, etot); 
#else
    calculate_energies(get_root(), eps2, epot, ekin, etot);	//dyn function
#endif
    PRC(epot); PRC(ekin); PRL(etot);
    
    real de_total = etot - etot0;
    vector vcm = hdyn_something_relative_to_root(cm, &hdyn::get_vel);
    real de_kick = 0.5 * cm->mass * (vcm*vcm - square(vcm-dv));
    real de_int = -etot_int;
    
    PRC(de_total); PRC(de_int); PRL(de_kick);
    PRL(cpu_time());
    
    // Update statistics on mass and energy change components:
    
    kc->dm_massloss -= dm;	// dm < 0, note

    kc->de_total += de_total;
    kc->de_merge += de_int;
    kc->de_massloss += de_total - de_int - de_kick;
    kc->de_kick += de_kick;
    
    hdyn* root = cm->get_root();
    xreal time = cm->time;
    
    // Special treatment of case cm->mass = 0.
    
    if (cm->mass <= 0) {
	
	// NOTE: Must copy cm's log entry before deleting the node.	<---
	//       Place it in the root log.				<---
	
	remove_node_and_correct_upto_ancestor(cm->get_parent(), cm);
    }

    predict_loworder_all(root, time);

    // Clean up various lists...

    // Check and correct for the presence of either merged node
    // on any perturber list.  Also check and correct neighbor
    // and coll pointers.
    //
    // Probably should perform these checks within merge_nodes.

    hdynptr del[2];
    del[0] = this;
    del[1] = bcoll;

    if (RESOLVE_UNPERTURBED_PERTURBERS || is_pert)
	correct_perturber_lists(get_root(), del, 2, cm);

    correct_nn_pointers(get_root(), del, 2, cm);

    // Also, if cm is on the perturbed binary list, better remove it.

    if (cm->on_perturbed_list()) {
	cerr << "Removing " << cm->format_label()
	     << " from perturbed list " << endl;
	cm->remove_from_perturbed_list();
    } else
	cerr << cm->format_label() << " not on perturbed list " << endl;

    return cm;
}
