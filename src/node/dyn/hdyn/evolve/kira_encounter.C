
#include "hdyn.h"
#include <star/dstar_to_kira.h>
#include <star/single_star.h>

real get_sum_of_radii(hdyn* bi, hdyn* bj) {

    // by default
    real sum_of_radii = bi->get_radius()+bj->get_radius();

    // For black hole use tidal horizon instead of radius
    if(bi->get_starbase()->get_element_type() == Black_Hole ||
       (find_qmatch(bi->get_log_story(), "black_hole") &&
        getiq(bi->get_log_story(), "black_hole")==1)) {

        if(!(bi->get_starbase()->get_element_type() == Black_Hole ||
             find_qmatch(bj->get_log_story(), "black_hole"))) {

            sum_of_radii = 2 * bj->get_radius()
                * (pow(bi->get_mass()/bj->get_mass(), ONE_THIRD));
	  }
      }
    else if(bj->get_starbase()->get_element_type() == Black_Hole ||
            (find_qmatch(bj->get_log_story(), "black_hole") &&
             getiq(bj->get_log_story(), "black_hole")==1)) {

        sum_of_radii = 2 * bi->get_radius()
            * (pow(bj->get_mass()/bi->get_mass(), ONE_THIRD));
      }

    return sum_of_radii;
  }

real print_encounter_elements(hdyn* bi, hdyn* bj,
			      char* s, 		    // default = "Collision"
			      bool verbose = true)
{
    kepler k;
    initialize_kepler_from_dyn_pair(k, bi, bj, true);	// minimal kepler
							// -- no phase info

    cerr << endl;
    cerr << s << " at time = " << bi->get_time();

    if (bi->get_use_sstar() && verbose) {

	cerr << " ("
	     << bi->get_starbase()->conv_t_dyn_to_star(bi->get_time())
	     << " [Myr])";
	cerr << " between \n"
	     << "     " << bi->format_label() << " (";
	put_state(make_star_state(bi), cerr);
	cerr << "; M = " << bi->get_starbase()->get_total_mass()
	     << " [Msun], "
	     << " R = " << bi->get_starbase()->get_effective_radius()
	     << " [Rsun]) and\n"
	     << "     " << bj->format_label()
	     << " (";
	put_state(make_star_state(bj), cerr);
	cerr << "; M = " << bj->get_starbase()->get_total_mass()
	     << " [Msun], "
	     << " R = " << bj->get_starbase()->get_effective_radius()
	     << " [Rsun]) "
	     <<"\n     at distance "
	     << k.get_periastron()
	     << " ("
	     << bi->get_starbase()->conv_r_dyn_to_star(k.get_periastron())
	     << " [Rsun]).";

    } else {

      cerr << " between " << bi->format_label();
      if(bi->get_starbase()->get_element_type() == Black_Hole ||
	 (find_qmatch(bi->get_log_story(), "black_hole") &&
	  getiq(bi->get_log_story(), "black_hole")==1))
	cerr << " [bh] ";
      
      cerr << " and " << bj->format_label();
      if(bj->get_starbase()->get_element_type() == Black_Hole ||
	 (find_qmatch(bj->get_log_story(), "black_hole") &&
	  getiq(bj->get_log_story(), "black_hole")==1))
	cerr << " [bh] ";
      
    }

    cerr << endl;

    if (k.get_energy() < 0) {
        cerr << "     Orbital parameters: ";

        print_binary_params(&k, bi->get_mass(), 0.0,
			    abs(bi->get_pos()), verbose, 10, 10);
	cerr << endl;
    }
    else
	cerr << "     E = " << k.get_energy() << endl;
      

    return k.get_energy();
}


void check_print_close_encounter(hdyn *bi)
{
    if (bi && bi->is_valid()) {

        hdyn *nn = bi->get_nn();

	// Use nearest neighbor for encounter diagnostics.
	// Only consider encounters between nodes at the top level
	// or the next level down (added by Steve 4/99 to prevent
	// multiple output in bound hierarchical systems).

	if (nn && bi->is_leaf()
	    && (bi->is_top_level_node()
		|| bi->get_parent()->is_top_level_node())
	    && bi->get_kepler() == NULL
	    && bi->get_nn() != NULL
	    && nn->is_valid()
	    && nn->is_leaf()
	    && (nn->is_top_level_node()
		|| nn->get_parent()->is_top_level_node())
	    && nn->get_kepler() == NULL
	    && bi->get_d_nn_sq()
	    		<= bi->get_stellar_encounter_criterion_sq())

		print_close_encounter(bi, nn);
    }
}

// Flag close encounter distances between stars.  Output occurs on any
// unbound encounter, and on the first encounter in a bound system.
// Self-contained function, using the log story to propogate data.

void print_close_encounter(hdyn* bi,
			   hdyn* bj)
{
    if (bi->get_mass() < bj->get_mass())  // Note: this "<" means that equal-
	return;				  // mass stars get printed twice!

    real d2cc_2 = -1;
    real d2cc_1 = -1;
    real d2cc   = -1;

    // Variables:	bi is primary, bj is current coll
    //			cc_name is label of previous coll
    //			cc_time is time of previous coll
    //			d2cc is squared distance to coll
    //			d2cc_1 is previous d2cc
    //			d2cc_2 is previous d2cc_1
    //			pcc_name is label of coll at last pericenter
    //			pcc_time is time of last pericenter
    //			pcc_cntr counts bound pericenters
    //
    // All but bi and bj are stored in the bi log story.

    char *cc_name = getsq(bi->get_log_story(), "cc_name");

#if 0
    cerr << "\nIn print_close_encounter with bi =  " << bi->format_label()
	 << "  at time " << bi->get_system_time() << endl;
    cerr << "     "; print_coll(bi,2);
    cerr << "     "; print_nn(bi,2); 
    cerr << "     (bj = " << bj->format_label();
    cerr << ":  coll = "; print_coll(bj);
    cerr << ",  nn = "; print_nn(bj); cerr << ")" << endl;

    cerr << "     cc_name = " << getsq(bi->get_log_story(), "cc_name")
	 << endl;
#endif
  
    if (cc_name && streq(cc_name, bj->format_label())) {

	// Retrieve coll information from the log story.

	d2cc_2 = getrq(bi->get_log_story(), "d2cc_1");
	d2cc_1 = getrq(bi->get_log_story(), "d2cc");
	d2cc   = square(bi->get_pos() - bj->get_pos());

	if (d2cc_2 > d2cc_1 && d2cc >= d2cc_1) {    // just passed pericenter

	    int pcc_cntr = 0;
	    real E = get_total_energy(bi, bj);

	    if (E > 0) {

	        // Always print unbound encounter elements.

		print_encounter_elements(bi, bj, "Close encounter");

	    }
	    else {

	        // Only print first bound encounter.

		char *pcc_name = getsq(bi->get_log_story(), "pcc_name");
	
		if (pcc_name && streq(pcc_name, bj->format_label())) {
		    pcc_cntr = getiq(bi->get_log_story(), "pcc_cntr");
		    if (pcc_cntr == -1) 
			pcc_cntr=0;
		}

		pcc_cntr++;

		if (pcc_cntr == 1) 
		    print_encounter_elements(bi, bj, "First bound encounter");
	    }

	    // Save data on last pericenter (bound or unbound).
            // Note that an unbound pericenter resets pcc_cntr to zero.

	    putiq(bi->get_log_story(), "pcc_cntr", pcc_cntr);
	    putsq(bi->get_log_story(), "pcc_name", bj->format_label());
	    putrq(bi->get_log_story(), "pcc_time", bi->get_time());  

	}

    }
    else if (getsq(bi->get_log_story(), "pcc_name"))

	if (cc_name == NULL
	    || strcmp(cc_name, getsq(bi->get_log_story(), "pcc_name"))) {

	    // Unlikely multiple encounter(?).  Has been known to occur
	    // when one incoming star overtakes another.

	    cerr << endl << "print_close_encounter:  "
	         << "bi = " << bi->format_label()
		 << " at time " << bi->get_system_time() << endl;
	    cerr << "     current coll = " << bj->format_label();
	    cerr << "  previous coll = ";
	    if (cc_name)
	        cerr << cc_name;
	    else
	        cerr << "(null)";
	    cerr << endl;
	    cerr << "     coll at last pericenter = "
	         << getsq(bi->get_log_story(), "pcc_name");
	    cerr << "  (periastron counter = "
	         << getiq(bi->get_log_story(), "pcc_cntr") << ")\n";
	}

    putsq(bi->get_log_story(), "cc_name", bj->format_label());
    putrq(bi->get_log_story(), "d2cc_1", d2cc_1);  
    putrq(bi->get_log_story(), "d2cc", d2cc);  
    putrq(bi->get_log_story(), "cc_time", bi->get_time());  

#if 0
    cerr << "\nAt end of print_close_encounter at time "
	 << bi->get_system_time() << endl;
    cerr << "bi = " << bi->format_label();
    cerr << ", cc_name = " << getsq(bi->get_log_story(), "cc_name") << endl;
#endif

}

