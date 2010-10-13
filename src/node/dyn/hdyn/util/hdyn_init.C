
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// hdyn_init:  Starlab hdyn-specific initialization functions.

#include "hdyn.h"

// Note that initialize_unperturbed is called by get_hdyn(), so be very
// careful about introducing kira-specific functionality here...

void hdyn::initialize_unperturbed()
{
    // Unperturbed binary motion will have unperturbed_timestep > 0.
    // Fully unperturbed motion will also have fully_unperturbed = true.

    // Unperturbed binary (kepler) details are not saved, so reconstruct
    // them from the hdyn.

//    cerr << "in hdyn_init for " << format_label()
//	 << "at " << time << "/" << system_time << "; ";
//    PRL(init);

    for_all_nodes(hdyn, this, b)
        if (b->is_low_level_node()
	    && !b->get_elder_sister()
	    && b->get_unperturbed_timestep() > 0) {

          //	    cerr << "Reestablishing unperturbed motion for "
          //		 << "(" << b->format_label() << ",";
          //	    cerr << b->get_binary_sister()->format_label() << ")\n";
          //	    PRL(b->get_fully_unperturbed());

	  set_kepler_tolerance(2);

	  // cerr << " to reinitialize_kepler_from_hdyn()" << endl;

	  b->reinitialize_kepler_from_hdyn();	// sister parameters will
						// mirror those of b
      }
}

void hdyn::initialize_slow()
{
    // Slow motion will be already set up for elder binary sister.

    for_all_nodes(hdyn, this, b) {
        if (b->is_low_level_node()
	    && !b->get_elder_sister()
	    && b->get_slow()) {

	    // Elder sister contains all relevant data.  Parameters
	    // for the younger sister simply follow those for b.

	    hdyn *s = b->get_younger_sister();
	    if (s && get_total_energy(b, s) < 0)
		s->slow = slow;
	    else {
		cerr << "*** ignoring slow binary data for unbound pair "
		     << b->get_parent()->format_label() << " ***" << endl;
		b->delete_slow();
	    }
	}
    }
}
