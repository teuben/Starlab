
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

void hdyn::initialize_unperturbed()
{
    // Unperturbed binary motion will have unperturbed_timestep > 0.
    // Fully inperturbed motion will also have fully_unperturbed = true.

    // Unperturbed binary (kepler) details are not saved, so reconstruct
    // them from the hdyn.

    for_all_nodes(hdyn, this, b)
        if (b->is_low_level_node()
	    && !b->get_elder_sister()
	    && b->get_unperturbed_timestep() > 0) {

          //	    cerr << "Reestablishing unperturbed motion for "
          //		 << "(" << b->format_label() << ",";
          //	    cerr << b->get_binary_sister()->format_label() << ")\n";
          //	    PRL(b->get_fully_unperturbed());

	  set_kepler_tolerance(2);
	  b->reinitialize_kepler_from_hdyn();	// sister parameters will
						// mirror those of b
      }
}

void hdyn::initialize_slow()
{
    // Slow motion will be already set up for elder binary sister.

    for_all_nodes(hdyn, this, b)
        if (b->is_low_level_node()
	    && !b->get_elder_sister()
	    && b->get_slow()) {

	    // Elder sister contains all relevant data.  Parameters
	    // for the younger sister simply follow those for b.

	    get_younger_sister()->slow = slow;
      }
}
