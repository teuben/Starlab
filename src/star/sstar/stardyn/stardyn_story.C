
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// dyn_story.C:  initialization of system parameters from input snapshots.
//
// Externally visible functions:
//
//	void check_addstar		// check and set stellar properties

#include "dyn.h"
#include "star/sstar_to_dyn.h"

void check_addstar(dyn* b)
{
    // Check for and add stellar properties.

    if(find_qmatch(b->get_oldest_daughter()
		    ->get_starbase()->get_star_story(), "Type")) {

	real T_start = 0;

	addstar(b,                            // Note that T_start and
		T_start,                      // Main_Sequence are
		Main_Sequence,                // defaults. They are
		true);                        // ignored if a star
	b->set_use_sstar(true);               // is there.
    }
}
