
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// addstar: Add star class to existing node structure.

#include "sstar_to_dyn.h"
#include "single_star.h"
#include "util_io.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  addstar  -- for all particles, add a star part using "new star()".
//-----------------------------------------------------------------------------

void  addstar(dyn * b, real t_current, stellar_type type, bool verbose) {

  addstar(dynamic_cast(node*, b), t_current, type, verbose);

  b->get_starbase()->set_use_hdyn(true);
  
  for_all_daughters(dyn, b, bi) 
    bi->set_radius(bi->get_starbase()
		     ->conv_r_star_to_dyn(
                   bi->get_starbase()->get_effective_radius()));

}

#endif
