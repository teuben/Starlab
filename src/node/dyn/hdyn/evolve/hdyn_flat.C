
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//
//  hdyn_flat.C: functions related to orbit integration within the hdyn class.
//    		     - member functions for FLAT trees only
//.............................................................................
//    version 1:  Jul 1996   Steve McMillan
//    version 2:
//.............................................................................

#include "hdyn.h"

//-----------------------------------------------------------------------------
//  flat_accumulate_acc_and_jerk -- calculates the contribution to the
//                                  acceleration & jerk on this node from
//                                  the node leaf.
//                                  note: for flat trees only
//-----------------------------------------------------------------------------

inline void
hdyn::flat_accumulate_acc_and_jerk(hdyn * leaf,    // attracting leaf node
				   real eps2)	   // softening length squared
{
    vec d_pos = leaf->get_pred_pos() - get_pred_pos();
    vec d_vel = leaf->get_pred_vel() - get_pred_vel();
    real r2inv = 1.0 / (d_pos * d_pos + eps2);
    real a3 = -3.0 * (d_pos * d_vel) * r2inv;
    real mrinv = leaf->get_mass() * sqrt(r2inv);
    pot -= mrinv;
    real mr3inv = mrinv * r2inv;
    acc += mr3inv * d_pos;
    jerk += mr3inv * (d_vel + a3 * d_pos);
}

//-----------------------------------------------------------------------------
//  flat_calculate_acc_and_jerk -- calculates acceleration & jerk on this node:
//                                 from all nodes under p (if p is root ptr);
//                                 or from the node p (if p is a leaf ptr).
//
//                                 NOTE: for flat trees only
//
//-----------------------------------------------------------------------------

#ifdef RECURSIVE_LOOP

void hdyn::flat_calculate_acc_and_jerk(hdyn * p,    // root or leaf
				       real eps2)   // softening length squared
{
    if (p->is_root()) {			// p is root, so
					// invoke for all leaves
	for_all_daughters(hdyn, p, d)	
	    flat_calculate_acc_and_jerk(d, eps2);

    } else {				// p is leaf, but
					// not this leaf
	if (p != this)
	    flat_accumulate_acc_and_jerk(p, eps2);
    }
}

#else

// Looping by recursion (as above) is relatively expensive on many systems.
// Remove the recursion for the sake of efficiency (and readability?) by
// using an explicit loop.

void hdyn::flat_calculate_acc_and_jerk(hdyn * p,    // root or leaf
				       real eps2)   // softening length squared
{
    if (p->is_root()) {			// p is root, so
					// invoke for all leaves
	for_all_daughters(hdyn, p, d)	
	    if (d != this)
		flat_accumulate_acc_and_jerk(d, eps2);

    } else {				// NOTE: else clause currently not used

	if (p != this)
	    flat_accumulate_acc_and_jerk(p, eps2);

    }
}

#endif
