
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  dyn_kepler.h: declarations of kepler-related functions that know about the
 *		  dyn class
 *.............................................................................
 *    version 1:  Nov 1993   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  ....
 *.............................................................................
 */

#ifndef  STARLAB_DYN_KEPLER_H
#  define  STARLAB_DYN_KEPLER_H

#include  "kepler.h"

void  new_kepler(dyn * cm, real t = 0);
void  dyn_to_kepler(dyn * cm, real t = 0);

void  new_child_kepler(dyn * cm, real t = 0, real circ_limit = 0);
						// Changed by Steve 6/98 to
						// allow circular_binary_limit

void  dyn_to_child_kepler(dyn * cm, real t = 0);

#endif
