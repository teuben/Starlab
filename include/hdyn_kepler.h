
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  hdyn_kepler.h: declarations of kepler-related functions that know about the
 *		   hdyn class
 *.............................................................................
 *    version 1:  Nov 1993   Piet Hut & Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  ....
 *.............................................................................
 */

#ifndef  STARLAB_HDYN_KEPLER_H
#  define  STARLAB_HDYN_KEPLER_H

#include  "dyn_kepler.h"

void  new_kepler(hdyn *);
void  hdyn_to_kepler(hdyn *);
void  new_child_kepler(hdyn *);
void  hdyn_to_child_kepler(hdyn *);

#endif
