
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//  hdyn_grape.C:  Simple switch between GRAPE-4 and GRAPE-6 functions.
//		   Default is to compile GRAPE-4 code, unless GRAPE-6
//		   is specifically selected by setting the variable
//		   STARLAB_HAS_GRAPE6 at compile time.
//
//    version 1:  Jul 2000   Steve McMillan

#if defined(STARLAB_HAS_GRAPE6)

#   include "hdyn_grape6.C"

#else

#   include "hdyn_grape4.C"

#endif
