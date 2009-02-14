
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// user_diag: User-provided diagnostic code.  This function is called
// at the end of sys_stats() to provide additional problem-specific
// output from the standalone sys_stats or kira tools.


#include "dyn.h"
#ifndef TOOLBOX

void user_diag(dyn *b)
{
  cerr << "  user_diag: No user-provided diagnostic code" << endl;
}

#endif
