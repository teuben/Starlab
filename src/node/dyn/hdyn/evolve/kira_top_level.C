
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

#include "hdyn.h"

// The n_top_level functions are a convenient way of keeping track
// of the total number of top-level nodes.  *Much* more efficient than
// using root->n_daughters().  They should probably be made part of
// the node class eventually...
//
// n_top_level is maintained via set_n_top_level in kira.C as the tree
// changes.  Currently, get_n_top_level is used only in kira_ev.C.

static int n_top_level = 0;

void set_n_top_level(hdyn* b)
{
    n_top_level = b->n_daughters();
}

int get_n_top_level()
{
    return n_top_level;
}
