
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// pretty_print_tree:  print out the tree structure of the input snapshot(s).
////
//// Options:     none

//   version 1:  Dec 1994   Piet Hut

#include "node.h"

//===========================================================================

#ifdef TOOLBOX

//---------------------------------------------------------------------------
//  main  --  driver to directly print out a tree structure
//---------------------------------------------------------------------------

main(int argc, char ** argv)
{
    node *root;    // root node

    check_help();

    while (root = get_node()) {
	root->pretty_print_tree();
	rmtree(root);
    }
}

#endif
