
//// pretty_print_tree:  print out the tree structure of the input snapshot(s).
////
//// Options:     none

//   version 1:  Dec 1994   Piet Hut

#include "node.h"

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver to directly print out a tree structure
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
    {
    node *root;    // root node

    check_help();

    while (root = get_node())
	{
	root->pretty_print_tree(cerr);
	rmtree(root);
	}
    }

#endif

/*===========================================================================*/

/* endof: pretty_print_tree.C */

