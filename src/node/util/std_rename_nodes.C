
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// std_rename_nodes:  Rename all CM nodes bottom up to reflect the standard
////                    Starlab naming scheme (a,b).
////
//// Options:      none

//		   Steve McMillan, Sep 2003

#ifdef TOOLBOX

#include "node.h"

// Should probably be moved to node and made global (Steve, 9/03).

local void name_parent_from_components(node *od, node *yd)
{
    char name1[256], name[256];
    strcpy(name1, od->format_label());
    sprintf(name, "(%s,%s)", name1, yd->format_label());
    od->get_parent()->set_name(name);
    od->get_parent()->set_index(-1);
}

int main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    node *b;
    while (b = get_node()) {
	b->log_history(argc, argv);

	for_all_leaves(node, b, bb)
	    if (bb->is_low_level_node()) {
		node *sis = bb->get_elder_sister();
		if (!sis || sis->is_parent()) {
		    if (!sis)
			name_parent_from_components(bb,
						    bb->get_younger_sister());
		    else
			name_parent_from_components(sis, bb);
		}
	    }

	put_node(b);
	rmtree(b);
    }
}
#endif
