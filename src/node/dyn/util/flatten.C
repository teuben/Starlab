
//// flatten:  flatten a dyn tree to a single-level linked list under
////           the root node.
////
//// Options:     -c    add a comment to the output snapshot [false]

//   version 1:  March 1994   Piet Hut

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  unbundle_node  --  ...
//                 note:
//                      the positions and velocities of the daughters of the
//                      node-to-be-unbundled are increased by the values of the
//                      corresponding quantities of the node-to-be-unbundled;
//                      the other physical quantities of the daughters remain
//                      unchanged, while those of the node-to-be-unbundled are
//                      be discarded.
//                 note:
//                      the parent, daughter and grand daughters really should
//                      be synchronized in time; to be implemented.
//-----------------------------------------------------------------------------

local void  unbundle_node(dyn * ud)
    {
    dyn * parent;
    parent = ud->get_parent();

    if (parent == NULL)
	err_exit("unbundle() : no parent for this node");

    dyn * od;
    od = ud->get_oldest_daughter();

    if (od == NULL)
	return;                         // nothing left to unbundle

    vector pos = ud->get_pos();
    vector vel = ud->get_vel();

    for_all_daughters(dyn, ud, dd)
	{
	dd->inc_pos(pos);
	dd->inc_vel(vel);
	}

    // pointer adjustments:

    dyn * elder_sister;
    elder_sister = ud->get_elder_sister();

    if (elder_sister == NULL)
	parent->set_oldest_daughter(od);
    else
	{
	elder_sister->set_younger_sister(od);
	od->set_elder_sister(elder_sister);
	}

    dyn * d = od;
    dyn * pd;

    while (d)
	{
	d->set_parent(parent);
	pd = d;
	d = d->get_younger_sister();
	}
    
    dyn * younger_sister;
    younger_sister = ud->get_younger_sister();

    if (younger_sister)
        {
	younger_sister->set_elder_sister(pd);
	pd->set_younger_sister(younger_sister);
	}

    delete ud;
    }

/*-----------------------------------------------------------------------------
 *  flatten_node  --  
 *-----------------------------------------------------------------------------
 */
void  dyn::flatten_node()
    {
    for_all_daughters(dyn, this, d)
	{
	if (d->is_grandparent())
	    d->flatten_node();

	unbundle_node(d);
	}
    }

/*===========================================================================*/

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  flatten_node() as a tool.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    dyn *b;

    while (b = get_dyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        b->flatten_node();
	put_dyn(cout, *b);	
	delete b;
    }
}

#endif

// endof: flatten.C

