
//// merge:  merge all low-level nodes
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -f    merge the specified fraction, chosen randomly [1]

//   version 1:  September 1999		Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void merge_low_level_nodes(dyn * b,	    // root node
			   real frac)	    // default = 1
{
    for_all_daughters(dyn, b, bb)
	if (bb->is_parent()
	    && randinter(0,1) <= frac) {    // N.B. really random, not "first
					    // frac" or "every 1/frac", ...

	    dyn * od = bb->get_oldest_daughter();
	    dyn * yd = od->get_younger_sister();

	    // Assume that CM node is consistent, so all we need to
	    // do is delete the daughters.

	    rmtree(od);
	    rmtree(yd);

	    bb->set_oldest_daughter(NULL);
	}
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  flatten_node() as a tool.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    real frac = 1;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:f:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'f': frac = atof(poptarg);
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

        merge_low_level_nodes(b, frac);
	put_dyn(cout, *b);	
	delete b;
    }
}

#endif
