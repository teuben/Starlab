
//// freeze:  set all velocities to zero, leaving the positions unchanged.
////
//// Options:     -c    add a comment to the output snapshot [false]

//   version 1:  Dec 1992   Piet Hut

#include "dyn.h"

#ifdef TOOLBOX

local void freeze(dyn * b)
    {
    dyn * bi;

    for (bi=b->get_oldest_daughter(); bi != NULL; bi=bi->get_younger_sister())
	bi->scale_vel(0);
    }

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

    while (b = get_dyn(cin))
	{
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        freeze(b);
	put_dyn(cout, *b);	
	rmtree(b);
	}
}

#endif

// endof: freeze.C
