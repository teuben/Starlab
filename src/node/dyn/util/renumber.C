
//// renumber:  renumber all particles sequentially.
////
//// Options:     -c    add a comment to the output snapshot [false]

//   version 1:  Jan 1993   Piet Hut

#include "dyn.h"

#ifdef TOOLBOX

local void renumber(dyn * b)
{
    int  i = 0;

    for_all_leaves(dyn, b, bi)
	bi->set_label(++i);
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

    while (b = get_dyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        renumber(b);
	put_dyn(cout, *b);	
	delete b;
    }
}

#endif

/* endof: renumber.c */
