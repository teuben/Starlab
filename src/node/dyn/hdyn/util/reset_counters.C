
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// reset_counters:  Clear all time/force counters in the input snapshot(s).
////
//// Options:     -c    add a comment to the output snapshot [false]

// Externally visible function:
//
//	void reset_counters

#include "hdyn.h"

#ifndef TOOLBOX

void reset_counters(hdyn *bi)
{
    bi->inc_steps(-bi->get_steps());		// only have inc and get...
    bi->inc_direct_force(-bi->get_direct_force());
    bi->inc_indirect_force(-bi->get_indirect_force());
}

#else

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

    hdyn *b;

    while (b = get_hdyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	for_all_nodes(hdyn, b, bi)
	    reset_counters(bi);
	put_dyn(cout, *b);	
	delete b;
    }
}

#endif
