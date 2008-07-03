
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// sync_time:  Force all times of all nodes to be the same as
////             system_time, and optionally set system_time.
////
////             This is the avenue of last resort in fixing a
////             problematic data set!
////
////             Experimental: May require some fine tuning...
////
//// Steve McMillan (July 2008).

#include "hdyn.h"

#ifdef TOOLBOX

main(int argc, char **argv)
{
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:t:";
    char *comment;
    real syst;
    bool c_flag = false, t_flag = false;
    

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

 	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;

 	    case 't': t_flag = true;
		      syst = atof(poptarg);
		      break;

            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    hdyn *b;

    while (b = get_hdyn()) {
        if (c_flag) b->log_comment(comment);
        b->log_history(argc, argv);

	if (t_flag) b->set_system_time(syst);
	for_all_nodes(hdyn, b, bb) bb->set_time(b->get_system_time());

	put_node(b);
    }
}

#endif
