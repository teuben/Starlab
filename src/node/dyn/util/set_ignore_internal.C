
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// set_ignore_internal.C:  Add a flag specifying that internal
////                         interactions should be suppressed
////                         i.e. create a system of test particles.
////
//// Options:      none

//   version 1:  Dec 2001   Steve McMillan

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char *argv[])
{
    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "";

    dyn *b = get_dyn(cin);
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);

    // Parse the argument list:

//      while ((c = pgetopt(argc, argv, param_string)) != -1) {
//  	switch (c) {
//  	    default:
//  	    case '?':	params_to_usage(cerr, argv[0], param_string);
//  			get_help();
//  			return false;
//  	}
//      }

    putrq(b->get_log_story(), "ignore_internal", 1);
    put_node(cout, *b);
}

#endif
