
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// dumbp:  dump out an N-body system in a dumb format suitable
////         for digestion by NBODY1-5 and starcluster.
////         This is the inverse function to readp.
////
//// Options:     -p    specify precision of output [6 sig. fig.]
////              -t    include time in output [no]

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "p:t";

    int p = STD_PRECISION;
    bool time = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'p': p = atoi(poptarg);
		      break;
	    case 't': time = !time;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

    cout.precision(p);

    dyn *root;
    while (root = get_dyn()) {
	put_col(root);
	rmtree(root);
    }
}

#endif

// end of: dumbp.C
