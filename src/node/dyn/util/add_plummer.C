
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// add_plummer.C:  Add Plummer parameters to an input snapshot and write
////                 it out again.  Do not change the input data.
////                 See add_power_law.
////
//// Options:      -c    add comment [none]
////               -C    specify center [(0,0,0)]
////               -M    specify mass [1]
////               -a/R  specify scale [1]
////               -n    force interpretation of parameters in N-body units [no]

//   version 1:  Aug/Sep 2001   Steve McMillan

// This function is just a special case of add_power_law, and is now
// handled as such.

#include "dyn.h"

#ifndef TOOLBOX

void add_plummer(dyn *b,
		 real coeff, real scale,
		 vector center,			// default = (0, 0, 0)
		 bool n_flag,			// default = false
		 bool verbose)			// default = false
{
    add_power_law(b, coeff, 0, scale, center, n_flag, verbose, false);
}

#else

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    real mass = 1, scale = 1;
    vector center = 0;

    bool n_flag = false;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "a:c:C:::M:nR:";

    dyn *b = get_dyn(cin);
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string)) != -1) {
	switch (c) {
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'C':	center = vector(atof(poparr[0]),
					atof(poparr[1]),
					atof(poparr[2]));
			break;
	    case 'M':	mass = atof(poptarg);
			break;
	    case 'a':
	    case 'R':	scale = atof(poptarg);
			break;

	    case 'n':	n_flag = true;
	    		break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    add_plummer(b, mass, scale, center, n_flag, true);
    put_node(cout, *b);
}

#endif
