
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
////                 Use before scaling.
////
//// Options:      -c    add comment [none]
////               -C    specify center [(0,0,0)]
////               -M    specify mass [1]
////               -a/R  specify scale [1]

#ifdef TOOLBOX

#include "dyn.h"

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    real mass = 1, scale = 1;
    vector center = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "a:c:C:::M:R:";

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
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    if (c_flag)
	b->log_comment(comment);

    putrq(b->get_log_story(), "kira_plummer_mass", mass);
    putrq(b->get_log_story(), "kira_plummer_scale", scale);
    putvq(b->get_log_story(), "kira_plummer_center", center);

    cerr << "add_plummer:  M = " << mass << ", R = " << scale << endl;
    cerr << "              center = " << center << endl;

    put_node(cout, *b);
}
#endif
