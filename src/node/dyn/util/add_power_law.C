
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// add_plummer.C:  Add power-law parameters to an input snapshot and
////                 write it out again.  Do not change the input data.
////                 Use before scaling.
////
//// Options:      -A    specify coefficient [1]
////               -c    add comment [none]
////               -C    specify center [(0,0,0)]
////               -e/E/x/X  specify exponent [0]
////               -a/R  specify scale [1]
////
//// Note that -e 0 gives a Plummer field.

#ifdef TOOLBOX

#include "dyn.h"

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    real coeff = 1, scale = 1, exponent = 0;
    vector center = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "A:a:c:C:::e:E:R:x:X:";

    dyn *b = get_dyn(cin);
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string)) != -1) {
	switch (c) {
	    case 'A':	coeff = atof(poptarg);
			break;
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'C':	center = vector(atof(poparr[0]),
					atof(poparr[1]),
					atof(poparr[2]));
			break;
	    case 'e':
	    case 'E':
	    case 'x':
	    case 'X':
			exponent = atof(poptarg);
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

    putrq(b->get_log_story(), "kira_pl_coeff", coeff);
    putrq(b->get_log_story(), "kira_pl_exponent", exponent);
    putrq(b->get_log_story(), "kira_pl_scale", scale);
    putvq(b->get_log_story(), "kira_pl_center", center);

    cerr << "add_power_law:  A = " << coeff << ", R = " << scale
	 << " x = " << exponent << endl;
    cerr << "                center = " << center << endl;

    put_node(cout, *b);
}
#endif
