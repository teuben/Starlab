
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// (Re)scale an N-body system to specified M, Q (=T/U), R, or E.
//// This hdyn version is identical to the dyn "scale", except that
//// GRAPE-4 or GRAPE-6 is used if available.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//  Steve McMillan, July 2000; revised July 2001, June 2004

#ifdef TOOLBOX

#include "hdyn.h"

// Code is identical to the  main program in dyn/util/scale.C, except
// that get_top_level_energies() is replaced by kira_top_level_energies(),
// which uses GRAPE if possible.

main(int argc, char ** argv)
{
    check_help();

    // Look for "local" command-line arguments first.  Note that the
    // dyn version is explicitly coded to accept and ignore this option.
    // Using pgetopt seems to break the code in parse_scale_main(), so
    // loop directly over the argument list here.

    bool force_nogrape = false;

    for (int i = 1; i < argc; i++)
	if (streq(argv[i], "-0")) {
	    force_nogrape = true;
	    break;
	}

    // Parse the remaining options using the dyn version parser.

    real m = 0, q = -1, e = 0, r = 0;
    real eps = 0;

    bool debug = false;
    bool c_flag = false;
    bool e_flag = false;
    bool m_flag = false;
    bool q_flag = false;
    bool r_flag = false;

    if (!parse_scale_main(argc, argv, eps,
			  c_flag,
			  e_flag, e, m_flag, m,
			  q_flag, q, r_flag, r,
			  debug,
			  "$Revision$", _SRC_)) {
	get_help();
	exit(1);
    }

    hdyn *b = get_hdyn();
    b->log_history(argc, argv);

    unsigned int config = kira_config(b);	// default settings
    kira_print_config(config);

    if (config && force_nogrape) {
        kira_config(b, 0);
        cerr << "GRAPE suppressed" << endl;
    }

    scale(b, eps, c_flag, e_flag, e, m_flag, m, q_flag, q, r_flag, r,
	  debug, kira_top_level_energies);

    put_node(b);
}

#endif
