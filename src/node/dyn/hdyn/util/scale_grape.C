
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// scale_grape:
////
////         (re)scale an N-body system to specified M, Q (=T/U), R, or E.
////         This version is identical to the dyn "scale", except that
////         GRAPE-4 or GRAPE-6 is used if available.

//  Steve McMillan, July 2000; revised July 2001

#ifdef TOOLBOX

#include "hdyn.h"

#define USE_GRAPE
#include "../evolve/kira_grape_include.C"	// GRAPE-specific functions:
						// invoke by setting USE_GRAPE

// Code is identical to the  main program in dyn/util/scale.C, except
// that get_top_level_energies() is replaced by kira_top_level_energies(),
// which uses GRAPE if possible.

main(int argc, char ** argv)
{
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
			  debug)) {
	get_help();
	exit(1);
    }

    hdyn *b = get_hdyn(cin);
    b->log_history(argc, argv);

    scale(b, eps, c_flag, e_flag, e, m_flag, m, q_flag, q, r_flag, r,
	  debug, kira_top_level_energies);

    put_node(cout, *b);
}

#endif
