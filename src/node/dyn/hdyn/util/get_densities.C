
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// Compute densities of top-level nodes.
//// Densities are saved in node dyn stories. The program only actually
//// does the computation if a GRAPE is attached.
////
//// Usage: get_densities
////
//// Options:
//// None.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.


//  Steve McMillan, June 2004

#ifdef TOOLBOX

#include "hdyn.h"

main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "0";

    bool force_nogrape = false;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case '0':	force_nogrape = true;
			break;

	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }

    hdyn *b = get_hdyn();
    b->log_history(argc, argv);

    unsigned int config = kira_config(b);	// default settings
    kira_print_config(config);

    if (config && force_nogrape) {
        kira_config(b, 0);
        cerr << "GRAPE suppressed" << endl;
    }

    vec cod_pos, cod_vel;
    kira_calculate_densities(b, cod_pos, cod_vel);
    PRL(cod_pos);
    PRL(cod_vel);

    put_node(b);
}

#endif
