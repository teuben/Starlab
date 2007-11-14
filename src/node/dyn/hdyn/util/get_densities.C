
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
////		  -i    output individual bound specifications [false]
////		  -v    verbose [false]
////              -c    transpose the cluster to CM frame [yes]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.


//  Steve McMillan, June 2004

#ifdef TOOLBOX

#include "hdyn.h"

static real cod_pot = 0;
local void set_cod_pot(hdyn *b, vec cod_pos)
{
    cod_pot = 0;
    if (b->get_external_field()) {
	vec tmp;
	get_external_acc(b, cod_pos, tmp, cod_pot, tmp, tmp, true);
    }
}

main(int argc, char ** argv)
{
    check_help();

    bool verbose = false;
    bool individual = true;
    bool cm_flag = true;

    extern char *poptarg;
    int c;
    char* param_string = "ci0";

    bool force_nogrape = false;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case '0':	force_nogrape = true;
			break;
	    case 'i': individual = !individual;
	              break;
	    case 'c': cm_flag = !cm_flag;
	              break;
	    case 'v': verbose = !verbose;
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
    //cod_pos = b->get_pos();
    //cod_vel = b->get_vel();

    if (individual) {

      if (cm_flag) {
	for_all_nodes(hdyn, b, bi) {
	  bi->inc_pos(-cod_pos);
	  bi->inc_vel(-cod_vel);
	}
      }

      set_cod_pot(b, cod_pos);

      int bound;
      for_all_daughters(hdyn, b, bi) {
	bound = 0;
	if (0.5*square(bi->get_vel()-cod_vel) + bi->get_pot() - cod_pot < 0) 
	  bound = 1;
	      
	cerr << bound << " " 
	     << bi->get_index() << " "
	     << bi->get_mass() << " "
	     << bi->get_pos() << " "
	     << bi->get_vel() << endl;
      }
    }

    if(verbose)
      put_node(b);
}

#endif
