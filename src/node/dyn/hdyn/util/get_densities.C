
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// get_densities: compute densities of top-level nodes.
////                Densities are saved in node dyn stories.
////                Only runs on GRAPE systems.

//  Steve McMillan, June 2004

#ifdef TOOLBOX

#include "hdyn.h"

main(int argc, char ** argv)
{
    check_help();

    hdyn *b = get_hdyn();
    b->log_history(argc, argv);

    unsigned int config = kira_config(b);	// default settings
    kira_print_config(config);

    vec cod_pos, cod_vel;
    kira_calculate_densities(b, cod_pos, cod_vel);
    PRL(cod_pos);
    PRL(cod_vel);

    put_node(b);
}

#endif
