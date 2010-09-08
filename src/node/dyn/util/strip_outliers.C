
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Strip outlying centers of mass from an N-body system.  If a
//// current density center is found in the root dyn story, use it as
//// the center.  Otherwise, if a valid center of mass is found, it is
//// used.  If neither center is found, the geometric center is used.
////
//// Usage: strip_outliers [OPTIONS] < input > output
////
//// Options:     
////	          -c    add a comment to the output snapshot [false]
////              -r    specify stripping radius [100]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//-----------------------------------------------------------------------------
//   version 1:  Sep 2010   Steve McMillan	steve@physics.drexel.edu
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifdef TOOLBOX

int strip_outliers(dyn *b, real r2, bool verbose = false)
{
    vec cpos, cvel;
    get_std_center(b, cpos, cvel);
    r2 = r2*r2;

    int nstrip = 0;
    dyn *bb = b->get_oldest_daughter();
    while (bb) {
	dyn *next = bb->get_younger_sister();
	if (square(bb->get_pos()-cpos) > r2) {
	    if (verbose)
		cerr << "removing node " << bb->format_label()
		     << " at time " << b->get_system_time() << endl << flush;
	    detach_node_from_general_tree(bb);
	    rmtree(bb);
	    nstrip++;
	}
	bb = next;
    }

    // Correct the root node mass.  Do NOT modify its position or
    // velocity, or change any center of mass information.

    real m = 0;
    for_all_daughters(dyn, b, bbb) m += bbb->get_mass();
    b->set_mass(m);

    return nstrip;
}

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = false;      // if TRUE: a comment given on command line
    real r = 100;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:r:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'r': r = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);

        b->log_history(argc, argv);
	int nstrip = strip_outliers(b, r);
	cerr << "stripped " << nstrip << " nodes at time "
	     << b->get_system_time() << endl;

	put_dyn(b);
	rmtree(b);
    }
}

#endif
