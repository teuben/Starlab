
//// to_com:  Bring all positions and velocities to center-of-mass frame.
////          Uses compute_com and writes to root dyn story.
////          Does not correct the virial radius in the case of a
////          tidal field.
////
//// Options:     -c    add a comment to the output snapshot [false]

//   version 1:  Dec 1992   Piet Hut

#include "dyn.h"

#ifndef TOOLBOX

void dyn::to_com()
{
    vector com_pos;
    vector com_vel;

    compute_com(this, com_pos, com_vel);

//    pos += com_pos;	// Note: The position and velocity of the root node
//    vel += com_vel;	//       have no particular dynamical significance.
			//	 They are not integrated forward in time,
			// 	 and play no role in any integrator.

    for_all_daughters(dyn, this, bj) {
	bj->pos -= com_pos;
	bj->vel -= com_vel;
    }

    // Correct entries in top-level dyn story.

    com_pos = vector(0,0,0);
    com_vel = vector(0,0,0);

    putvq(get_dyn_story(), "com_pos", com_pos);
    putvq(get_dyn_story(), "com_vel", com_vel);
}

#else

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    dyn *b;

    while (b = get_dyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        b->to_com();
	put_dyn(cout, *b);	
	delete b;
    }
}

#endif

// endof: to_com

