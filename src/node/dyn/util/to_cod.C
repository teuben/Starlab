
//// to_cod:  Bring all positions and velocities to center-of-density frame.
////          Uses compute_cod and writes to root dyn story.
////          Does not correct the virial radius in the case of a
////          tidal field.
////
//// Options:     -c    add a comment to the output snapshot [false]

//   version 1:  Dec 1992   Simon Portegies Zwart (adusted from to_com) 

#include "dyn.h"

#ifdef TOOLBOX

local void to_cod(dyn *b)
{
    vector cod_pos;
    vector cod_vel;

    compute_max_cod(b, cod_pos, cod_vel);

//    pos += cod_pos;	// Note: The position and velocity of the root node
//    vel += cod_vel;	//       have no particular dynamical significance.
			//	 They are not integrated forward in time,
			// 	 and play no role in any integrator.

    for_all_daughters(dyn, b, bj) {
	bj->inc_pos(-cod_pos);
	bj->inc_vel(-cod_vel);
    }	

    // Correct entries in top-level dyn story.

    cod_pos = vector(0,0,0);
    cod_vel = vector(0,0,0);

    putvq(b->get_dyn_story(), "density_center_pos", cod_pos);
    putvq(b->get_dyn_story(), "density_center_vel", cod_vel);
}

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

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        to_cod(b);
	put_dyn(b);
	delete b;
    }
}

#endif

/* endof: to_cod */

