
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// reflect_velocities:  multiply all velocity components by -1,
////                      thereby effectively reversing the direction
////                      of time in the input N-body system(s).
////
//// Options:      -c     add a comment to the output snapshot [false]

//   Piet Hut, Nov. 1994

#include "dyn.h"

#ifdef TOOLBOX

local void  flip_velocities(dyn * b)
{
    dyn * bi;

    b->scale_vel(-1);

    for (bi=b->get_oldest_daughter(); bi != NULL; bi=bi->get_younger_sister())
	bi->scale_vel(-1);
}


main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
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

        flip_velocities(b);
	put_dyn(b);
	rmtree(b);
    }
}

#endif

// end of: reflect_velocities.C
