
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// readp:  convert ASCII "dumbp" format (mass1, pos1, vel1,
////                                       mass2, pos2, vel2,
////                                       mass3, pos3, vel3,
////                                       etc.)
////
////         data into a Starlab snapshot (flat tree).
////         This is the inverse function to dumbp.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -i    number the particles sequentially [don't number]

//	     Steve McMillan, July 1999

#include "dyn.h"

#ifdef TOOLBOX

int main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:i";

    char *comment;
    bool c_flag = false;
    bool i_flag = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    dyn *root;
    while (root = get_col()) {
	dyn::set_col_output(false);	// force dyn output (default is col)
	put_dyn(root);
	rmtree(root);
    }
}

#endif
