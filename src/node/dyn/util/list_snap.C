
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// list_snap:  Print times of all snapshots in the input stream.
////
//// Options:    none

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    check_help();

    dyn *b = NULL;
    int count = 0;

    while (b = get_dyn()) {

	cerr << "snap " << ++count;
	real time = getrq(b->get_dyn_story(),"t");

	if (time > -VERY_LARGE_NUMBER)

	    cerr << "  time = " << time << endl;

	else

	    cerr << "  time unknown" << endl;

	rmtree(b);
    }
}

#endif
