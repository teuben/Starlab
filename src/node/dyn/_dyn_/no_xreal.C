
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// no_xreal:  Rewrite the input with real (not xreal) data.
////
//// Options: none

#include "_dyn_.h"

#ifdef TOOLBOX

main(int argc, char** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision$", _SRC_);

    _dyn_  * b;

    while (b = get__dyn_()) {
        put_node(b, cout, false);
	delete b;
    }
}

#endif
