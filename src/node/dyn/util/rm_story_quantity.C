
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Remove a line of the form
////
////	keyword = value
////
//// to the root dyn story an input snapshot.  Keyword is treated as a
//// string.
////
//// Usage: rm_story_quantity keyword < input > output
////
//// Written by Steve McMillan, 8/09.  This is a dyn function to allow
//// it to deal with "col" format data.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

int main(int argc, char ** argv)
{
    bool rm_story = (argc > 1);
    dyn *b;

    while (b = get_dyn()) {

	if (rm_story) rmq(b->get_dyn_story(), argv[1]);

	put_dyn(b);
	rmtree(b);
    }
    return 0;
}
