
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add a line of the form
////
////	keyword = value
////
//// to the root dyn story an input snapshot.  Keyword and value are
//// treated as strings.
////
//// Usage: add_story_quantity keyword value < input > output
////
//// Written by Steve McMillan, 8/09.  This is a dyn function to allow
//// it to deal with "col" format data.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

int main(int argc, char ** argv)
{
    bool add_story = (argc > 2);
    dyn *b;

    while (b = get_dyn()) {

	if (add_story) putsq(b->get_dyn_story(), argv[1], argv[2]);

	put_dyn(b);
	rmtree(b);
    }
    return 0;
}
