
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// kira_check.C:  Perform various miscellaneous kira checks.
//
// Externally visible functions:
//
//	bool check_kira_flag
//	bool check_allowed

#include "hdyn.h"

bool check_kira_flag(hdyn* b, char* kira_flag)
{
    // Return true iff the specified kira flag exists and is set to 1.

    if (!find_qmatch(b->get_log_story(), kira_flag)) return false;
    if (getiq(b->get_log_story(), kira_flag) == 1)
	return true;
    else
	return false;
}

bool check_allowed(bool allow_kira_override,
		   char * what_is_allowed,
		   bool verbose, bool& need_skip)
{
    if (allow_kira_override) {

	if (verbose) {
	    cerr << endl
		 << "*** Turning on " << what_is_allowed
		 << endl;
	    need_skip = true;
	}
	return true;

    } else {

	if (verbose) {
	    cerr << endl
		 << "*** Warning: " << what_is_allowed
		 << " is now turned off"
		 << endl;
	    need_skip = true;
	}
    }
    return false;
}

