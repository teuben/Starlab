
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Sequentially merge all snapshots in the input stream.  Do *not*
//// rescale masses, positions, or velocities, but always set the
//// system center of mass to 0.  Root node information is retained
//// for all snapshots read.
////
//// Usage: merge_snaps [OPTIONS] < input > output
////
//// Options:
//// None.
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision$", _SRC_);

    dyn *b, *root = NULL;
    int count = 0;

    while (b = get_dyn()) {

	if (root) {

	    // Attach top-level nodes of b to root.

	    dyn *last = root->get_oldest_daughter();
	    while (last->get_younger_sister())
		last = last->get_younger_sister();

	    for_all_daughters(dyn, b, bb) {
		last->set_younger_sister(bb);
		bb->set_elder_sister(last);
		bb->set_parent(root);
		last = bb;
	    }

	    // Merge the root log stories, with a prefix indicating which
	    // initial snapshot this was.  This code is not completely general,
	    // but should be sufficient for simple (typical) stories.

	    count++;

	    story *sr = root->get_log_story();
	    story *s = b->get_log_story();
	    for (story * d = s->get_first_daughter_node(); d != NULL;
		 d = d->get_next_story_node())
		if (!d->get_chapter_flag()) {
		    char tmp[1024];
		    sprintf(tmp, "  +%4.4d:  ", count);
		    strncat(tmp, d->get_text(), 1013);
		    tmp[1023] = '\0';
		    add_story_line(sr, tmp);
		}

	    delete b;

	} else {
	    b->log_history(argc, argv);
	    root = b;
	}
    }

    // Recompute total mass and force the center of mass to 0.

    real mass = 0;
    for_all_daughters(dyn, root, bb) mass += bb->get_mass();
    root->set_mass(mass);

    // Set com using set_com() to update dyn story too.

    root->set_com();

    if (root->get_system_time() == 0.0) {

	// Rewrite initial_mass if it exists.  Delete initial_total_energy
	// and initial_rvirial.

	putrq(root->get_log_story(), "initial_mass", mass);
	rmq(root->get_log_story(), "initial_total_energy");
	rmq(root->get_log_story(), "initial_rvirial");
    }

    put_dyn(root);
    rmtree(root);
}

#endif
