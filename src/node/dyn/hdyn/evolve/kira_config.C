
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// Handle basic configuration of the available hardware/software.
// See also the other soft switches in kira_ev.C, kira_energy.C,
// and kira_density.C.
//
// Externally visible functions:
//
//	unsigned int kira_config
//	void kira_print_config
//	void check_release_grape
//	void clean_up_hdyn_grape
//
// A non-zero return from kira_config means we have some "extra help"
// available to us, e.g. in the form of a GRAPE or (in future) a treecode.

#include "hdyn.h"
#include "grape4.h"
#include "grape6.h"

// (Add further relevant include files here.)

// Return values:
//
//	no bits set (i.e. 0)			no GRAPE, no treecode, etc.
//	bit 1 set				GRAPE-4 present
//	bit 2 set				GRAPE-6 present
//
//	(etc.)
//
// Note that for now we generally assume that the system doesn't have
// both a GRAPE-4 and a GRAPE-6 attached...


// kira_config: determine whether we have GRAPE hardware available, and
//		set up the function pointer to compute top-level accs
//		and jerks.

unsigned int kira_config(hdyn *b,
			 int force_config)	// default = -1
{
    unsigned int true_config = 0;

    //----------------------------------------------------------------
    // We currently use a positive number of pipes as an indicator
    // that we have a GRAPE system.  This actually only tells us that
    // we have a GRAPE library and an accessible configuration file.
    //
    // Could also check for the existence of the PCIbus interface
    // card, but we currently have no definitive test that a working
    // GRAPE is actually attached to it.  To check for the card on a
    // linux system, use "lspci -v" and look for a "9060" (GRAPE-4)
    // or "9080" (GRAPE-6) device.  Probably better to do this as a
    // configuration option.

    if (h3npipe_() > 0) true_config += 1;	// GRAPE-4 check
    if (g6_npipes_() > 0) true_config += 2;	// GRAPE-6 check

    // PRL(true_config);

    // (Add further determining factors here.)

    //----------------------------------------------------------------

    unsigned int config = true_config;
    if (force_config >= 0) config = force_config;

    b->set_config(config);
//    PRC(force_config); PRC(config); PRL(b->get_config());

#if 0

    // Set up acc_and_jerk function pointer.
    // Precedence (if not forced) is: GRAPE-6, GRAPE-4, front-end.

    acc_function_ptr fn = NULL;

    // Current options for kira_calculate_top_level_acc_and_jerk
    // are as follows:

    if (config == 0)
	fn = top_level_acc_and_jerk_for_list;
    else if (config%2 == 0) 
	fn = grape6_calculate_acc_and_jerk;
    else
	fn = grape4_calculate_acc_and_jerk;

    b->set_kira_calculate_top_level_acc_and_jerk(fn);

#endif

    return true_config;
}

void kira_print_config(unsigned int config)
{
    if (config <= 0)
        cerr << "no hardware acceleration" << endl;
    else if (config%2 == 0)
        cerr << "found GRAPE-6 acceleration" << endl;
    else
        cerr << "found GRAPE-4 acceleration" << endl;
}

void check_release_grape(unsigned int config,
			 kira_options *ko,
			 xreal time,
			 bool verbose)		// default = true
{
    if (config > 0) {
	if (config%2 == 0)
	    check_release_grape6(ko, time, verbose);
	else
	    check_release_grape4(ko, time, verbose);
    }
}

void clean_up_hdyn_grape(unsigned int config)
{
    if (config > 0) {
	if (config%2 == 0)
	    clean_up_hdyn_grape6();
	else
	    clean_up_hdyn_grape4();
    }
}
