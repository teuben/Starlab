
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// hdyn_kepler.C:  kepler-specific hdyn functions.

// Externally visible functions:
//
//	void hdyn::update_kepler_from_hdyn
//	void hdyn::reinitialize_kepler_from_hdyn

#include "hdyn.h"

void new_kepler(hdyn * com)	// com is the center-of-mass hdyn
{
    new_kepler((dyn *) com, com->get_time());
}


void hdyn_to_kepler(hdyn * com)
{
    dyn_to_kepler(com, com->get_time());
}


void hdyn::update_kepler_from_hdyn()
{
    set_kepler_tolerance(2);	// avoid termination on kepler error

    if (kep == NULL) kep = new kepler;
    hdyn *sister = get_binary_sister();

    kep->set_time(time);
    kep->set_total_mass(parent->get_mass());
    kep->set_rel_pos(pos - sister->pos);
    kep->set_rel_vel(vel - sister->vel);
    kep->initialize_from_pos_and_vel();

    sister->kep = kep;
    sister->unperturbed_timestep = unperturbed_timestep;
    sister->fully_unperturbed = fully_unperturbed;
}


void hdyn::reinitialize_kepler_from_hdyn()
{
    update_kepler_from_hdyn();

    // This function is only called on restart -- from get_hdyn (hdyn.h)
    // and from correct_multiples (kira_init.C).  Since the unperturbed
    // timestep was OK when written out, and should have been picked up
    // from the input snapshot, there should be no need to change anything
    // here... (SLWM 3/98).  This function is now also called on system
    // reinitialization, for compatibility with restart.

#if 0

    // Never advance unperturbed motion beyond parent's next step,
    // but be sure that the step doesn't end at a bad place..

    //int usteps = get_unperturbed_steps(true);
    real usteps = get_unperturbed_steps(true);
    //unsigned long usteps = get_unperturbed_steps(true);

    if (usteps > 0) {
	unperturbed_timestep = timestep*usteps;
	get_binary_sister()->unperturbed_timestep = unperturbed_timestep;
    }

#endif

}
