
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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
    kep->set_rel_pos(sister->pos - pos);
    kep->set_rel_vel(sister->vel - vel);
    kep->initialize_from_pos_and_vel();

    sister->kep = kep;
    sister->unperturbed_timestep = unperturbed_timestep;
    sister->fully_unperturbed = fully_unperturbed;
}


void hdyn::reinitialize_kepler_from_hdyn()
{
    update_kepler_from_hdyn();

    // This function is called on restart -- from get_hdyn (hdyn.h), from
    // correct_multiples (kira_init.C), and from initialize_unperturbed(),
    // which is called from full_reinitialize().  (The function is called
    // on system reinitialization for compatibility with restart.)  Since
    // the unperturbed timestep was OK when written out, and should have
    // been picked up from the input snapshot, there should be no need to
    // change anything here... (SLWM 3/98)  

}
