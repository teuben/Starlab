
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// kira_diag:  Starlab kira/diag-specific functions.

#include "hdyn.h"

#ifndef TOOLBOX

kira_diag::kira_diag() {		    // (Steve's 6/00
	    				    //  preferences...)

    name					= NULL;
    t1						= -VERY_LARGE_NUMBER;
    t2						= VERY_LARGE_NUMBER;

    kira_runtime_flag				= 0;

    kira_main					= false;
    check_heartbeat				= true;
    n_check_heartbeat				= 50000;
    n_check_runtime				= 2500;

    unpert_function_id				= false;
    report_start_unperturbed			= false;	// SPZ mod
    report_continue_unperturbed			= false;
    report_end_unperturbed			= false;	// SPZ mod
    report_pericenter_reflection		= false;
    report_impending_multiple_status		= false;
    report_zero_unpert_steps			= false;
    report_multiple				= true;
    unpert_report_level				= 2;
    end_unpert_report_level			= 2;
    multiple_report_level			= 2;

    tree					= true;
    tree_level					= 0;

    ev_function_id				= false;
    ev						= false;
    grape					= false;	// SLWM 5/02
    grape_level					= 0;
    timestep_check				= false;	// SLWM 5/02
    correct					= false;
    slow_perturbed				= false;

    kira_ev					= false;

    slow					= true;
    slow_level					= 0;

    report_stellar_evolution			= false;
    report_stellar_mass_loss			= false;
    report_binary_mass_loss			= false;

    report_adjust_perturbed_list		= false;

    report_correct_perturber_list		= false;

    report_close_encounters			= 1;

    report_kepler_trig_error			= false;
}

bool kira_diag::check_diag(hdyn *b)
{
    // Acceptance criteria involve both system time and node label.
    // Only a single time range and label are allowed, for now.

    if (b->get_system_time() < t1 || b->get_system_time() > t2) return false;

    if (name && !clump_contains(b, name)) return false;
	
    return true;
}

#define PRS(x) s << "    " << #x << " = " << x << endl

void kira_diag::print(ostream &s)
{
    s << endl << "kira_diag settings:" << endl;

    PRS(kira_main);

    PRS(kira_runtime_flag);
    PRS(check_heartbeat);
    PRS(n_check_heartbeat);
    PRS(n_check_runtime);

    PRS(unpert_function_id);
    PRS(report_start_unperturbed);
    PRS(report_continue_unperturbed);
    PRS(report_end_unperturbed);
    PRS(report_pericenter_reflection);
    PRS(report_impending_multiple_status);
    PRS(report_zero_unpert_steps);
    PRS(report_multiple);
    PRS(unpert_report_level);
    PRS(end_unpert_report_level);
    PRS(multiple_report_level);

    PRS(tree);
    PRS(tree_level);

    PRS(ev_function_id);
    PRS(ev);
    PRS(grape);
    PRS(grape_level);
    PRS(timestep_check);
    PRS(correct);
    PRS(slow_perturbed);

    PRS(kira_ev);

    PRS(slow);
    PRS(slow_level);

    PRS(report_stellar_evolution);
    PRS(report_stellar_mass_loss);
    PRS(report_binary_mass_loss);

    PRS(report_adjust_perturbed_list);

    PRS(report_correct_perturber_list);
}

#endif
