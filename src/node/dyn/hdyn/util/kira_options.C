
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// kira_options:  Starlab kira/options-specific functions.

#include "hdyn.h"

#ifndef TOOLBOX

kira_options::kira_options() {		    // (Steve's 6/00
					    //  preferences...)

    print_xreal					= true;

    perturber_criterion				= 2;	// = hybrid

    optimize_scheduling				= true;
    optimize_block				= false;
    allow_unperturbed				= true;
    allow_multiples				= false;

    min_unpert_steps				= 5;
    full_merge_tolerance			= 1.e4;
    relax_factor				= 10;
    partial_merge_tolerance			= 0.01;
    full_merge_tol_for_close_binary		= 1.e-4;
    multiple_merge_tolerance			= 1.e6;
    unconditional_stable_fac			= 5;	// fairly conservative
    partial_stable_fac				= 30;

    use_aarseth_criterion			= true;
    aarseth_stable_fac				= 2.8;	// don't change!

    close_criterion				= 2;	// = force

    allow_keplstep				= true;

    use_old_correct_acc_and_jerk		= false;

    grape_check_count				= 15000;
    grape_max_cpu				= 15;	// seconds
    grape_last_cpu				= 0;

    use_perturbed_list				= true;

}

#define PRS(x) s << "    " << #x << " = " << x << endl

void kira_options::print(ostream &s)
{
    s << endl << "kira_options settings:" << endl;

    PRS(print_xreal);

    PRS(perturber_criterion);
    PRS(optimize_scheduling);
    PRS(optimize_block);
    PRS(allow_unperturbed);
    PRS(allow_multiples);

    PRS(min_unpert_steps);
    PRS(full_merge_tolerance);
    PRS(relax_factor);
    PRS(partial_merge_tolerance);
    PRS(full_merge_tol_for_close_binary);
    PRS(multiple_merge_tolerance);
    PRS(unconditional_stable_fac);
    PRS(partial_stable_fac);

    PRS(use_aarseth_criterion);
    PRS(aarseth_stable_fac);

    PRS(close_criterion);

    PRS(allow_keplstep);

    PRS(use_old_correct_acc_and_jerk);

    PRS(grape_check_count);
    PRS(grape_max_cpu);
    PRS(grape_last_cpu);

    PRS(use_perturbed_list);

}

#endif
