
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// kira_options:  Starlab kira/options-specific functions.

#include "hdyn.h"
#include "../evolve/kira_defaults.h"

#ifndef TOOLBOX

kira_options::kira_options() {

    print_xreal		=	DEFAULT_PRINT_XREAL;

    perturber_criterion	=	DEFAULT_PERTURBER_CRITERION;

    optimize_scheduling	=	DEFAULT_OPTIMIZE_SCHEDULING;
    optimize_block	=	DEFAULT_OPTIMIZE_BLOCK;
    allow_unperturbed	=	DEFAULT_ALLOW_UNPERTURBED;
    allow_multiples	=	DEFAULT_ALLOW_MULTIPLES;

    min_unpert_steps	=	DEFAULT_MIN_UNPERT_STEPS;
    full_merge_tolerance =	DEFAULT_FULL_MERGE_TOLERANCE;
    relax_factor	=	DEFAULT_RELAX_FACTOR;
    partial_merge_factor =	DEFAULT_PARTIAL_MERGE_FACTOR;
    full_merge_tol_for_close_binary =
				DEFAULT_FULL_MERGE_TOL_FOR_CLOSE_BINARY;
    multiple_merge_tolerance =	DEFAULT_MULTIPLE_MERGE_TOLERANCE;
    unconditional_stable_fac =	DEFAULT_UNCONDITIONAL_STABLE_FAC;
    partial_stable_fac	=	DEFAULT_PARTIAL_STABLE_FAC;

    use_aarseth_criterion =	DEFAULT_USE_AARSETH_CRITERION;
    aarseth_stable_fac	=	DEFAULT_AARSETH_STABLE_FAC;

    close_criterion	=	DEFAULT_CLOSE_CRITERION;

    allow_keplstep	=	DEFAULT_ALLOW_KEPLSTEP;

    use_old_correct_acc_and_jerk =
				DEFAULT_USE_OLD_CORRECT_ACC_AND_JERK;

    grape_check_count	=	DEFAULT_GRAPE_CHECK_COUNT;
    grape_max_cpu	=	DEFAULT_GRAPE_MAX_CPU;
    grape_last_cpu	=	DEFAULT_GRAPE_LAST_CPU;

    use_perturbed_list	=	DEFAULT_USE_PERTURBED_LIST;

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
    PRS(partial_merge_factor);
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
