
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// The following macros control most of the diagnostic output
// coming from the dynamical integration routines in kira.
//
// Recommended defaults are indicated [in square brackets].
//
// This file is NO LONGER USED, as the functionality has been
// absorbed into hdyn::kira_diag.

//------------------------------------------------------------------------
// Top-level (kira.C):

#define DEBUG_MAIN				false		// [false]

//------------------------------------------------------------------------
// Unperturbed motion (hdyn_unpert.C):

#define REPORT_START_UNPERTURBED        	false		// [true]

#define REPORT_CONTINUE_UNPERTURBED		false		// [false]

//#define REPORT_CONTINUE_UNPERTURBED \
//	(streq(format_label(), "743"))
//#define REPORT_CONTINUE_UNPERTURBED \
//	(time >  68.4374 && (streq(format_label(), "308a") \
//	 		   || streq(format_label(), "21a")))
//#define REPORT_CONTINUE_UNPERTURBED \
//	(time > 1215.8348 && time < 1215.8351)

#define REPORT_END_UNPERTURBED			false		// [true]

#define REPORT_PERICENTER_REFLECTION		false		// [false]

//#define REPORT_PERICENTER_REFLECTION \
//	(streq(format_label(), "308a") \
//	 		   || streq(format_label(), "21a"))

#define REPORT_IMPENDING_MULTIPLE_STATUS	false		// [false]
#define REPORT_ZERO_UNPERT_STEPS		false		// [false]
#define REPORT_MULTIPLE				false		// [true]

#define UNPERT_REPORT_LEVEL			0		// [0]
#define END_UNPERT_REPORT_LEVEL			0		// [0]
#define MULTIPLE_REPORT_LEVEL			0		// [0]

#define UNPERT_FUNCTION_ID			false		// [false]

//#define UNPERT_FUNCTION_ID \
//	(time > 1215.8348 && time < 1215.8351)

// Note:  constructs such as the following are legal, as all of the
//	  above macros are only used in member functions:
//
// #define UNPERT_FUNCTION_ID \
//	(time > 99 && system_time > 100 && system_time < 101 \
//	 && streq(format_label(), "42"))

//------------------------------------------------------------------------
// Tree maintenance (hdyn_tree.C):

#define DEBUG_TREE				false		// [true]
#define TREE_DEBUG_LEVEL			1		// [0]

//------------------------------------------------------------------------
// Force evaluation (hdyn_ev.C):

#define DEBUG_EV				false		// [false]
#define EV_FUNCTION_ID				false		// [false]
#define HARP_DEBUG				true		// [true]
#define CORRECT_DEBUG				false		// [false]
#define SLOW_PERTURBED_DEBUG			false		// [false]

//------------------------------------------------------------------------
// Force evaluation (kira_ev.C):

#define DEBUG_KIRA_EV				false		// [false]

// Not really debugging, but define here for now...

#define USE_OLD_CORRECT_ACC_AND_JERK		false		// [false]

//------------------------------------------------------------------------
// Slow binary motion (hdyn_slow.C):

#define DEBUG_SLOW				true		// [true]
#define SLOW_DEBUG_LEVEL			0		// [0]

//------------------------------------------------------------------------
// Stellar evolution (kira_stellar.C):

#define REPORT_STELLAR_EVOLUTION		false		// [false]
#define REPORT_STELLAR_MASS_LOSS		false		// [false]
#define REPORT_BINARY_MASS_LOSS			false		// [false]

//------------------------------------------------------------------------
// Perturbed binary list (perturbed_list.C):

#define REPORT_ADJUST_PERTURBED_LIST		false		// [false]

//------------------------------------------------------------------------
// Perturber correction (correct_perturbers.C):

#define REPORT_CORRECT_PERTURBER_LIST		false		// [false]
