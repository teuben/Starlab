
// Turn on/off various pieces of timing code.

// In calculate_acc_and_jerk_for_list (kira_ev.C):  Set TIME_INTERNAL
// to perform timing of force-evaluation functions.

#define TIME_INTERNAL
#undef  TIME_INTERNAL

// In integrate_list (kira.C):  Define TIME_LIST to time various force
// calculations.

#define TIME_LIST
#undef  TIME_LIST

// In integrate_list (kira.C) and calculate_acc_and_jerk_for_list (kira_ev.C):
// Define CPU_COUNTERS to time various parts of the integration procedure.

#define CPU_COUNTERS
#undef  CPU_COUNTERS
