
// Compiled-in kira debugging options (for efficiency).  Still in progress.
// Various new options should be integrated with existing code.

//-------------------------------------------------------------------------
//
// General switch for non-essential debugging/diagnostic checks:

#define KIRA_DEBUG
#undef  KIRA_DEBUG

//-------------------------------------------------------------------------
//
// Track progress through the code (main integration loop and elsewhere,
// as defined below), after time T_DEBUG, if defined.

// Main compile-time switch is the definition of T_DEBUG, which also
// sets a start time.

#define T_DEBUG 0
#undef  T_DEBUG

// Debugging time range is [T_DEBUG, T_DEBUG_END].

#define T_DEBUG_END VERY_LARGE_NUMBER

// Debugging level may be used in some functions:

#define T_DEBUG_LEVEL 0

// Control which functions actually turn on debugging (we undef T_DEBUG
// in the function if the function name is undefined here).  The following
// currently include this file:

#define IN_DEBUG_RANGE(t) ((t) >= T_DEBUG && (t) <= T_DEBUG_END)

#undef	T_DEBUG_hdyn_ev
#define	T_DEBUG_hdyn_grape6
#define	T_DEBUG_kira
#define	T_DEBUG_kira_ev
