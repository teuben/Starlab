
// Default parameter settings for parameters we might want to change.

// Used in util/hdyn_io.C and evolve/kira_init.C:

#define	DEFAULT_ETA				0.1
#define	DEFAULT_EPS				0.0
#define	DEFAULT_GAMMA				1.e-7
#define	DEFAULT_D_MIN_FAC			0.25
#define	DEFAULT_LAG_FACTOR			2.5

#define	DEFAULT_DT				10.0
#define	DEFAULT_DT_REINIT			1.0
#define	DEFAULT_DT_FULLDUMP			1.0
#define	DEFAULT_DT_LOG				1.0
#define	DEFAULT_DT_SSTAR			0.015625

// Used in util/kira_options.C and evolve/kira_runtime.C:

#define DEFAULT_PRINT_XREAL			true

#define DEFAULT_PERTURBER_CRITERION		2	// hybrid

#define DEFAULT_OPTIMIZE_SCHEDULING		true
#define DEFAULT_OPTIMIZE_BLOCK			true
#define DEFAULT_ALLOW_UNPERTURBED		true
#define DEFAULT_ALLOW_MULTIPLES			false

#define DEFAULT_MIN_UNPERT_STEPS		5
#define DEFAULT_FULL_MERGE_TOLERANCE		1.e-8	// new 1/14/02
#define DEFAULT_RELAX_FACTOR			10.0
#define DEFAULT_PARTIAL_MERGE_FACTOR		0.01
#define DEFAULT_FULL_MERGE_TOL_FOR_CLOSE_BINARY	1.e-4
#define DEFAULT_MULTIPLE_MERGE_TOLERANCE	1.e-8	// new 1/14/02
#define DEFAULT_UNCONDITIONAL_STABLE_FAC	5.0	// conservative
#define DEFAULT_PARTIAL_STABLE_FAC		30.0

#define DEFAULT_USE_AARSETH_CRITERION		true
#define DEFAULT_AARSETH_STABLE_FAC		2.8	//don't change!

#define DEFAULT_CLOSE_CRITERION			2	// force

#define DEFAULT_ALLOW_KEPLSTEP			true

#define DEFAULT_USE_OLD_CORRECT_ACC_AND_JERK	false

#define DEFAULT_GRAPE_CHECK_COUNT		15000
#define DEFAULT_GRAPE_MAX_CPU			15.0	// seconds
#define DEFAULT_GRAPE_LAST_CPU			0.0

#define DEFAULT_USE_PERTURBED_LIST		true

// Used in util/hdyn_io.C and evolve/kira_runtime.C:

#define DEFAULT_MAX_SLOW_FACTOR			1	// suppress
#define DEFAULT_MAX_SLOW_PERTURBATION		1.e-3
#define DEFAULT_MAX_SLOW_PERTURBATION_SQ	1.e-6	// keep consistent!
