
// Counters for kira:

typedef unsigned long step_t;

#define SLOW_BINS 20

class kira_counters {

    public:

	real cpu_time;			// CPU time since this process started
	real total_cpu_time;		// total CPU time for the calculation

	// These should eventually be #ifdefed to avoid
	// taking too much CPU time.

	real cpu_time_predict;		// updated in kira_ev.C
	real cpu_time_top_level_force;	// updated in kira_ev.C
	real cpu_time_low_level_force;	// updated in kira_ev.C
	real cpu_time_external_force;	// updated in kira_ev.C
	real cpu_time_total_force;	// updated in kira.C
	real cpu_time_correct;		// updated in kira.C
	real cpu_time_unperturbed;	// updated in kira.C
	real cpu_time_final_step;	// updated in kira.C
	real cpu_time_tree_check;	// updated in kira.C
	real cpu_time_integrate;	// updated in kira.C
	real cpu_time_other;		// updated in kira.C

	step_t step_top_single;		// top-level single steps
	step_t step_top_cm;		// top-level CM steps
	step_t step_low_level;		// low-level steps
	step_t force_correction;	// force corrections

	step_t pert_step;		// low-level perturbed steps
	step_t pert_with_list;		// low-level *interparticle force*
    					// calculations using perturber list
	step_t pert_without_list;	// low-level *force* calculations
					// without perturber list
	step_t perturber_overflow;	// perturber list overflows
	step_t neighbor_sync;		// synchronize nearest neighbor

	step_t full_unpert_step;	// full unperturbed steps
	step_t full_unpert_orbit;	// full unperturbed orbits
	step_t partial_unpert_orbit;	// partial unperturbed orbits

	step_t tree_change;		// changes in tree structure

	step_t top_level_combine;	// top-level combines within tree
	step_t low_level_combine;	// low-level combines within tree
	step_t top_level_split;		// top-level splits within tree

	step_t leaf_merge;		// physical mergers between stars
	step_t total_kick;		// number of kicks
	step_t step_stellar;		// number of stellar evolution steps
	step_t step_dmslow;		// number of slow mass loss events
	step_t step_dmfast;		// number of fast mass loss events
	step_t step_correct;		// number stellar evolution corrections

	step_t step_slow[SLOW_BINS];	// slow binary steps [slowdown factor]

	real dm_escapers;		// total mass loss in escapers (>0)
	real dm_massloss;		// total mass loss by stellar evolution

	real initial_etot;		// total energy at t = 0
	real de_total;			// total energy change (subdivided:)

	real de_merge;			// energy change by mergers
	real de_tidal_diss;		// energy change by tidal dissipation
	real de_massloss;		// energy change by stellar mass loss
	real de_kick;			// energy change by supernova kicks

	kira_counters () {

	    cpu_time                 =
	    total_cpu_time           =
	    cpu_time_predict	     =
	    cpu_time_top_level_force =
	    cpu_time_low_level_force =
	    cpu_time_external_force  =
	    cpu_time_total_force     =
	    cpu_time_correct	     =
	    cpu_time_unperturbed     =
	    cpu_time_final_step	     =
	    cpu_time_tree_check	     =
	    cpu_time_integrate	     =
	    cpu_time_other	     = 0;

	    step_top_single          =
	    step_top_cm              = 
	    step_low_level           = 
	    force_correction         = 
	    pert_step                = 
	    pert_with_list           = 
	    pert_without_list        = 
	    perturber_overflow       = 
	    neighbor_sync	     =
	    full_unpert_step         = 
	    full_unpert_orbit        = 
	    partial_unpert_orbit     = 
	    tree_change              = 
	    top_level_combine        = 
	    low_level_combine        = 
	    top_level_split          = 
	    leaf_merge               =  
	    step_stellar             =
	    step_dmslow              =
	    step_dmfast              =
	    step_correct             =
	    total_kick               = 0;
	    
	    for (int i = 0; i < SLOW_BINS; i++) step_slow[i] = 0;

	    dm_escapers              =
	    dm_massloss              =
	    de_kick                  =
	    de_massloss              =
	    de_merge                 =
	    de_tidal_diss            =
	    de_total                 =
	    initial_etot             = 0;
	}

	void inc_slow(int kappa) {
	    if (kappa > 1) {
		int k = 0, k2 = 2;
		while (k2 < kappa && k < SLOW_BINS-1) {k++; k2 *= 2;}
		step_slow[k]++;
	    }
	}
};

// Runtime options and switches:

class kira_options {

    public:

	// Used in several references to put_node:

	bool print_xreal;

	// Used in hdyn_inline.C:

	int  perturber_criterion;

	// Used in hdyn_unpert.C (see source for description):

	bool optimize_scheduling;
	bool optimize_block;
	bool allow_unperturbed;
	bool allow_multiples;

	int  min_unpert_steps;
	real full_merge_tolerance;
	real relax_factor;
	real partial_merge_factor;
	real full_merge_tol_for_close_binary;
	real multiple_merge_tolerance;
	real unconditional_stable_fac;
	bool use_aarseth_criterion;
	real aarseth_stable_fac;
	real partial_stable_fac;

	// Used in hdyn_tree.C:

	int  close_criterion;

	// Used in hdyn_ev.C:

	bool allow_keplstep;

	// Used in kira_ev.C:

    	bool use_old_correct_acc_and_jerk;

	// Used in hdyn_grape4/6.C:

	int  grape_check_count;
	real grape_max_cpu;
	real grape_last_cpu;	// (not really an option, but
				//  convenient to keep it here)

	int  grape_coll_freq;
	int  grape_pert_freq;

	// Used in perturbed_list.C:

	bool use_perturbed_list;

	// In kira_options.C:

	kira_options();
	void print(ostream & s = cerr);
};

// Diagnostic/debugging flags (may be changed during a run):
// Warning: these can produce a *lot* of output!
//
// Currently, only the following diag elements allow the possibility
// of diag->check_diag():
//
//	ev_function_id
//	ev
//	timestep_check

class kira_diag {

    private:

	char *name;
	real t1, t2;

    public:

	// Generic debugging (general -- for short-term use only):

	int kira_runtime_flag;

	// Used in kira.C:

	bool kira_main;
	bool check_heartbeat;
	int  n_check_heartbeat;
	int  n_check_runtime;

	// Used in hdyn_unpert.C:

	bool unpert_function_id;
	bool report_start_unperturbed;
	bool report_continue_unperturbed;
	bool report_end_unperturbed;
	bool report_pericenter_reflection;
	bool report_impending_multiple_status;
	bool report_zero_unpert_steps;
	bool report_multiple;
	int  unpert_report_level;
	int  end_unpert_report_level;
	int  multiple_report_level;

	// Used in hdyn_tree.C and hdyn_merge.C:

	bool tree;
	int tree_level;

	// Used in hdyn_ev.C:

	bool ev_function_id;
	bool ev;
	bool grape;
	int  grape_level;			// (also in kira.C)
	bool timestep_check;
	bool correct;
	bool slow_perturbed;			// (also in hdyn_slow.C)

	// Used in kira_ev.C:

	bool kira_ev;

	// Used in hdyn_slow.C:

	bool slow;
	int  slow_level;

	// Used in kira_stellar.C:

	bool report_stellar_evolution;
	bool report_stellar_mass_loss;
	bool report_binary_mass_loss;

	// Used in perturbed_list.C:

	bool report_adjust_perturbed_list;

	// Used in correct_perturbers.C:

	bool report_correct_perturber_list;

	// Used in kira_encounter.C:

	int report_close_encounters;

	// Passed to kepler.C:

	bool report_kepler_trig_error;

	void set_name(char *n)			{name = n;}
	char *get_name()			{return name;}
	void clear_name()			{if (name) delete [] name;
						 name = NULL;}

	void set_t1(real t)			{t1 = t;}
	void set_t2(real t)			{t2 = t;}
	void set_range(real tt1, real tt2)	{t1 = tt1, t2 = tt2;}
	real get_t1()				{return t1;}
	real get_t2()				{return t2;}

	// In kira_diag.C:

	kira_diag();
	bool check_diag(hdyn *b);
	void print(ostream & s = cerr);
};
