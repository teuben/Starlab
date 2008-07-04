
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Functions associated with kira counters.  Special treatment is needed
// here because we are largely circumventing the standard I/O mechanisms.
//
// Externally visible functions:
//
//	void initialize_counters_from_log
//	void write_counters_to_log
//	void print_counters

#include "hdyn.h"
#include "kira_timing.h"

// For initialize_counters_from_log:

local void getq(story* s, const char* l, int &i)	    {i = getiq (s, l);}
local void getq(story* s, const char* l, unsigned long &i)  {i = getulq(s, l);}
local void getq(story* s, const char* l, unsigned long long &i)
							    {i = getullq(s, l);}
local void getq(story* s, const char* l, real &r)	    {r = getrq (s, l);}
local void getq(story* s, const char* l, vec &v)	    {v = getvq (s, l);}

#  define GETLOG(x) if (find_qmatch(b->get_log_story(), #x)) \
		        getq(b->get_log_story(), #x, \
			     b->get_kira_counters()->x); \
		    else \
		        b->get_kira_counters()->x = 0;

// For write_counters_to_log:

local void putq(story* s, const char* l, int i)		    {putiq (s, l, i);}
local void putq(story* s, const char* l, unsigned long i)   {putiq(s, l, i);}
local void putq(story* s, const char* l, unsigned long long i)
							    {putiq(s, l, i);}
local void putq(story* s, const char* l, real r)	    {putrq (s, l, r);}
local void putq(story* s, const char* l, vec v)		    {putvq (s, l, v);}

#define PUTLOG(x) putq(b->get_log_story(), #x, b->get_kira_counters()->x);

// For print_counters:

#define PRLOGC(x) cerr << #x << " = " << b->get_kira_counters()->x << ",  "
#define PRLOGL(x) cerr << #x << " = " << kc->x; \
		  if (kc_prev) cerr << "  (" <<  kc->x - kc_prev->x << ")"; \
		  cerr << endl

void initialize_counters_from_log(hdyn* b)
{
    // Counters may be saved in the root log story.

    b->get_kira_counters()->cpu_time = 0;
    GETLOG(total_cpu_time);

    GETLOG(cpu_time_predict);
    GETLOG(cpu_time_top_level_force);
    GETLOG(cpu_time_low_level_force);
    GETLOG(cpu_time_external_force);
    GETLOG(cpu_time_total_force);
    GETLOG(cpu_time_correct);
    GETLOG(cpu_time_unperturbed);
    GETLOG(cpu_time_final_step);
    GETLOG(cpu_time_tree_check);
    GETLOG(cpu_time_integrate);
    GETLOG(cpu_time_other);

    GETLOG(step_top_single);
    GETLOG(step_top_cm);
    GETLOG(step_low_level);
    GETLOG(force_correction);
    GETLOG(limit_timestep_correction);

    GETLOG(pert_step);
    GETLOG(pert_with_list);
    GETLOG(pert_without_list);
    GETLOG(perturber_overflow);
    GETLOG(neighbor_sync);

    GETLOG(full_unpert_step);
    GETLOG(full_unpert_orbit);

    // Latter number has been known to overflow (bug?), apparently
    // corrupting the data structures.  Check and reset if this
    // occurs.  Apply an arbitrary limit of 2^62.  See also
    // hdyn_unpert.C.

    if (b->get_kira_counters()->full_unpert_orbit > ((step_t) 1)<<62) {
	cerr << "kira_counters: overflow in counter 'full_unpert_orbit': "
	     << "resetting to 0" << endl;
	b->get_kira_counters()->full_unpert_orbit = 0;
    }

    GETLOG(partial_unpert_orbit);

    GETLOG(tree_change);
    GETLOG(top_level_combine);
    GETLOG(low_level_combine);
    GETLOG(top_level_split);

    GETLOG(leaf_merge);
    GETLOG(total_kick);

    // Special case...

    int nss = 0;

    if (find_qmatch(b->get_log_story(), "n_step_slow"))
	getiq(b->get_log_story(), "n_step_slow", nss);

    if (nss > 0 && find_qmatch(b->get_log_story(), "step_slow"))
	getia(b->get_log_story(), "step_slow",
	      b->get_kira_counters()->step_slow, nss);

    GETLOG(dm_escapers);
    GETLOG(dm_massloss);

    GETLOG(initial_etot);
    GETLOG(de_total);
    GETLOG(de_merge);
    GETLOG(de_tidal_diss);
    GETLOG(de_massloss);
    GETLOG(de_kick);
}

void write_counters_to_log(hdyn* b)
{
    // cerr<<"in write_counters_to_log"<<endl;

    PUTLOG(cpu_time);
    PUTLOG(total_cpu_time);

    PUTLOG(cpu_time_predict);
    PUTLOG(cpu_time_top_level_force);
    PUTLOG(cpu_time_low_level_force);
    PUTLOG(cpu_time_external_force);
    PUTLOG(cpu_time_total_force);
    PUTLOG(cpu_time_correct);
    PUTLOG(cpu_time_unperturbed);
    PUTLOG(cpu_time_final_step);
    PUTLOG(cpu_time_tree_check);
    PUTLOG(cpu_time_integrate);
    PUTLOG(cpu_time_other);

    PUTLOG(step_top_single);
    PUTLOG(step_top_cm);
    PUTLOG(step_low_level);
    PUTLOG(force_correction);
    PUTLOG(limit_timestep_correction);

    PUTLOG(pert_step);
    PUTLOG(pert_with_list);
    PUTLOG(pert_without_list);
    PUTLOG(perturber_overflow);
    PUTLOG(neighbor_sync);

    PUTLOG(full_unpert_step);
    PUTLOG(full_unpert_orbit);
    PUTLOG(partial_unpert_orbit);

    PUTLOG(tree_change);
    PUTLOG(top_level_combine);
    PUTLOG(low_level_combine);
    PUTLOG(top_level_split);

    // Special case...

    int nss = SLOW_BINS;
    while (nss > 0 && b->get_kira_counters()->step_slow[nss-1] <= 0) nss--;

    if (nss > 0) {
	putiq(b->get_log_story(), "n_step_slow", nss);
	putia(b->get_log_story(), "step_slow",
	      b->get_kira_counters()->step_slow, nss);
    }

    PUTLOG(leaf_merge);
    PUTLOG(total_kick);

    PUTLOG(step_stellar);
    PUTLOG(step_dmslow);
    PUTLOG(step_dmfast);
    PUTLOG(step_correct);

    PUTLOG(dm_escapers);
    PUTLOG(dm_massloss);

    PUTLOG(initial_etot);
    PUTLOG(de_total);
    PUTLOG(de_merge);
    PUTLOG(de_tidal_diss);
    PUTLOG(de_massloss);
    PUTLOG(de_kick);
}

void print_counters(kira_counters* kc,
		    bool print_all,		// default = false
		    kira_counters* kc_prev)	// default = NULL
{
    cerr << "\n  Counters: \n";

    cerr << "    CPU time = " << kc->cpu_time;
    if (kc_prev)
	cerr << "  (delta = " << kc->cpu_time - kc_prev->cpu_time << ")";
    cerr << endl;

    if (print_all) {

#ifdef CPU_COUNTERS
	PRI(4); PRLOGL(cpu_time_predict);
	PRI(4); PRLOGL(cpu_time_top_level_force);
	PRI(4); PRLOGL(cpu_time_low_level_force);
	PRI(4); PRLOGL(cpu_time_external_force);
	PRI(4); PRLOGL(cpu_time_total_force);
	PRI(4); PRLOGL(cpu_time_correct);
	PRI(4); PRLOGL(cpu_time_unperturbed);
	PRI(4); PRLOGL(cpu_time_final_step);
	PRI(4); PRLOGL(cpu_time_tree_check);
	PRI(4); PRLOGL(cpu_time_integrate);
	PRI(4); PRLOGL(cpu_time_other);
	cerr << endl;
#endif

	PRI(4); PRLOGL(step_top_single);
	PRI(4); PRLOGL(step_top_cm);
	PRI(4); PRLOGL(step_low_level);

	PRI(4); PRLOGL(pert_step);
	PRI(4); PRLOGL(pert_with_list);
	PRI(4); PRLOGL(pert_without_list);
	PRI(4); PRLOGL(perturber_overflow);
	PRI(4); PRLOGL(neighbor_sync);

	PRI(4); PRLOGL(force_correction);
	PRI(4); PRLOGL(limit_timestep_correction);

	PRI(4); PRLOGL(full_unpert_step);
	PRI(4); PRLOGL(full_unpert_orbit);
	PRI(4); PRLOGL(partial_unpert_orbit);
	cerr << endl;

	PRI(4); PRLOGL(tree_change);
	PRI(4); PRLOGL(top_level_combine);
	PRI(4); PRLOGL(low_level_combine);
	PRI(4); PRLOGL(top_level_split);

	// Special case...

	int nss = SLOW_BINS;
	while (nss > 0 && kc->step_slow[nss-1] <= 0) nss--;

	cerr << "    step_slow =";

	if (nss > 0) {
	    for (int i = 0; i < nss; i++)
		cerr << " " << kc->step_slow[i];
	    cerr << endl;
	    if (kc_prev) {
		cerr << "              (";
		for (int i = 0; i < nss; i++)
		    cerr << " " << kc->step_slow[i] - kc_prev->step_slow[i];
		cerr << " )" << endl;
	    }
	} else
	    cerr << " 0" << endl;

	PRI(4); PRLOGL(leaf_merge);
	PRI(4); PRLOGL(total_kick);
	PRI(4); PRLOGL(step_stellar);
	PRI(4); PRLOGL(step_dmslow);
	PRI(4); PRLOGL(step_dmfast);
	PRI(4); PRLOGL(step_correct);

	PRI(4); PRLOGL(dm_escapers);
	PRI(4); PRLOGL(dm_massloss);

	PRI(4); PRLOGL(de_total);
	PRI(4); PRLOGL(de_merge);
	PRI(4); PRLOGL(de_tidal_diss);
	PRI(4); PRLOGL(de_massloss);
	PRI(4); PRLOGL(de_kick);
    }

    if (kc_prev)
	*kc_prev = *kc;
}
