
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// Functions associated with log output, snapshots, etc.
//
// Externally visible functions:
//
//	void log_output
//	void snap_output
//	void check_runtime

#include "node.h" 
#include "double_star.h"
#include "main_sequence.h"
#include "dstar_to_kira.h"
#include "seba.h"

local void putq(story* s, char* l, int i)		{putiq (s, l, i);}
local void putq(story* s, char* l, unsigned long i)	{putulq(s, l, i);}
local void putq(story* s, char* l, real r)		{putrq (s, l, r);}
local void putq(story* s, char* l, vec v)		{putvq (s, l, v);}

#define PUTLOG(x) putq(sb->get_star_story(), #x, sb->get_seba_counters()->x);

#define PRLOGC(x) cerr << #x << " = " << sb->get_seba_counters()->x << ",  "
#define PRLOGL(x) cerr << #x << " = " << sbc->x; \
		  if (sbc_prev) cerr << "  (" <<  sbc->x - sbc_prev->x << ")"; \
		  cerr << endl

void print_counters(seba_counters* sbc, seba_counters* sbc_prev)
{
    cerr << "\n  Counters: \n";

    cerr << "    CPU time = " << sbc->cpu_time;
    if (sbc_prev)
	cerr << "  (delta = " << sbc->cpu_time - sbc_prev->cpu_time << ")";
    cerr << endl;

    PRI(4); PRLOGL(add_dstar);
    PRI(4); PRLOGL(del_dstar);

    PRI(4); PRLOGL(step_seba);
    PRI(4); PRLOGL(step_sstar);
    PRI(4); PRLOGL(step_dstar);

    PRI(4); PRLOGL(detached);
    PRI(4); PRLOGL(semi_detached);
    PRI(4); PRLOGL(contact);

    PRI(4); PRLOGL(dynamic);
    PRI(4); PRLOGL(thermal);
    PRI(4); PRLOGL(nuclear);
    PRI(4); PRLOGL(aml_driven);

    PRI(4); PRLOGL(supernovae);
    PRI(4); PRLOGL(first_rlof);
    PRI(4); PRLOGL(common_envelope);
    PRI(4); PRLOGL(spiral_in);
    PRI(4); PRLOGL(double_spiral_in);
    PRI(4); PRLOGL(mergers);

    PRI(4); PRLOGL(aml_mergers);
    PRI(4); PRLOGL(gwr_mergers);

    PRI(4); PRLOGL(recursive_overflow);

    if (sbc_prev)
	*sbc_prev = *sbc;
}

local void write_counters_to_log(starbase* sb)
{
    // cerr<<"in write_counters_to_log"<<endl;

    PUTLOG(cpu_time);

    PUTLOG(add_dstar);
    PUTLOG(del_dstar);
 
    PUTLOG(step_seba);
    PUTLOG(step_sstar);
    PUTLOG(step_dstar);

    PUTLOG(detached);
    PUTLOG(semi_detached);
    PUTLOG(contact);

    PUTLOG(dynamic);
    PUTLOG(thermal);
    PUTLOG(nuclear);
    PUTLOG(aml_driven);

    PUTLOG(supernovae);
    PUTLOG(first_rlof);
    PUTLOG(common_envelope);
    PUTLOG(spiral_in);
    PUTLOG(double_spiral_in);
    PUTLOG(mergers);

    PUTLOG(aml_mergers);
    PUTLOG(gwr_mergers);

    PUTLOG(recursive_overflow);

}



typedef  struct {
    real  dt;
    star*  s;
} dt_pair, *dt_pair_ptr;

local int compare_dt(const void * pi, const void * pj)	  // increasing dt
{
    if (((dt_pair_ptr) pi)->dt < ((dt_pair_ptr) pj)->dt)
        return -1;
    else if (((dt_pair_ptr)pi)->dt > ((dt_pair_ptr)pj)->dt)
        return +1;
    else
        return 0;
}

typedef  struct {
    real  count;
    star*  s;
} count_pair, *count_pair_ptr;

local int compare_steps(const void * pi, const void * pj) // decreasing count
{
    if (((count_pair_ptr) pi)->count > ((count_pair_ptr) pj)->count)
        return -1;
    else if (((count_pair_ptr)pi)->count < ((count_pair_ptr)pj)->count)
        return +1;
    else
        return 0;
}
