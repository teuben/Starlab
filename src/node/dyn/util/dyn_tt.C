//
// dyn_tt.C: dyn-specific node-handling functions.
//

#include "dyn.h"

void dyn::null_pointers()
{
    // Clear all pointers (don't touch what they point to)
    // -- for use in cleaning up temporary nodes...  Careful!!

    kep = NULL;

    node::null_pointers();
}

bool dyn::nn_stats(real energy_cutoff, real kT,
		   vector center, bool verbose,
		   bool long_binary_output,	// default = true
		   int which)			// default = 0
{
    // Kludgy dummy function to enable hdyn::nn_stats().

    // Print header information on first entry.

    if (verbose && long_binary_output
	&& is_top_level_node() && elder_sister == NULL && which == 0)
        cerr << "\n  (No NN information available for dyn class)" << endl;

    return true;        // Avoid "(none)" output in sys_stats.
}

real dyn::print_pert(bool long_binary_output,	// default = true
		     int indent)		// default = BIN_INDENT
{
    // Similarly kludgy function to enable hdyn::print_pert().

    if (!long_binary_output) cerr << endl;
    return 0;
}

real dyn::get_radius()
{
    if (oldest_daughter == NULL) return 0;

    real rmax2 = 0;

    for (dyn* b = this; b != NULL; b = (dyn*) b->next_node(this))
	if (b->get_oldest_daughter() == NULL) {
	    
	    vector dx = b->get_pos();
	    real r2 = dx*dx;
	    if (r2 > rmax2) rmax2 = r2;
	}
    return sqrt(rmax2);
}

vector something_relative_to_root(dyn* bi,
				  dyn_VMF_ptr get_something)
{
#ifdef DEBUG
    cerr << "something_relative_to_ancestor\n";
#endif    
    vector  d_something = 0.0;
    for(dyn * b = bi;b->get_parent() != NULL; b = b->get_parent()){
	if(b == NULL){
	    cerr << "something_relative_to_root: Error, bj is not the "
	         << "ancestor of bi\n";
	    exit(1);
	}
	d_something += (b->*get_something)();
    }
    return d_something;
}
