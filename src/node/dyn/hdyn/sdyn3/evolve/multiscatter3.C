
#include "sigma3.h"

void initialize_processors(int nproc, int debug){}	// Dummy functions
void terminate_processors(int debug){}

// multiscatter3: Perform a specified number of scattering experiments
//		  in a given impact parameter range, accumulating statistics
//		  as we go.  Return the total number of hits in this range.

int multiscatter3(scatter_profile & prof, sigma_out & out,
		  real rho_sq_min, real rho_sq_max, int rho_zone,
		  real dt_snap, real snap_cube_size,
		  real cpu_time_check, real cpu_init, real &cpu_save,
		  int& scatt_total, real& cpu_total, stat_fp acc_stats,
		  int debug, int scatter_summary_flag)
{
    int total_hits = 0;

    for (int trial = 0; trial < out.trials_per_zone; trial++) {
	
	int scatter_summary = scatter_summary_flag;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Get the next scattering.

	// The function single_scatter_init will randomize the angles
	// in the same way as the standalone tool scatter3.  Save the
	// random counter n_rand in case we want a one-line summary below.

	initial_state3 init;
	int n_rand;

	single_scatter_init(prof, rho_sq_min, rho_sq_max, init, n_rand,
			    scatter_summary_flag, dt_snap, snap_cube_size);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Perform the scattering, keeping track of the outcome.

	intermediate_state3 inter;
	final_state3 final;

	real cpu_scatter = cpu_time();
	int single_result = single_scatter(init, inter, final,
					   cpu_time_check,
					   dt_snap, snap_cube_size);
	total_hits += single_result;
	cpu_scatter = cpu_time() - cpu_scatter;

	scatt_total++;
	cpu_total += cpu_scatter;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Hand-holding: check the CPU time.  Note that the printed
	// CPU time is the time since this routine was entered.

	if (cpu_time() - cpu_save > cpu_time_check) {
	    cpu_save = cpu_time();
	    cerr << "\nsigma3:  CPU time = " << cpu_save - cpu_init
		 << ",  " << out.total_trials + 1 << " trials,  "
		 << out.n_hit_tot + single_result
		 << " hits,  i_max = " << out.i_max
		 << endl << flush;
	}

	if (abs(debug) > 2) {
	    int p = cerr.precision(STD_PRECISION);
	    cerr << "single_scatter: rho_max^2 = " << rho_sq_max
		 << " ; returning with n_hit = " << single_result
		 << " and n_hit_tot = " << out.n_hit_tot + single_result
		 << endl;
	    cerr.precision(p);
	}

	if (final.descriptor == error) {
	    if (scatter_summary == 0)
		summarize_scattering_initial(init, n_rand,
					     dt_snap, snap_cube_size);
	    scatter_summary = 2;
	}

	if (scatter_summary > 0)
	    summarize_scattering_final(inter, final,
				       scatter_summary, cpu_scatter);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Accumulate statistics on the results.

	single_scatter_stats(prof, init, inter, final, rho_zone,
			     out, acc_stats, single_result);
    }

    return total_hits;

}

