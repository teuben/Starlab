
#include "sigma_MPI.h"

void initialize_processors(int nproc, int debug){}	// Dummy functions
void terminate_processors(int debug){}

// multiscatter:  Perform a specified number of scattering experiments
//		  in a given impact parameter range, accumulating statistics
//		  as we go.  Return the total number of hits in this range.


local void pp(sdyn* b, ostream & s, int level = 0) {

    s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    if(b != b->get_root()) {
      b->pretty_print_node(s);
      s << " \t"<< b->get_mass() << " \t"
		<< b->get_pos() << " (r= " << abs(b->get_pos()) << ")   " 
		<< b->get_vel() << " (v= " << abs(b->get_vel()) << ")" << endl;
      //	<< "r= " << abs(b->get_pos()) << "    " 
      //	<< "v= " << abs(b->get_vel()) << endl;
    }

    for (sdyn * daughter = b->get_oldest_daughter();
	 daughter != NULL;
	 daughter = daughter->get_younger_sister())
	pp(daughter, s, level + 1);	
}

// Adjusted for n-body scattering experiments.

int  multiscatter(sigma_out &out, 
		  sigma_input &input, 
		  MPI_Datatype inputtype,
		  scatter_exp &experiment, 
		  MPI_Datatype scatter_exp_type,
		  real &cpu_save, int& scatt_total, real& cpu_total) { 


  real eta = input.eta;
  real delta_t = input.delta_t;
  real dt_out = input.dt_out;
  real dt_snap = input.dt_snap;
  real ttf = input.tidal_tol_factor;
  real snap_cube_size = input.snap_cube_size;
  real cpu_time_check = input.cpu_time_check;  
  int debug = input.debug;
  int scatter_summary_flag = false;

  int total_hits = 0;

  int single_result = 0;

  real cpu_init = 0;
  total_hits += master_process(out, input, inputtype,
			       experiment, scatter_exp_type);

  cpu_total += cpu_time() - cpu_init;

  scatt_total++;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Hand-holding: check the CPU time.  Note that the printed
	// CPU time is the time since this routine was entered.

  //      real cpu_init = 0;
	if (cpu_time() - cpu_save > cpu_time_check) {
	    cpu_save = cpu_time();
	    cerr << "\nsigma:  CPU time = " << cpu_save - cpu_init
		 << ",  " << out.total_trials + 1 << " trials,  "
		 << out.n_hit_tot + single_result
		 << " hits,  i_max = " << out.i_max
		 << endl << flush;
	}

	if (abs(input.debug) > 2) {
	    cerr.precision(6);
	    cerr << "single_scatter: rho_max^2 = " << input.rho_sq_max
		 << " ; returning with n_hit = " << single_result
		 << " and n_hit_tot = " << out.n_hit_tot + single_result
		 << endl;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Accumulate statistics on the results.
	//	single_scatter_stats(&experiment, out); //, acc_stats);

    return total_hits;

}



