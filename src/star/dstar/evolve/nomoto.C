//
//  binev.C: evolves a single binary borrowing the tree structure from node 
//
//	       Simon Portegies Zwart, Nov 1996
//

#include "node.h"
#include "double_star.h"
#include "main_sequence.h"
//#include "dstar_to_dyn.h"

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  binev  --
 *-----------------------------------------------------------------------------
 */
local bool  evolve_binary(node * bi,
                          real start_time, real end_time) { 

  double_star* ds = dynamic_cast(double_star*, 
				 bi->get_starbase());
  
  //		Setup star from input data.
  real dt, time=start_time;
  ds->evolve_element(time);
  
  ds->dump("init.dat", true);

  if (!bi->is_root() &&
      bi->get_parent()->is_root()) 

    do {

      //dt = bi->get_starbase()->get_evolve_timestep();
      dt = 100;
      time = min(time+dt, end_time);

      PRC(bi->get_starbase()->get_evolve_timestep());PRC(time);
      PRC(end_time);PRL(ds->get_bin_type());
      PRL(((double_star*)bi->get_starbase())->get_donor_timescale());
      
      ds->evolve_element(time);

      if (ds->get_bin_type() == Merged || 
	  ds->get_bin_type() == Disrupted)
	return false;

      if (ds->get_primary()->remnant() || ds->get_secondary()->remnant())
	return false;

    }
    while (time<end_time);
    
  return true;

}

void main(int argc, char ** argv) {

    bool e_flag = false;
    bool R_flag = true;

    bool reandom_initialization = false;
    int  n = 1;
    int  n_init = 0;
    int id=1;
    stellar_type type = Main_Sequence;
    binary_type bin_type = Detached;
    real binary_fraction = 1.0;

    real m_tot = 1;
    real r_hm = 1;
    real t_hc = 1;

    real m_min = 0.5;
    real m_max = 100;
    real lm_min = log10(m_min);
    real lm_max = log10(m_max);
    real a_min = 10;
    real a_max = 1.e+4;
    real la_min = log10(a_min);
    real la_max = log10(a_max);
    real q_min = 0;
    real q_max = 1;
    real sma    = 138;
    real ecc    = 0.0;
    real m_prim = 13.1;
    real m_sec  =  9.8;

    real start_time = 0;
    real end_time   = 35;

    int input_seed=0, actual_seed;
    char seedlog[64];
    extern char *poptarg;
    int c;
    char* param_string = "A:a:e:M:m:N:n:Q;q:Rs:T:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case 'A': a_max = atof(poptarg);
	              la_max = log10(a_max);
                      break;
            case 'a': a_min = atof(poptarg);
	              la_min = log10(a_min);
                      break;
            case 'e': e_flag = true;
                      break;
            case 'M': m_max = atof(poptarg);
	              lm_max = log10(m_max);
                      break;
            case 'm': m_min = atof(poptarg);
	              lm_min = log10(m_min);
                      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 'N': n_init = atoi(poptarg);
		      break;
            case 'Q': q_max = atof(poptarg);
                      break;
            case 'q': q_min = atof(poptarg);
                      break;
            case 'R': R_flag = false;
                      break;
            case 'T': end_time = atof(poptarg);
                      break;
            case 's': input_seed = atoi(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (n <= 0) err_exit("mknodes: N > 0 required!");

    actual_seed = srandinter(input_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

    // Create flat tree 
    node *root  = mknode(n);
    root->log_history(argc, argv);
    root->log_comment(seedlog);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    node* the_binary = root->get_oldest_daughter();

    if (!R_flag) {

      m_prim = m_max;
      m_sec  = m_min;
      sma    = a_min;
      ecc    = 0;
      
      the_binary->set_mass(m_prim);
      
      //		Add the stars and secondary.
      add_secondary(the_binary, m_sec);
      addstar(root, start_time, type);
      double_star* ds =	new_double_star(the_binary, sma, ecc, 
					start_time, n_init, bin_type);
      ds->set_use_hdyn(false);
      ds->get_primary()->set_identity(0);
      ds->get_secondary()->set_identity(1);
	
      node *b      = root->get_oldest_daughter();
      starbase *s  = b->get_starbase();
      star *st     = dynamic_cast(star*, b->get_starbase());

      evolve_binary(the_binary, 0, end_time);

    }
    else
      for (int i=0; i<n; i++) {

      // mkrandom_binary(lm_min, lm_max, la_min, la_max,
      //		 m_prim, m_sec, sma, ecc);

      mkrandom_binary(lm_min, lm_max,
		      2.35,	      // (Added by Steve to allow compilation)
		      la_min, la_max,
		      false,	      // (same here...)
		      m_prim, m_sec, sma, ecc);

      the_binary->set_mass(m_prim);

      //		Add the stars and secondary.
      add_secondary(the_binary, m_sec);
      addstar(dynamic_cast(node*, root), start_time, type);
      double_star* ds =	new_double_star(the_binary, sma, ecc, 
					start_time, i + n_init, bin_type);

      ds->set_use_hdyn(false);
      ds->get_primary()->set_identity(0);
      ds->get_secondary()->set_identity(1);
	
      node *b      = root->get_oldest_daughter();
      starbase *s  = b->get_starbase();
      star *st     = dynamic_cast(star*, b->get_starbase());
   
      evolve_binary(the_binary, 0, end_time);

      delete ds->get_primary();
      delete ds->get_secondary();
      delete ds;

    }
}

#endif
