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

local void  evolve_binary(node* the_binary, real end_time, int n_step) {

      real time = 0;   
      real dt=0;
      double_star* bin_star = dynamic_cast(double_star*, 
					   the_binary->get_starbase());
	bin_star->dump(cerr);
	cerr<<"time = "<<time<<endl;
	bin_star->put_state();
	cerr << "Porb = " << bin_star->get_period() << endl;

	real dt_min = end_time/(real)n_step;
	
      if (!the_binary->is_root())
         if (the_binary->get_parent()->is_root()) {
             cerr<<"binary "<<the_binary<<endl;
             do {
                dt = min(dt_min,
			 the_binary->get_starbase()->get_evolve_timestep());
                time += max(0., min(end_time-time, dt));
                cerr<<"dt, t= " << dt<<" " << time<<endl;
                the_binary->get_starbase()->evolve_element(time);
		//the_binary->get_starbase()->dump(cerr);
		bin_star->dump(cerr);
		cerr<<"time = "<<time<<endl;
		bin_star->put_state();
		cerr << "Porb = " << bin_star->get_period() << endl;
	    }
	     while (time<end_time);
	     // perform evolution until exact end_time.
	     bin_star->evolve_the_binary(end_time);

         }
   }

void main(int argc, char ** argv) {

    bool m_flag = false;
    bool P_flag = false;
    int  n = 1;
    int id=1;
    stellar_type type = Main_Sequence;
    binary_type bin_type = Detached;
    real binary_fraction = 1.0;
//		Initial conditions from Tutukov & Yunglesson 1993.
    real m_tot = 1;
    real r_hm = 1;
    real t_hc = 1;;
    real sma    = 138;
    real period = 39.3; // days.
    real ecc    = 0.0;
    real m_prim = 13.1;
    real m_sec  =  9.8;
    real mass_ratio = 0.75;
    real lower_limit = 0.0;
    int random_seed = 0;
    char seedlog[64];
    real start_time = 0;
    real end_time   = 35;
    extern char *poptarg;
    int c;
    char* param_string = "a:P:e:M:m:n:q:s:T:R:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case 'a': sma = atof(poptarg);
                      break;
            case 'P': P_flag = true;
	              period = atof(poptarg);
                      break;
            case 'e': ecc = atof(poptarg);
                      break;
            case 'M': m_prim = atof(poptarg);
                      break;
            case 'm': m_flag = true;
		      m_sec = atof(poptarg);
                      break;
	    case 'n': n = atoi(poptarg);
		      break;
            case 'q': mass_ratio = atof(poptarg);
                      break;
            case 'R': r_hm = atof(poptarg);
                      break;
            case 'T': end_time = atof(poptarg);
                      break;
            case 's': random_seed = atoi(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    int actual_seed = srandinter(random_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

    if (m_flag) 
       mass_ratio = m_sec / m_prim;
    else
       m_sec = mass_ratio*m_prim;
    m_tot = m_prim;// + m_sec;

    if (P_flag)
      sma = period_to_semi(period, m_prim, m_sec);
    
    cerr.precision(INT_PRECISION);
    // Create flat tree 
    node *root  = mknode(1);
    root->log_history(argc, argv);
    root->log_comment(seedlog);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    node* the_binary = root->get_oldest_daughter();

//		Add the stars and secondary.
    add_secondary(the_binary, mass_ratio);

    int i = 0;
    for_all_leaves(node, the_binary, bi)
      bi->set_index(i++);
		  
    addstar(root, start_time, type);
    double_star* new_double = 
		 new_double_star(the_binary, sma, ecc, 
				 start_time, id, bin_type);

    root->get_starbase()->set_use_hdyn(false);
    put_node(root);

//	Test pointer structure
    node *b = root->get_oldest_daughter();
    starbase *s = b->get_starbase();
    star *st     = (star*)b->get_starbase();
   
    evolve_binary(the_binary, end_time, n);
    double_star* bin_star = dynamic_cast(double_star*, 
					 the_binary->get_starbase());
    bin_star->dump(cerr);
    put_node(root);
}

#endif
