//// starev: evolve a single star.
////         creates a single star and evolves it in time.
////
//// Options:    -c    comment to put in the starbase log structure.
////             -M    initial mass of the star [in solar units].
////             -n    number of output timesteps (timesteps are taken
////                   with constant time intervals) 
////             -R    Dynamical size scaling for the star
////                   [in units of the virial radius].
////             -S    Random seed.
////             -s    Initial stellar type [default is main_sequence].
////             -T or -t end time of the stellar evolution [in Million year].
////
//// Latest version (SPZ:1.0) February 1993.

//++ the single stars are part of the dynamical tree.
//++ Since the assymetry in supernovae are taken care of by the 
//++   dynamical model no kicks are applied.

//   version 1.0:  Februari 1993   Simon F. Portegies Zwart
//                                 spz@grape.c.u-tokyo.ac.jp 

//#include "dyn.h"
#include "node.h"
#include "single_star.h"
#include "main_sequence.h"
//#include "sstar_to_dyn.h"
#define EPSILON 1.e-10

#ifdef TOOLBOX

local void evolve_star_until_next_time(node* bi, const real out_time) {

          real current_time = ((star*)bi->get_starbase())->get_current_time();
          real time_step    =  bi->get_starbase()->get_evolve_timestep();

          while (out_time>current_time+time_step) {
             bi->get_starbase()->evolve_element(current_time+time_step);
             bi->get_starbase()->evolve_element(
                 Starlab::min(current_time+time_step+EPSILON, out_time));
             current_time = ((star*)bi->get_starbase())->get_current_time();
             time_step    =  bi->get_starbase()->get_evolve_timestep();

	     star_state ss(dynamic_cast(star*, bi->get_starbase()));
	     put_state(ss, cerr);

	     // print_star(bi->get_starbase(), cerr);

	     int p = cerr.precision(HIGH_PRECISION);
	     bi->get_starbase()->dump(cerr, false);
	     cerr.precision(p);
          }
         bi->get_starbase()->evolve_element(out_time);

         print_star(bi->get_starbase(), cerr);

/*
cerr<< "ov: " << bi->get_starbase()->get_element_type() 
        << " " <<bi->get_starbase()->get_total_mass() << " t= "
        << ((star*)bi->get_starbase())->get_current_time() << " "
        << ((star*)bi->get_starbase())->get_current_time()
           +bi->get_starbase()->get_evolve_timestep() << " -> "
        << out_time << endl;
*/
}

/*-----------------------------------------------------------------------------
 *  main  --
 *	usage:
 *		addstar -t # [options]  ,
 *
 *		where # is the initial age of the cluster.
 *	options:
 *            	The following options are allowed:
 *	cluster age:
 *		-t #	Where # stands for the initial age of the
 *			in Myear.
 *
 *		At present the running time of the integrator correspnds
 *		to the stellar age an a one by 10e6year basis.
 *		This however should be scaled to the cluster parameters.
 *-----------------------------------------------------------------------------
 */
int main(int argc, char ** argv)
    {
    stellar_type type = Main_Sequence;
    char * star_type_string;
    int  c;

    bool  t_flag = FALSE;
    bool  S_flag = FALSE;
    bool  c_flag = FALSE;
    bool  M_flag = FALSE;
    bool  n_flag = FALSE;
    bool  R_flag = FALSE;
    real  m_tot = 10;
    real  r_hm = 100;
    real  t_hc = 1;
    real  t_start = 0;           // default value;
    real  t_end = 100;
    int n_steps = 1;
    char  *comment;
    int input_seed=0, actual_seed;
    extern char *poptarg;
    int  pgetopt(int, char **, char *);
    char *param_string = "n:M:R:T:t:S:s:c:";

    check_help();
    
    if (argc <= 1)
        {
        cerr <<"usage: starev -M # -R # -T # -t # -S # -s # [-c \"..\"]\n";
        exit(1);
        }
 
    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
            case 'n': n_steps = atoi(poptarg);
                      break;
            case 'M': M_flag = TRUE;
	              m_tot = atof(poptarg);
                      break;
            case 'R': r_hm = atof(poptarg);
                      break;
            case 'T':                      
            case 't': t_end = atof(poptarg);
                      break;
            case 'S': S_flag = TRUE;
                      input_seed = atoi(poptarg);
                      break;
            case 's': strcpy(star_type_string, poptarg);
                      type = extract_stellar_type_string(star_type_string);
                      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	    }

    cerr.precision(HIGH_PRECISION);

    if(!S_flag) actual_seed = 0;
    actual_seed = srandinter(input_seed);

    // make flat tree 
    node *root  = mknode(1);
    root->log_history(argc, argv);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    node *the_star = root->get_oldest_daughter();

    addstar(root, t_start, type);

    // Starev does not include hdyn.h nor the hdyn library
    // get_use_hdyn is therefore not defined.
    // The result is that kick velocities will be scaled
    // spuriously......
    //    root->set_use_hdyn(false);
    cerr.precision(STD_PRECISION);

    //    put_node(root);

    real dt, time = 0;
    real delta_t = t_end/((real)n_steps);
    real out_time, current_time; 
    real time1, time2;
    real previous_time;
    int  nstps=0;
    for_all_daughters(node, root, bi) {
       out_time = 0;
       do {
          out_time = Starlab::min(out_time+delta_t, t_end);
          evolve_star_until_next_time(bi, out_time);
	  bi->get_starbase()->dump(cerr, false);
       }
       while(out_time < t_end);
    }

    for_all_daughters(node, root, bi) {
       cerr << "Time = " << bi->get_starbase()->get_current_time()
	    << " [Myear],  mass = " << bi->get_starbase()->get_total_mass() 
	    << " [Msun],  radius = " << bi->get_starbase()
                                          ->get_effective_radius() 
            << "   " << type_string(bi->get_starbase()->get_element_type())
	    << endl;
    }
//    put_node(root);
    delete root;
    return 0;
}

#endif
