//// nstarev: evolve a cluster of single stars.
////          the single stars should be provided in the input stream.
////
//// Options:    -c    comment to put in the starbase log structure.
////             -M    Mass scaling [total cluster mass in solar units].
////             -n    number of output timesteps (timesteps are taken
////                   with constant time intervals) 
////             -R    Dynamical size scaling for the star
////                   [in units of the virial radius].
////             -S    Random seed.
////             -s    Initial stellar type [default is main_sequence].
////             -T    Dynamical time scaling
////                   [in units of the NBODY time].
////             -t    end time of the stellar evolution [in Million year].
////
//++ Notes:
//++  As in starev no kicks are applied in the supernova event.
//++
//++ Example of usage:      
//++  mknode -n 10 | mkmass -u 100 -l 10 | addstar -T 1 | nstarev
//++
//++ See also: addstar
//++           mkmass
//++           mknode
//++           starev


//// Latest version (SPZ:1.0) February 1993.

//   version 1:  Februari 1993   Simon F. Portegies Zwart
//                               spz@grape.c.u-tokyo.ac.jp

//++ Due to overhull currently not operational.

//#include "node.h"
#include "dyn.h"
#include "single_star.h"
#include "main_sequence.h"
//#include "sstar_to_dyn.h"

#define EPSILON 1.e-10

#ifdef TOOLBOX

#define  SEED_STRING_LENGTH  60



local void alter_star_from_random_time(dyn* b) {

    star* ss = ((star*)b->get_starbase());
    b->set_mass(b->get_starbase()
		 ->conv_m_star_to_dyn(b->get_starbase()
				        ->get_total_mass()));

    ss->set_relative_age(b->get_starbase()
			   ->get_current_time());
    ss->set_current_time(b->get_starbase()
			   ->conv_t_dyn_to_star(b->get_system_time()));
    ss->set_current_time(0);
}

local void evolve_star_until_next_time(dyn* bi, const real out_time) {

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

//             print_star(bi->get_starbase(), cerr);
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

local void  evolve_the_stellar_system(dyn* b, real time) {

      dyn * bi;

      if (b->get_oldest_daughter() != NULL)
         for (bi=b->get_oldest_daughter(); bi != NULL;
                                           bi=bi->get_younger_sister())
         evolve_the_stellar_system(bi, time);
      else {
	evolve_star_until_next_time(bi, time);
//         b->get_starbase()->evolve_element(time);
      }
   }

int main(int argc, char ** argv) {

    bool s_flag = false;
    bool c_flag = false;
    bool m_flag = false;
    int  n = 1;
    int n_steps = 1;

    int  input_seed, actual_seed;
    int id=1;
    real m_tot=1,r_hm=1, t_hc=1;
    real m_rel=5, m_core=0.01;
    stellar_type type = Main_Sequence;
    real t_start = 0;
    real t_end   = 100;
    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];
    extern char *poptarg;
    int c;
    char* param_string = "M:N:n:r:s:t:T:";

    check_help();

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

            case 'r': r_hm = atof(poptarg);
                      break;
            case 'M': m_flag = true;
		      m_tot = atof(poptarg);
                      break;
	    case 'N': n = atoi(poptarg);
		      break;
	    case 'n': n_steps = atoi(poptarg);
		      break;
            case 'T': t_hc = atof(poptarg);
                      break;
            case 't': t_end = atof(poptarg);
                      break;
            case 'S': type = (stellar_type)atoi(poptarg);
                      break;
            case 's': s_flag = TRUE;
                      input_seed = atoi(poptarg);
                      break;
            case 'c': c_flag = TRUE;
                      comment = poptarg;
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}

    cerr.precision(HIGH_PRECISION);

    if(!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    // make flat tree 
    dyn* root = get_dyn(cin);
    //    node *root  = mknode(1);
    root->log_history(argc, argv);
//    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

//    node *the_star = root->get_oldest_daughter();

    addstar(root, t_start, type);

    cerr.precision(STD_PRECISION);

    //    put_node(cout, *root);

    real dt, time = 0;
    real delta_t = t_end/((real)n_steps);
    real out_time, current_time; 
    real time1, time2;
    real previous_time;
    int  nstps=0;
    if(t_end>0)
    for_all_daughters(dyn, root, bi) {
       out_time = 0;
       do {
          out_time = Starlab::min(out_time+delta_t, t_end);
          evolve_star_until_next_time(bi, out_time);
//	  bi->get_starbase()->dump(cerr, false);
       }
       while(out_time < t_end);
    }
    else if(t_end == 0) {
      real m, t_ZAMS;
      for_all_daughters(dyn, root, bi) {
	m = bi->get_starbase()->get_total_mass();
        t_ZAMS =  main_sequence_time(m);
        evolve_star_until_next_time(bi, randinter(0, t_ZAMS));
        bi->set_mass(bi->get_starbase()
                       ->conv_m_star_to_dyn(bi->get_starbase()
                                              ->get_total_mass()));
      }
    }
    else
      for_all_daughters(dyn, root, bi) {
        evolve_star_until_next_time(bi, randinter(0, -t_end));

	alter_star_from_random_time(bi);

      }
#if 0
    for_all_daughters(dyn, root, bi) {
       cerr << "Time = " << bi->get_starbase()->get_current_time()
	    << " [Myear],  mass = " << bi->get_starbase()->get_total_mass() 
	    << " [Msun],  radius = " << bi->get_starbase()
                                          ->get_effective_radius() 
            << "   " << type_string(bi->get_starbase()->get_element_type())
	    << endl;
    }
#endif

    put_node(cout, *root);
    delete root;

  }

#if 0

////&&&&&&&&&&&&&&&&&&&&&&
    if (s_flag == FALSE)
        input_seed = 0;                         /* default */
    actual_seed = srandinter(input_seed);

    node *b;

    b = get_node(cin);
    if (c_flag == TRUE)
        b->log_comment(comment);
    b->log_history(argc, argv);
    b->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    addstar(b, t_start, type);

    //put_node(cout, *b);

    real t_out, dt = t_end/n;
    for (int i=0; i<n; i++) {
      t_out = dt * (i+1);

      evolve_the_stellar_system(b, t_out);
      //      for_all_daughters(node, b, bi) 
      //	bi->get_starbase()->evolve_element(t_out);

      put_node(cout, *b);
    }
    return 0;
}

#endif
#endif
