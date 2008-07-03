//// nstarev: evolve a cluster of single stars.
////          the single stars should be provided in the input stream.
////
//// Options:    -c    comment to put in the starbase log structure
////             -n    number of stars to evolve to t_end [all]
////             -s    Random seed
////             -S    evolve until TAMS and not further [false]
////             -t    start time of the stellar evolution in Million year [0]
////             -T    end time of the stellar evolution in Million year [-]
////             -v    verbose output [false]
////
//++ Notes:
//++  As in starev no kicks are applied in the supernova event.
//++
//++ Example of usage:      
//++  mknode -n 10 | mkmass -u 100 -l 10 | addstar -R 1 | nstarev 
//++
//++ See also: addstar
//++           mkmass
//++           mknode
//++           starev


//// Latest version (SPZ:1.0) February 1993.

//   version 1:  Februari 1993   Simon F. Portegies Zwart
//                               spz@grape.c.u-tokyo.ac.jp
//   version 2:  May 2006        Simon F. Portegies Zwart

#include "dyn.h"
#include "single_star.h"
#include "main_sequence.h"

#define EPSILON 1.e-10

#ifdef TOOLBOX

#define DEBUG false

real get_effective_radius(dyn *b) {
     return b->get_starbase()->conv_r_star_to_dyn( 
            b->get_starbase()->get_effective_radius());
}

real get_total_mass(dyn *b) {
      return b->get_starbase()->conv_m_star_to_dyn(
             b->get_starbase()->get_total_mass());
}

stellar_type get_element_type(dyn *b) {
      return (stellar_type)b->get_starbase()->get_element_type();
}

local void evolve_the_stars(node* bi, const real end_time) {

  if(DEBUG) 
      cerr << "evolve_the_stars: " <<bi->format_label() <<" "<<end_time;
  
    real current_time = ((star*)bi->get_starbase())->get_current_time();
    real time_step    =  bi->get_starbase()->get_evolve_timestep();
    
    // PRC(end_time); PRC(current_time); PRC(time_step);

    // SPZ 20 Aug 2004: Just trying this, as it did not seem to work properly...

    while (end_time>current_time+time_step) {

      //  while (end_time>current_time) {
      //      cerr << "running loop at " << time_step << " for star "
      //           << bi->get_index() << endl;

       bi->get_starbase()->evolve_element(current_time+time_step);

       bi->get_starbase()->evolve_element(
             Starlab::min(current_time+time_step+EPSILON, end_time));
       current_time = ((star*)bi->get_starbase())->get_current_time();
       time_step    =  bi->get_starbase()->get_evolve_timestep();

      if(DEBUG) 
	  cerr <<" -dt="<<current_time<<" "<<time_step<<endl;
  }

  bi->get_starbase()->evolve_element(end_time);

  // bi->get_starbase()->get_seba_counters()->step_sstar++;

    if(DEBUG) 
	cerr<<" and leave"<<endl;
}

bool stellar_evolution(dyn *b, 
		       int n_outp,
		       real fraction_of_evolved_stars, 
		       bool TAMS = false,
		       bool verbose = false)
{

  PRL(fraction_of_evolved_stars);
    bool update_all_masses = false;
    
    real stellar_evolution_time = b->get_starbase()->
		conv_t_dyn_to_star(b->get_system_time());

    // Special treatment of black-hole "radii" implemented by Steve (1/05).
    // For efficiency in the integrator, it is useful to assign a negative
    // radius to each black hole.  The value returned by get_radius() will
    // always be positive, but member functions that "need to know" can
    // test the value directly and take appropriate action.

    cerr << "in stellar_evolution() at time " << b->get_system_time()
	 << endl << flush;

    real sum = 0;
    for_all_leaves(dyn, b, bi) {
	if (! bi->is_low_level_node()
	    || !((star*)(bi->get_starbase()))->is_binary_component()) {

	  if(TAMS)
	    stellar_evolution_time = ((single_star*)(bi->get_starbase()))
                                                       ->get_next_update_age();

	  real dt = stellar_evolution_time/(real)(n_outp+1);
	  sum += fraction_of_evolved_stars;
	  if (sum >= 1) {
	    sum -= 1;		

	    //	    cerr << "evolve_star: i= " << bi->get_starbase()->get_identity() <endl;

	    real time = 0;
	    do {
	      time += dt;

	    starbase *sb = bi->get_starbase();
	    if (DEBUG) ((star*)sb)->dump(cerr);


	    real old_dyn_mass = bi->get_mass();
	    evolve_the_stars(bi, time);

	    sb = bi->get_starbase();		// may have changed
	    
	    real new_dyn_mass_from_star = get_total_mass(bi);
	    real new_radius = sb->conv_r_star_to_dyn(sb->get_effective_radius());

	    if (old_dyn_mass <= 0 || new_dyn_mass_from_star <= 0) {
	      PRC(old_dyn_mass); PRC(new_dyn_mass_from_star); PRL(new_radius);
	      ((star*)sb)->dump(cerr);
	    }

	    // Special treatment of black hole dynamical radius.

	    if (sb->get_element_type() == Black_Hole) new_radius *= -1;
	    bi->set_radius(new_radius);

	    if (DEBUG) ((star*)sb)->dump(cerr);
 
	    if (abs(new_dyn_mass_from_star-old_dyn_mass)/old_dyn_mass
		>=cnsts.star_to_dyn(stellar_mass_update_limit)) {
	      update_all_masses = true;
	    }

	    if(verbose)
	      put_node(b);

	    }
	    while(time<stellar_evolution_time);
	  }
	  else {
	    // do not evolve star but set its time to t_end
	    ((star*)bi->get_starbase())->set_current_time(stellar_evolution_time);
	  }
	}
    }

    return update_all_masses;
}

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

#undef EPSILON

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

int main(int argc, char ** argv)
{
    bool s_flag = false;
    bool c_flag = false;
    bool T_flag = false;
    bool verbose = false;
    bool TAMS = false;

    int n_outp = 0;   // number of output steps [0]

    int n_evolve = -1;   // number of output steps [all]
    real fraction_of_evolved_stars = 1;
    real t_end   = 100;  // end time of evolution in Myr [100]

    int  input_seed, actual_seed;
    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];
    extern char *poptarg;
    int c;
    const char *param_string = "N:n:f:S:s:T:t:v:c:";

    check_help();

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'n': n_evolve = atoi(poptarg);
		      break;
	    case 'N': n_outp = atoi(poptarg);
		      break;
	    case 'f': fraction_of_evolved_stars = atof(poptarg);
		      break;
	    case 't': t_end = atof(poptarg);
                      break;
   	    case 'T': TAMS = atof(poptarg);
                      break;
   	    case 'S': TAMS = !TAMS;
                      break;
            case 's': s_flag = TRUE;
                      input_seed = atoi(poptarg);
                      break;
	    case 'v': verbose = !verbose;
                      break;
            case 'c': c_flag = TRUE;
                      comment = poptarg;
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}

    if(!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    dyn* b = get_dyn();

    if(n_evolve>0) 
      fraction_of_evolved_stars = (real)n_evolve/b->n_leaves();

    if(!find_qmatch(b->get_oldest_daughter()
		    ->get_starbase()->get_star_story(), "Type"))
      exit(-1);
	
    b->log_history(argc, argv);

    addstar(b,                            // Note that T_start and
	    b->get_system_time(),          // Main_Sequence are
	    Main_Sequence,                 // defaults. They are
	    true);                         // ignored if a star
    b->set_use_sstar(true);                
    b->get_starbase()->set_use_hdyn(false);

    b->set_system_time(b->get_starbase()->conv_t_star_to_dyn(t_end));
    stellar_evolution(b, n_outp, fraction_of_evolved_stars, TAMS, verbose);

    if(verbose) {
      for_all_daughters(dyn, b, bi) {
	cerr << "Time = " << bi->get_starbase()->get_current_time()
	     << " [Myear],  mass = " << bi->get_starbase()->get_total_mass() 
	     << " [Msun],  radius = " << bi->get_starbase()
	                                   ->get_effective_radius() 
	     << "   " << type_string(bi->get_starbase()->get_element_type())
	     << endl;
      }
    }

    put_node(b);
    rmtree(b);
}

#endif
