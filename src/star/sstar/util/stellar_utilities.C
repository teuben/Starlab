//// stellar_utilities: determines radius of star with mass M
////                    at specified time t and time at which moment
////                    star first reaches radius R.
////
//// Options:    -m    Mass of the star [in solar units]
////             -r    Radius of the star [in solar units]
////             -t    Time to which star is evolved [in Million year].
////
//++ Notes:
//++  As in starev no kicks are applied in the supernova event.
//++
//++ Example of usage:      
//++  stellar_utilities -m 10 -R 500
//++
//++ See also: addstar
//++           mkmass
//++           mknode
//++           starev


//// Latest version (SPZ:1.0) February 1993.
//
// stellar_radius.C
//

#include "node.h"
#include "star_state.h"
#include "main_sequence.h"

#ifdef TOOLBOX

node* get_complete_star(const real mass, const real time) {

      int number = 1;
      node *n = mknode(number);

      stellar_type type = Main_Sequence;
      int id = 0;
      real start_time = 0;
      real relative_mass, total_mass;
      relative_mass = total_mass = mass;
      real core_mass=0.01;
      real co_core=0;
      real p_rot=0, b_fld=0; 

      node* bi = n->get_oldest_daughter();
        single_star* element = 
		new_single_star(type, id, time, start_time, relative_mass, 
				total_mass, core_mass, co_core, 
				p_rot, b_fld, bi);

      bi->get_starbase()->evolve_element(time);
      return bi;
   }

real get_stellar_radius(const real mass, const real time, stellar_type& type) {

      int number = 1;
      node *n = mknode(number);

      type = Main_Sequence;
      int id = 0;
      real start_time = 0;
      real relative_mass, total_mass;
      relative_mass = total_mass = mass;
      real core_mass=0.01;
      real co_core = 0;
      real p_rot=0, b_fld=0; 

      node* bi = n->get_oldest_daughter();
        single_star* element =
                new_single_star(type, id, time, start_time, relative_mass,
				total_mass, core_mass, co_core,
				p_rot, b_fld, bi);

      bi->get_starbase()->evolve_element(time);
      type = bi->get_starbase()->get_element_type();
      return bi->get_starbase()->get_effective_radius();
   }

real obtain_maximum_radius(const real mass, 
		           real& time, stellar_type& type) {

      type = Main_Sequence;
      time = 0;

      int number = 1;
      node *n = mknode(number);

      int id = 0;
      real relative_mass, total_mass;
      relative_mass = total_mass = mass;
      real core_mass=min(0.01, total_mass);
      real co_core = 0;
      real p_rot=0, b_fld=0;

      node* bi = n->get_oldest_daughter();
        single_star* element =
                new_single_star(type, id, time, time, relative_mass,
				total_mass, core_mass, co_core
				p_rot, b_fld, bi);

      real radius=0;
      real r_max = 0;
      real dt = 0;
      real time_max=0;
      stellar_type current_type;
      do {
         dt = max(cnsts.safety(minimum_timestep),
                  bi->get_starbase()->get_evolve_timestep());
         time += dt;
         bi->get_starbase()->evolve_element(time);
         radius = bi->get_starbase()->get_effective_radius();
         current_type = bi->get_starbase()->get_element_type();
         
         if (radius>=r_max) {
            r_max = radius;
            time_max = time;
            type = current_type;
         }
      } 
      while (current_type!=Black_Hole   &&
             current_type!=Neutron_Star &&
             current_type!=White_Dwarf  &&
             current_type!=Brown_Dwarf  &&
             current_type!=Disintegrated);
         
      time = time_max;
      return r_max;
   }

bool close_enough(const real a, const real b) {

      if (abs((a-b)/b)<=0.0001)
         return false;
      return true;
   }

real obtain_time_at_radius(const real mass,
                           real& radius, stellar_type& type) {

      type = Main_Sequence;
      real time = 0;

      stellar_type zams_type = Main_Sequence;
      real zero = 0;
      real r_zams = get_stellar_radius(mass, zero, zams_type);
      if (radius<=r_zams)
         return zero;

      real time_max;
      stellar_type type_max;
      real r_max = obtain_maximum_radius(mass, time_max, type_max);
      if (radius>=r_max) {
         type = type_max;
         return time_max;
      }

      real time_min = time;
      real dt;
      real previous_radius, current_radius=0;
      do {
         dt = 0.5*(time_max-time_min);
         time += dt;
         previous_radius = current_radius;
         current_radius = get_stellar_radius(mass, time, type);
         if (current_radius>radius) {
            time_max = time;
            time -= dt;
         }
         else
            time_min = time;
      }
      while (close_enough(radius, current_radius) && 
             close_enough(previous_radius, current_radius));
      radius = current_radius;

      return time;
   }

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv) {

    int  c;
    real mass = 10;
    real time = 10;
    real radius = 250;

    char  *comment;
    extern char *poptarg;
    int  pgetopt(int, char **, char *);

    check_help();

    char* param_string = "m:r:t:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
        switch(c) {
            case 'm': mass = atof(poptarg);
                      break;
            case 'r': radius = atof(poptarg);
                      break;
            case 't': time = atof(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
            }

      node* str = get_complete_star(mass, time);
      stellar_type type = NAS;
      real rad = get_stellar_radius(mass, time, type);

      cout << rad 
           <<" = get_stellar_radius(mass=" <<mass<<", time=" 
           <<time<<", type="<<type<<")"<<endl;

      rad = obtain_maximum_radius(mass, time, type);
      cout << rad << " = obtain_maximum_radius(mass=" <<mass<<", time=" 
                  <<time<<", type="<<type<<")"<<endl;

      time = obtain_time_at_radius(mass, radius, type);
      cout << time  << " = obtain_time_at_radius(mass=" <<mass<<", radius="  
                    <<radius<<", type="<<type<<")"<<endl;
    }

#endif
