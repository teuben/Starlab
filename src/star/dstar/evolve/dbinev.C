//
//  mkdouble.C: construct a linked list of unit-mass nodes.
//
//	       Simon Portegies Zwart, July 1996
//

#include "hdyn.h"
#include "single_star.h"
#include "main_sequence.h"
#include "double_star.h"
#include "sstar_to_dyn.h"
#include "dstar_to_dyn.h"


void add_secondary(dyn* original, real mass_ratio) {

    dyn* primary = new dyn;
    dyn* secondary = new dyn;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    primary->set_mass(original->get_mass());
    secondary->set_mass(mass_ratio*original->get_mass());
    original->inc_mass(secondary->get_mass());

    // Naming convention:

    if (original->get_name() == NULL)
        if (original->get_index() >= 0) {
            char tmp[64];
            sprintf(tmp, "%d", original->get_index());
            original->set_name(tmp);
        }

    primary->set_name(original->get_name());
    secondary->set_name(original->get_name());
    strcat(primary->get_name(), "a");
    strcat(secondary->get_name(), "b");

   }

void mksecondary(dyn* b, real binary_fraction, real lower_limit) {

    // For now, use a flat distribution in secondary mass ratio.
    // Assume that the probability of a star being the primary of
    // a binary is independent of mass.

    real sum = 0;
    b->set_mass(0);

    for_all_daughters(dyn, b, bi) {
        sum += binary_fraction;
        if (sum >= 1) {
            sum -= 1;

            real mass_ratio = randinter(lower_limit, 1);        // Quick fix...
            add_secondary(bi, mass_ratio);

        }
        b->inc_mass(bi->get_mass());
    }
}

local void  evolve_binary(dyn* the_binary, real end_time) {

      real time = 0;   
      real dt=0;
      if (!the_binary->is_root())
         if (the_binary->get_parent()->is_root()) {
             cerr<<"binary "<<the_binary<<endl;
             do {
                dt = the_binary->get_starbase()->get_evolve_timestep();
                time += max(0., min(end_time-time, dt));
                cerr<<"dt, t= " << dt<<" " << time<<endl;
                the_binary->get_starbase()->evolve_element(time);
             }
	     while (time<end_time);
         }
   }

local  bool  evolve_the_stellar_system(dyn* b, real time) {

      for_all_nodes(dyn, b, bi) {
         cerr<<bi<< " is_root?: "<<bi->is_root()<<endl;
         cerr<<bi<< " is_parent?: "<<bi->is_parent()<<endl;
         if (!bi->is_root())
            if (bi->get_parent()->is_root()) {
                cerr<<"binary "<<bi<<endl;
                bi->get_starbase()->evolve_element(time);
            }

     }
      return true;
   }

void make_new_kepler_to_dyn(dyn* b, real m_total, real sma, real ecc,
			    int planar=0) {

     kepler *johannes = new kepler;

     real peri = 1; // Default value (unimportant unless ecc = 1).

     // For now, binary phase is random.

     real mean_anomaly = randinter(-PI, PI);

     make_standard_kepler(*johannes, 0, m_total, -0.5 * m_total / sma, ecc,
                          peri, mean_anomaly);
     set_random_orientation(*johannes, planar);

     b->set_kepler(johannes);
   }

void main(int argc, char ** argv) {

    bool m_flag = false;
    bool c_flag = false;

    bool reandom_initialization = false;
    int  n = 1;
    int id=1;
    real m_tot=1,r_hm=1, t_hc=1;
    binary_type type = Detached;
    binary_type bin_type = Detached;
    real binary_fraction = 1.0;
//		Initial conditions from Tutukov & Yunglesson 1993.
    real sma    = 138;
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
    char* comment;
    char* param_string = "a:e:M:m:n:q:s:T:t:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case 'a': sma = atof(poptarg);
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
            case 'T': end_time = atof(poptarg);
	              break;
            case 't': start_time = atof(poptarg);
                      break;
            case 's': random_seed = atoi(poptarg);
                      break;
            case 'c': c_flag = TRUE;
	              comment = poptarg;
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (n <= 0) err_exit("mkdyn: N > 0 required!");
    int actual_seed = srandinter(random_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

    dyn *b;

    b = get_dyn(cin);
    if (c_flag == TRUE)
        b->log_comment(comment);
    b->log_history(argc, argv);
    b->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    adddouble(b, start_time, type);

    put_dyn(cout, *b);

    evolve_binary(b, end_time);

    put_dyn(cout, *b);
}
