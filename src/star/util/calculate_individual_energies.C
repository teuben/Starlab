
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print out the energy of an N-body system.  Does not include 
//// the effects of any external tidal field.
////
//// Usage: energy [OPTIONS] < input > output
////
//// Options:   
////		  -e    specify softening parameter [0]
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"

//#ifdef TOOLBOX

enum {unbound_single=0, bound_single, bound_binary, unbound_binary};

// compute_individual_energies: calculate the appropriate color code for particle bi
//		                relative to "particle" bj (usually the root node).

void accumulate_potential_energy(dyn* bj, dyn*bi,
				 real& epot, real& rmin, dynptr& bmin)

// Determine the potential energy of bi with respect to bj, along
// with bi's nearest neighbor and minimum neighbor distance.

{
    if (bj->get_oldest_daughter() != (dyn*)NULL)
	for (dyn * bb = bj->get_oldest_daughter(); bb != (dyn*)NULL;
	    bb = bb->get_younger_sister()) {
	    accumulate_potential_energy(bb, bi, epot, rmin, bmin);
	}
    else
	if (bi != bj) {
	    vec d_pos = bi->get_pos() - bj->get_pos();
	    real mi = bi->get_mass();
	    real mj = bj->get_mass();
	    real r = sqrt(d_pos * d_pos);
	    epot += -mi*mj/r;
	    if (r < rmin) {rmin = r; bmin = bj;}
	}
}

int compute_individual_energies(dyn* bj, dyn* bi)
{
    int c = unbound_single;

    real   epot = 0, rmin = VERY_LARGE_NUMBER;
    dynptr bmin = bi;

    accumulate_potential_energy(bj, bi, epot, rmin, bmin);

    real   mi = bi->get_mass();
    real   mj = bmin->get_mass();
    vec d_vel = bi->get_vel() - bmin->get_vel();
    vec d_pos = bi->get_pos() - bmin->get_pos();
    real   r = sqrt(d_pos * d_pos);
    real   e0 = (0.5 * d_vel * d_vel - (mi + mj)/r);

    if (e0 < 0) {
	real   epot1 = 0, rmin1 = VERY_LARGE_NUMBER;
	dynptr bmin1 = bi;

	accumulate_potential_energy(bj, bmin, epot1, rmin1, bmin1);

	if (bi == bmin1) {
	    real e  = - mi*mj / r;
	    vec R_vel = (mi*bi->get_vel()+mj*bmin->get_vel())/(mi+mj);
	    real ekin = 0.5*(mi+mj)*R_vel*R_vel;

	    if (epot + epot1 - 2*e + ekin < 0) c = bound_binary;
	    else c = unbound_binary;

	} else c = bound_single;

    } else {
	vec vel = bi->get_vel();
	real ekin = 0.5*bi->get_mass()*vel*vel;
	
	if (ekin + epot > 0.0) c = unbound_single;
    }
    return c;
}


main(int argc, char **argv)
{
    dyn *root;
    int  e_flag = 0;
    real eps = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "e:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'e': e_flag = 1;
                      eps = atof(poptarg);
	              break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    // Loop over input until no more data remain.

    while ( (root = get_dyn()) != NULL) {

      int bound;
      for_all_daughters(dyn, root, bi) {
	bound = compute_individual_energies(root, bi);
	PRC(bound);PRC(bi->get_pos());PRC(bi->get_vel());
      }
    }
}

//#endif

// endof: energy.C
