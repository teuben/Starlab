
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Replace the star closest to com by a black hole.
////                
//// Usage:  makeblack_hole [OPTIONS]
////
//// Options:
////      -M    select black hole mass; if <1 mass is read as a fraction
////      -r    select black hole position
////      -v    select black hole velocity
////      -f    flag for seting black hole prameters
////      -i    select black hole identity
////
//// Written by Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//                 Simon Portegies Zwart, MIT November 2000

#include "dyn.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

local void makeblack_hole(dyn* b, int id, real m_bh,
			  bool r_flag, vec r_bh, bool v_flag, vec v_bh) {

    b->to_com();

    real r_min = VERY_LARGE_NUMBER;
    dyn *bh = NULL;
    if(id<0) {
      for_all_daughters(dyn, b, bi) {
	if(abs(bi->get_pos()) < r_min) {
	  r_min = abs(bi->get_pos());
	  bh = bi;
	}
      }
    }
    else {
      for_all_daughters(dyn, b, bi) {
	if(bi->get_index() == id) {
	  bh = bi;
	  break;
	}
      }
    }

    PRL(m_bh);
    if(m_bh<=1) {
	cerr << "fractional bh mass" << endl;
	m_bh *= b->get_mass()-bh->get_mass();
	PRL(m_bh);
    }

    real prev_mass; 
    if(bh) {
        prev_mass = bh->get_mass();
	bh->set_mass(m_bh);
	  cerr << "Set Black hole position and velocity" << endl;
	  if(r_flag) 
	    bh->set_pos(r_bh);
	  if(v_flag) 
	    bh->set_vel(v_bh);
	  else
	    bh->set_vel(bh->get_vel() * sqrt(prev_mass/m_bh));
    }
    else 
	err_exit("makeblack_hole: selected star not found");

    putiq(bh->get_log_story(), "black_hole", 1);

    real m_sum = 0;
    for_all_leaves(dyn, b, bi) {
	m_sum += bi->get_mass();
    }

    b->set_mass(m_sum);

    sprintf(tmp_string,
	    "         black hole added, total mass = %8.2f", m_sum); 
    b->log_comment(tmp_string);
    cerr << "Black hole mass is " << m_bh << endl;
}

local inline real mass(dyn *b, real r) 
{

  real m_sum = 0;
  for_all_leaves(dyn, b, bi) {
    if(abs(bi->get_pos())<=r)
      m_sum += bi->get_mass();
  }
  return m_sum;
}

local inline real vcirc2(dyn *b, real r)	// circular orbit speed:
						// recall vc^2 = r d(phi)/dr
{
    if (r > 0)
	return mass(b, r)/r;
    else
	return 0;
}

int main(int argc, char ** argv) {

    bool  c_flag = FALSE;
    char  *comment;

    real m_bh = 0;
    int id = -1;    // most central particle is selected
    bool r_flag = false;
    vec r_bh = 0;
    bool v_flag = false;
    vec v_bh = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];

    int c;
    char* param_string = "c:i:M:r:::v:::";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
	    case 'M': m_bh = atof(poptarg);
		      break;
	    case 'r':	r_bh = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	                r_flag = true;
	    		break;
	    case 'v':	v_bh = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	                v_flag = true;
	    		break;
	    case 'i': id = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    dyn* b;
    while (b = get_dyn()) {

      if (c_flag == TRUE)
	b->log_comment(comment);
      b->log_history(argc, argv);

      // Check if we have to reinterpret r and v in light of possible
      // external fields and physical parameters.

      //    real mass, length, time;
      //    bool phys = get_physical_scales(b, mass, length, time);

      if (id>b->n_leaves())
	err_exit("selected id exceeds particle number");

      if(v_flag)
	v_bh *= sqrt(vcirc2(b, abs(r_bh)));

      makeblack_hole(b, id, m_bh, r_flag, r_bh, v_flag, v_bh);

      real initial_mass = getrq(b->get_log_story(), "initial_mass");

      if (initial_mass > -VERY_LARGE_NUMBER)
        putrq(b->get_log_story(), "initial_mass", b->get_mass(),
	      HIGH_PRECISION);

      real m_sum = b->get_mass();
      real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
      if(old_mtot!=m_sum) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_sum,
							 old_r_vir,
							 old_t_vir);
      }

      put_dyn(b);
      rmtree(b);
    }
}
#endif
