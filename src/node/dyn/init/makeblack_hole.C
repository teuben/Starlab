
//// mkblack_hole: Replaces the star closest to com with a black hole
////                
////
//// Options:      -M    select black hole mass
////                     if <1 mass is read as a fraction
//// Options:      -i    select black hole identity

//                 Simon Portegies Zwart, MIT November 2000

#include "dyn.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

local void mkblack_hole(dyn* b, int id, real m_bh) {

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
        bh->set_vel(bh->get_vel() * sqrt(prev_mass/m_bh));
    }
    else 
	err_exit("mkblack_hole: selected star not found");

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

int main(int argc, char ** argv) {

    real m_bh = 0;
    int id = -1;    // most central particle is selected

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "i:M:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'M': m_bh = atof(poptarg);
		      break;
	    case 'i': id = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    dyn* b;
    b = get_dyn();
    b->log_history(argc, argv);

    if (id>b->n_leaves())
      err_exit("selected id exceeds particle number");

    mkblack_hole(b, id, m_bh);

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
    return 0;
}
#endif
