
//// mkblack_hole: Replaces the star closest to com with a black hole
////                
////
//// Options:      -M    select black hole mass
////                     if <1 mass is read as a fraction

//                 Simon Portegies Zwart, MIT November 2000

#include "dyn.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

local void mkblack_hole(dyn* b, real m_bh) {

    b->to_com();

    real r_min = VERY_LARGE_NUMBER;
    dyn *bh = NULL;
    for_all_daughters(dyn, b, bi) {
	if(abs(bi->get_pos()) < r_min) {
	    r_min = abs(bi->get_pos());
	    bh = bi;
	}
    }

    PRL(m_bh);
    if(m_bh<=1) {
	cerr << "fractional bh mass" << endl;
	m_bh *= b->get_mass()-bh->get_mass();
	PRL(m_bh);
    }

    if(bh) 
	bh->set_mass(m_bh);
    else 
	err_exit("mkblack_hole: central star not found");

    real m_sum = 0;
    for_all_leaves(dyn, b, bi) {
	m_sum += bi->get_mass();
    }

    b->set_mass(m_sum);


    putiq(bh->get_log_story(), "black_hole", 1);
}

void main(int argc, char ** argv) {

    real m_bh = 0;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "M:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'M': m_bh = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}


    dyn* b;
    b = get_dyn(cin);
    b->log_history(argc, argv);

    mkblack_hole(b, m_bh);

    put_dyn(cout, *b);
    rmtree(b);

}
#endif
