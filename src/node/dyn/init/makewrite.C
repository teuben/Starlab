
//// makewrite:  Turn a text file into a snapshot.  Used in some
////           demonstrations only.
////
//// Options:  -c    add a comment to the output snapshot [false]
////           -C    output data in 'col' format [no]
////           -r    random velocities [0]
////           -s    random seed
////
//// Example:  /usr/games/banner -w 40 STARLAB | makewrite -c "STARLAB"

#include "dyn.h"
#include <ctype.h>

#ifdef TOOLBOX

#define  LENGTH_SCALE_FACTOR   0.1
#define  LEFT_OFFSET           -20
#define  TOP_OFFSET             -5

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = false;
    bool  C_flag = false;

    real v_rand = 0;

    int input_seed, actual_seed;
    bool s_flag = false;

    check_help();
 
    extern char *poptarg;
    int c;
    char* param_string = "c:Cv:s:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'v': v_rand = atof(poptarg);
	    	      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}            
    
    if (!s_flag) input_seed = 0;
    actual_seed = srandinter(input_seed);

    dyn *b, *by, *bo;
    b = new dyn();
    bo = new dyn();
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    int j = TOP_OFFSET;		// x
    int i = LEFT_OFFSET;	// y
    int indx = 0;
    int n = 0;

    // Banner writes from top to bottom.  Want to translate this
    // into left to right.

    c = getchar();
    while(c != EOF) {
	while (c != '\n' && c != EOF) {
	    if (!isspace(c)) {
 
		n++;
		bo->set_pos(vec(j*LENGTH_SCALE_FACTOR,
				   i*LENGTH_SCALE_FACTOR, 0));
		bo->set_vel(v_rand*vec(randinter(-1,1),
					  randinter(-1,1),
					  randinter(-1,1)));

		bo->set_label(++indx);

 		by = new dyn();
		bo->set_younger_sister(by);
		by->set_elder_sister(bo);
		bo = by;
	    }
	    i++;		// i (= y) increases with each character
	    c = getchar();
	}
	j++;			// j (= x) increases with each new line
	i = LEFT_OFFSET;
	if (c != EOF)
	    c = getchar();
    }

    bo = bo->get_elder_sister();
    delete by;
    bo->set_younger_sister(NULL);

    for_all_daughters(dyn, b, bb)
	bb->set_mass(1.0/n);
    
    if (c_flag == TRUE)
        b->log_comment(comment);
    b->log_history(argc, argv);

    if (C_flag) b->set_col_output(true);

    b->to_com();
    put_dyn(b);
}

#endif

