
//// readstoa:  convert NEMO "stoa" format (N, ndim, time,
////                                        mass[i], i = 1,...,N,
////                                        pos[i],  i = 1,...,N,
////                                        vel[i],  i = 1,...,N)
////            data into a Starlab snapshot.
////
//// Options:     -i    number the particles sequentially [don't number]

//	      Jun Makino, Aug 1996

#include "dyn.h"

#ifdef TOOLBOX

void main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "i";

    bool i_flag = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'i': i_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}


    dyn *root, *by, *bo;

    // Create root node.

    root = new dyn();
    if (i_flag) root->set_label(0);

    int n; scanf("%d", &n); PRL(n);
    int ndim; scanf("%d", &ndim); PRL(ndim);
    real time; scanf("%lf", &time); PRL(time);
    root->set_system_time(time);
    
    // Create first daughter node.

    bo = new dyn();
    root->set_oldest_daughter(bo);
    bo->set_parent(root);
    if (i_flag) bo->set_label(1);

    // Create other daughter nodes.

    for (int i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag) by->set_label(i+1);
	by->set_parent(root);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
	by->set_mass(bo->get_mass());
        bo = by;
    }

    real total_mass = 0;

    for_all_daughters(dyn, root, b) {
	real mass; cin >> mass; 
	b->set_mass(mass);
	total_mass += mass;
    }

    root->set_mass(total_mass);

    for_all_daughters(dyn, root, b) {
	vector pos; cin >> pos; 
	b->set_pos(pos);
    }

    for_all_daughters(dyn, root, b){
	vector vel; cin >> vel; 
	b->set_vel(vel);
    }

    root->log_history(argc, argv);
    put_node(cout, *root);
}

#endif

/* end of: readstoa.c */
