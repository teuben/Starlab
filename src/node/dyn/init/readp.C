
//// readp:  convert ASCII "dumbp" format (mass1, pos1, vel1,
////                                       mass2, pos2, vel2,
////                                       mass3, pos3, vel3,
////                                       etc.)
////
////         data into a Starlab snapshot (flat tree).
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -i    number the particles sequentially [don't number]

//	     Steve McMillan, July 1999

#include "dyn.h"

#ifdef TOOLBOX

void main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:i";

    char *comment;
    bool c_flag = false;
    bool i_flag = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Create the root node.

    dyn * root = new dyn();
    root->set_system_time(0);
    if (i_flag) root->set_index(0);

    // Create the system node by node.

    dyn *b, *bo;
    real total_mass = 0;
    int n = 0;

    while (!cin.eof()) {

	real mass; cin >> mass; if (cin.eof()) break;	// note: only minimal
	vector pos; cin >> pos; if (cin.eof()) break;	// checking for
	vector vel; cin >> vel;				// corrupted data

	dyn * b = new dyn();

	b->set_mass(mass);
	b->set_pos(pos);
	b->set_vel(vel);

	b->set_parent(root);

	if (n++ == 0)
	    root->set_oldest_daughter(b);
	else {
	    b->set_elder_sister(bo);
	    bo->set_younger_sister(b);
	}

	if (i_flag) b->set_index(n);
	bo = b;
    }
	    
    root->set_mass(total_mass);

    root->log_history(argc, argv);
    if (c_flag) root->log_comment(comment);

    put_node(cout, *root);
}

#endif
