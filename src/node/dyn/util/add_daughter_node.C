
//// add_daughter_node:  add one extra node to the tree structure in the
////                     input snapshot, under a specified node.
////
//// Options:     -c     add a comment to snapshot [false]
////              -e     echo tree structure [false]
////              -i     specify index of node to add to [root]
////              -j     specify index for new node [none]
////              -m     specify node mass [1]
////              -s     specify random seed [take from system clock]
////              -r     specify node radial position (angle random) [0]
////              -v     specify node speed (direction random) [0]

//   version 1:  Dec 1994   (node version) Piet Hut
//   version 2:  Jun 2002   (dyn version)  Steve McMillan

#include "dyn.h"

//===========================================================================

#ifdef TOOLBOX

//-----------------------------------------------------------------------------
//  main  --  driver to directly add one extra daughter node.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool  c_flag = false;
    bool  e_flag = false;     // echo flag: if true, echo tree structure
    bool  i_flag = false;
    bool  j_flag = false;
    bool  s_flag = false;

    int  i, j;
    int input_seed = 0;
    real m = 1;               // default mass: unity
    real r = 0;
    real v = 0;
    char  *comment;

    check_help();

    extern char *poptarg;
    int  c;
    char* param_string = "c:ei:j:m:r:s:v:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
        switch(c) {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'e': e_flag = true;
		      break;
	    case 'i': i_flag = true;
		      i = atoi(poptarg);
		      break;
	    case 'j': j_flag = true;
		      j = atoi(poptarg);
		      break;
	    case 'm': m = atof(poptarg);
		      break;
	    case 'r': r = atof(poptarg);
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'v': v = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            

    if (!s_flag) input_seed = 0;                         	// default
    int actual_seed = srandinter(input_seed);

    char  seedlog[128];
    sprintf(seedlog, "       random number generator seed = %d", actual_seed);
    
    dyn * root;    // root node
    dyn * p;       // parent node
    dyn * d;       // older daughter node
    dyn * y;       // younger daughter node
    dyn * n;       // new daughter node

    root = get_dyn();
    root->log_comment(seedlog);

    if (!root)
	err_exit("add_daughter_node: no input nodes provided");

    if (c_flag) root->log_comment(comment);
    root->log_history(argc, argv);

    if (!i_flag)   // default parent: root
	p = root;
    else
	p = (dyn*)node_with_index(i, root);

    if (p == NULL)
	err_exit("add_daughter_node: no such parent");

    real costheta = randinter(-1, 1);
    real sintheta = sqrt(1-costheta*costheta);
    if (randinter(-1,1) > 0) sintheta = -sintheta;
    real phi = randinter(0, 2*M_PI);
    vec pos = r*vec(sintheta*cos(phi), sintheta*sin(phi), costheta);
    PRL(pos);
    
    costheta = randinter(-1, 1);
    sintheta = sqrt(1-costheta*costheta);
    if (randinter(-1,1) > 0) sintheta = -sintheta;
    phi = randinter(0, 2*M_PI);
    vec vel = v*vec(sintheta*cos(phi), sintheta*sin(phi), costheta);

    n = new dyn();
    n->set_mass(m);
    n->set_pos(pos);
    n->set_vel(vel);
    if (j_flag) n->set_label(j);

    d = p->get_oldest_daughter();
    if (d == NULL)
	p->set_oldest_daughter(n);
    else {
	y = d->get_younger_sister();
	while (y) {
	    d = y;
	    y = d->get_younger_sister();
	}
	d->set_younger_sister(n);
    }

    cerr << "Added daughter node " << n->format_label();
    if (!i_flag)
	cerr << " to root node" << endl;
    else
	cerr << " to node " << p->format_label() << endl;
    
    if (e_flag) root->pretty_print_tree(cerr);

    put_dyn(root);
    rmtree(root);
}

#endif

//===========================================================================

// endof: add_daughter_node.C

