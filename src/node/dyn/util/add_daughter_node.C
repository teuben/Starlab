
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add one extra node to the tree structure in the input snapshot,
//// under a specified node.
////
//// Usage: add_daughter_node [OPTIONS] < input > output
////
//// Options:     
////              -c     add a comment to snapshot [false]
////              -e     echo tree structure [false]
////              -i     specify index of node to add to [root]
////              -j     specify index for new node [none]
////              -m     specify node mass [1]
////              -s     specify random seed [take from system clock]
////              -r     specify node radial position (angle random) [0]
////              -R     specify node 3-d position [0,0,0]
////              -v     specify node speed (direction random) [0]
////              -V     specify node 3-d velocity [0,0,0]
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//   version 1:  Dec 1994   (node version) Piet Hut
//   version 2:  Jun 2002   (dyn version)  Steve McMillan
//   version 3:  NOC 2012   (dyn version)  Steve McMillan

#include "dyn.h"

//===========================================================================

#ifdef TOOLBOX

//-----------------------------------------------------------------------------
//  main  --  driver to directly add one extra daughter node.
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    bool c_flag = false;
    bool e_flag = false;     // echo flag: if true, echo tree structure
    bool i_flag = false;
    bool j_flag = false;
    bool s_flag = false;
    bool R_flag = false;
    bool V_flag = false;

    int  i, j;
    int input_seed = 0;
    real m = 1;               // default mass: unity
    real r = 0;
    real v = 0;
    real xx, xy, xz, vx, vy, vz;
    char  *comment;

    check_help();

    extern char *poptarg, *poparr[];
    int  c;
    const char *param_string = "c:ei:j:m:r:R:::s:v:V:::";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
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
	    case 'R': xx = atof(poparr[0]);
		      xy = atof(poparr[1]);
		      xz = atof(poparr[2]);
		      R_flag = true;
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'v': v = atof(poptarg);
		      break;
	    case 'V': vx = atof(poparr[0]);
		      vy = atof(poparr[1]);
		      vz = atof(poparr[2]);
		      V_flag = true;
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

    vec pos = 0;
    if (R_flag)
	pos = vec(xx,xy,xz);
    else {
	real costheta = randinter(-1, 1);
	real sintheta = sqrt(1-costheta*costheta);
	if (randinter(-1,1) > 0) sintheta = -sintheta;
	real phi = randinter(0, 2*M_PI);
	pos = r*vec(sintheta*cos(phi), sintheta*sin(phi), costheta);
    }
    PRL(pos);

    vec vel = 0;
    if (V_flag)
	vel = vec(vx,vy,vz);
    else {
	real costheta = randinter(-1, 1);
	real sintheta = sqrt(1-costheta*costheta);
	if (randinter(-1,1) > 0) sintheta = -sintheta;
	real phi = randinter(0, 2*M_PI);
	vel = v*vec(sintheta*cos(phi), sintheta*sin(phi), costheta);
    }

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

    // Believe and correct initial_mass if system_time = 0.

    if ((real)root->get_system_time() == 0
	&& find_qmatch(root->get_log_story(), "initial_mass")) {
	real mass = getrq(root->get_log_story(), "initial_mass") + m;
	putrq(root->get_log_story(), "initial_mass", mass, HIGH_PRECISION);
    }

    put_dyn(root);
    rmtree(root);
}

#endif

//===========================================================================

// endof: add_daughter_node.C

