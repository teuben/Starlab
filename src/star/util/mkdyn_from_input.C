
//// mkdyn_from_input:  create a linked list of dyns and read pos etc. 
////                    from input
////
//// Options:    -F      Specify input filename with the following format:
////                     number of stars
////                     id x m v t j
////          where:     id: identity
////                     x:  position vector
////                     m:  mass
////                     v:  velocity vector
////                     t:  time
////                     j:  angular momentum vector
////
//	      Simon Portegies Zwart, April 1999

#include "dyn.h"

#ifndef TOOLBOX
#else

local void mk_mrv(dyn* b, bool i_flag, int id, 
		          bool m_flag, real mass, 
                          bool r_flag, vector pos, 
                          bool v_flag, vector vel) { 

    b->set_index(id);
    b->set_mass(mass);
    b->set_pos(pos);
    b->set_vel(vel);

}

local dyn * initialize_dyn(istream &in, int n, bool i_flag, 
			   bool m_flag, bool r_flag,
			   bool v_flag, bool j_flag) {


  if(n<=0)
    in >> n;

  vector r, v, j;
  real m, time;
  int id;
  dyn *root = mkdyn(n);
  for_all_daughters(dyn, root, b) {
    if (i_flag && m_flag && r_flag)
      in >> id >> m >> r;
    else
      in >> id >> r >> m >> v >> time >> j;
//    PRC(m);PRC(x);

    mk_mrv(b, i_flag, id, 
	      m_flag, m, 
	      r_flag, r, 
	      v_flag, v); 
  }

  return root;
}

void main(int argc, char ** argv)
{
    bool F_flag = false;
    bool i_flag = false;
    bool m_flag = false;
    bool r_flag = false;
    bool v_flag = false;
    bool j_flag = false;
    char *filename;

    int n = -1;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "F:imrvjn:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'F': F_flag = true;
		      filename = poptarg;
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'm': m_flag = true;
		      break;
	    case 'r': r_flag = true;
		      break;
	    case 'v': v_flag = true;
		      break;
	    case 'j': j_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    dyn * root = NULL;
    if (F_flag) {
      ifstream infile(filename);
      if (!infile) cerr << "error: couldn't create file "
	                  << filename <<endl;
      cerr << "Reading input from file "<< filename <<endl;

      root = initialize_dyn(infile, n, i_flag, m_flag, r_flag, v_flag, j_flag);
    }
    else {
      root = initialize_dyn(cin, n, i_flag, m_flag, r_flag, v_flag, j_flag);
    }

    root->log_history(argc, argv);
    put_dyn(cout, *root);
    rmtree(root);
}

#endif

/* end of: mknode.c */

