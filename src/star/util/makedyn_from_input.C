
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
////
//	      Simon Portegies Zwart, April 1999

#include "dyn.h"

#ifndef TOOLBOX
#else

local void mk_mrv(dyn* b, int i_flag, int id, 
		          int m_flag, real mass, 
                          int r_flag, vec pos, 
                          int v_flag, vec vel) { 

    b->set_index(id);
    b->set_mass(mass);
    b->set_pos(pos);
    b->set_vel(vel);

}

local dyn * initialize_dyn(istream &in, int n, int i_flag, 
			   int m_flag, int r_flag,
			   int v_flag) {


  if(n<=0)
    in >> n;

  vec r, v;
  real m, time;
  real m_tot = 0;
  int id=0;
  dyn *root = mkdyn(n);
  for_all_daughters(dyn, root, b) {
    if(i_flag < 0 && m_flag < r_flag && r_flag < v_flag) {
      id++;
      in >> m >> r >> v;
    }
    else if (i_flag < m_flag && m_flag < r_flag && r_flag < v_flag)
      in >> id >> m >> r >> v;
    else {
      cerr << "order not not yet implemented."<<endl;
      PRC(i_flag);PRC(m_flag);PRC(r_flag);PRL(v_flag);
    }

    mk_mrv(b, i_flag, id, 
	      m_flag, m, 
	      r_flag, r, 
	      v_flag, v); 
    m_tot += m;
  }

  root->set_mass(m_tot);
  return root;
}

int main(int argc, char ** argv)
{
    bool F_flag = false;
    int i_flag = -1;
    int m_flag = -1;
    int r_flag = -1;
    int v_flag = -1;
    int order = 0;
    bool vctr[] = {false, false, false, false};
    char *filename;

    int n = -1;
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "F:imrvn:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'F': F_flag = true;
		      filename = poptarg;
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 'i': i_flag = order++;
		      break;
	    case 'm': m_flag = order++;
		      break;
	    case 'r': r_flag = order++;
	              vctr[r_flag] = true;
		      break;
	    case 'v': v_flag = order++;
	              vctr[v_flag] = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    cerr << "Selected order: " << endl;
    PRC(order);PRC(i_flag);PRC(m_flag);PRC(r_flag);PRL(v_flag);

    dyn * root = NULL;
    if (F_flag) {
      ifstream infile(filename);
      if (!infile) cerr << "error: couldn't create file "
	                  << filename <<endl;
      cerr << "Reading input from file "<< filename <<endl;

      root = initialize_dyn(infile, n, i_flag, m_flag, r_flag, v_flag);
    }
    else {
      root = initialize_dyn(cin, n, i_flag, m_flag, r_flag, v_flag);
    }

    root->log_history(argc, argv);
    put_dyn(root);
    rmtree(root);
    return 0;
}

#endif

/* end of: mknode.c */

