
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// renumber:            Renumber stars in a specific order
////
//// Options:   -c c      add a comment to the output snapshot [false]
////            -I/i      start numbering number               [1]
////            -M        Renumber the stars in order of mass
////                      (highest mass=I/lowest mass=i)       [false]
////

//// Options:

//   version 1:  Jan 1993   Piet Hut
//	 	 Feb 2001   Simon Portegies Zwart

#include "node.h"

#ifndef TOOLBOX

typedef  struct {
    node* str;
    real  mass;
} nm_pair, *nm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_mass  --  compare the masses of two particles
//-----------------------------------------------------------------------------

local int compare_mass(const void * pi, const void * pj)
{
    if (((nm_pair_ptr) pi)->mass < ((nm_pair_ptr) pj)->mass)
        return(1);
    else if (((nm_pair_ptr)pi)->mass > ((nm_pair_ptr)pj)->mass)
        return(-1);
    else
        return(0);
}

void renumber(node* b, int istart, bool mass_order,
	      bool name_nodes) {

    int i;
    if(!mass_order) {

      i = istart;
      for_all_leaves(node, b, bj)
	bj->set_label(i++);
    }
    else {

      // Renumber the stars in order of mass.
      // Highest mass gets smallest number (strange choise, but).

      int n = b->n_leaves();
      nm_pair_ptr nm_table = new nm_pair[n];
      if (nm_table == NULL) {
	cerr << "renumber: "
	     << "not enough memory left for nm_table\n";
	return;
      }

      i=0;
      for_all_daughters(node, b, bi) {
	nm_table[i].str = bi;
	nm_table[i].mass = bi->get_mass();
	i++;
      }

      qsort((void *)nm_table, (size_t)n, sizeof(nm_pair), compare_mass);

      for (i=0; i<n; i++) {
	nm_table[i].str->set_index(istart+i);
      }
      delete []nm_table;

    }

    char tmp[128];
    if(name_nodes)
      for_all_leaves(node, b, bj) {
      PRL(bj->get_index());
      if (bj->get_index() >= 0) {
	sprintf(tmp, "%d", bj->get_index());
	bj->set_name(tmp);
      }
    }
}

#else

void main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;

    bool M_flag = false;
    bool N_flag = false;
    int istart = 1;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "MmNI:i:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'i':
	    case 'I': istart = atoi(poptarg);
		      break;
	    case 'm':
	    case 'M': M_flag = true;
		      break;
	    case 'N': N_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}


    node* b;
    b = get_node(cin);
    b->log_history(argc, argv);

    renumber(b, istart, M_flag, N_flag);

    put_node(cout, *b);
    rmtree(b);

}
#endif

