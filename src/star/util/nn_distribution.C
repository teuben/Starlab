
//// nn:  find nearest neighbor.
////
////          bladiebla
////
//// Options:     -c    add a comment to the output snapshot [false]

//-----------------------------------------------------------------------------
//   version 1:  Aug 2001   Simon Portegies Zwart, MIT
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "dyn.h"
#include "assert.h"

#ifndef TOOLBOX

#else

void nn_distribution(dyn * b) {

    b->flatten_node();
    int n = b->n_leaves();
    if(n<=0) {
	for_all_leaves(dyn, b, bb)
	  n++;
    }

    real *nn = new real[n];
    int *indexA = new int[n];
    int *indexB = new int[n];
    for(int j=0; j<n; j++) {
      nn[j] = VERY_LARGE_NUMBER;
    }

    int i=-1;
    real r;
    for_all_daughters(dyn, b, bi) {
      i++;
      indexA[i] = bi->get_index();
      for_all_daughters(dyn, b, bj) {
	if(bi!=bj) {
	    r = abs(bi->get_pos()-bj->get_pos());
	  if(r<=nn[i]) {
	    nn[i] = r;
	    indexB[i] = bj->get_index();
	  }
	}
      }
    }

    if(indexA[0]==-1)
      for(int j=0; j<n; j++) 
	cout << nn[j] << endl;
    else
      for(int j=0; j<n; j++)
	cout << indexA[j] << "\t" 
	     << indexB[j] << "\t" << nn[j] << endl;

    delete [] nn;
    delete [] indexA;
    delete [] indexB;
}

//-----------------------------------------------------------------------------
//  main  --  driver to use  compute_luminosity_radii() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = false;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }    

    dyn* b;
    while (b = get_dyn()) {

      if (c_flag == TRUE)
	b->log_comment(comment);

      b->log_history(argc, argv);

      nn_distribution(b);
      rmtree(b);
    }
}

#endif


