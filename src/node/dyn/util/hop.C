
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// hop:       After the HOP group-finding algorithm for N-body systems 
////            Eisenstein, D.J. & Hut, P., 1998, ApJ 498, 137
////            
////            input N-body snapshot, output list of clump memberships
////
//// Options:     -i    number the particles sequentially [don't number]
//// Options:     -n    number the neighbors searched for density extremum
//// Options:     -v    Write snapshot at the end

//   Version 1.0:       Simon Portegies Zwart,          Hamilton Aug 2004, 
//   Adopted to starlab:
//   Version 1.1:       Simon Portegies Zwart,          Haarlem March 2005, 

#include "hop.h"

//#ifndef TOOLBOX

dyn* hop::densest_nth_nn(dyn *b) {

  nearest_neighbor *nn = new nearest_neighbor[nn_search];

  for_all_daughters(dyn, b->get_root(), bb) {
      real sep2 = square(bb->get_pos() - b->get_pos());
      for (int i=0; i<nn_search; i++) {
	if (sep2 < nn[i].d_nn_sq) {
	  for (int j=nn_search-1; j>i; j--) {
	    nn[j] = nn[j-1];
	  }
	  nn[i].nn = bb;
	  nn[i].d_nn_sq = sep2;
	  break;
	}
      }
    }
  real d, dmax = -VERY_LARGE_NUMBER;
  int imax;
  for (int i=0; i<nn_search; i++) {
    d = getrq(nn[i].nn->get_dyn_story(), "density");
    if (d>dmax) {
      dmax = d;
      imax = i;
    }
  }

  return nn[imax].nn;
}

void hop::find_clump_center(dyn *b) {

  cerr << " Clumb for i=" <<  b->get_index() <<endl;

  dyn *onnd, *nnd = b;
  do {
    onnd = nnd;
    nnd = densest_nth_nn(onnd);
  }
  while(onnd->get_index()!=nnd->get_index());

  cerr << " id= " << nnd->get_index() << endl;

  add_cluster_center(nnd);
  putiq(b->get_dyn_story(), "hop_clump_center_id", nnd->get_index());

}

void hop::add_cluster_center(dyn* bc) {

  vector<cluster>::iterator ic = cl.begin(); 
  bool known_before = false;
  for (ic = cl.begin(); ic<cl.end(); ic++) {
    if(ic->get_h() == bc) {
      known_before = true;
      break;
    }
  }

  if(known_before) {
    ic->increment(bc);
  }
  else {
    cluster ccl(bc);
    cl.push_back(ccl);
  }
}

void find_primary_cluster(dyn *b) {

  hop h;
  h.set_nn_search(int(3*sqrt((real)b->n_daughters())));
  h.find_primary_cluster(b);
}

void hop::find_primary_cluster(dyn *b) {

  if(!nn_search)
    nn_search = int(3*sqrt((real)b->n_daughters()));

  int nb;
  real mb;

  dyn *bdmax = NULL;
  real d, dmax = 0;
  for_all_daughters(dyn, b, bi) {
    d = getrq(bi->get_dyn_story(), "density");
    if(d>dmax) {
      dmax = d;
      bdmax = bi;
    }
  }

  for_all_daughters(dyn, b, bi) 
    find_clump_center(bi);

}

//#else 

int main(int argc, char ** argv)
{
  //    check_help();

    int nth_neighbor = -1;
    bool use_nsqrt = true;

    extern char *poptarg;
    int c; 
    char* param_string = "in:v";

    bool i_flag = false;
    bool v_flag = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'i': i_flag = true;
		      break;
	    case 'n': nth_neighbor = atoi(poptarg);
	              use_nsqrt = false;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (use_nsqrt)
      cerr << "No nn-th neighbor identified, chose sqrt(n)" << endl;

    int k=12;
    bool first = true;
    dyn* b;
    int i=0;
    hop *h;
    while((b = get_dyn())) {
      cerr << "Read snapshot #"<<i++<<endl;
      if(use_nsqrt)
	nth_neighbor = int(3*sqrt((real)b->n_daughters()));
      h = new hop();
      h->set_nn_search(nth_neighbor);
      
      h->find_primary_cluster(b);
      h->put();
      delete h;

      if(v_flag) {
	put_dyn(b);
      }
    }
}

//#endif
/* end of: hop.C */
