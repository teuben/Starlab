#ifndef     STARLAB_HOP_H
#   define  STARLAB_HOP_H

#include "dyn.h"
#include "vector"

real d_nn_sq, d_nn_sq2, d_nn_sq3;
static dyn *nn2, *nn3;

class nearest_neighbor {
public:
  dyn* nn;
  real d_nn_sq;
  nearest_neighbor() {
    nn = NULL;
    d_nn_sq = VERY_LARGE_NUMBER;
  }
};

class cluster {
 protected:

  dyn *h;

  int nstar;
  real mass;
  vec pos;
  vec vel;

 public:

  cluster() {
    h=NULL;
    nstar = 0;
    mass = 0;
  }
  cluster(dyn *b) {
    h=b;
    nstar = 1;
    mass = b->get_mass();
    pos = b->get_pos();
    vel = b->get_vel();
  }
  ~cluster() {}

  dyn *get_h() {return h;}
  void set_h(dyn *b) {h=b;}
  void set_nstar(int n) {nstar = n;}
  void set_mass(real m) {mass = m;}
  void set_pos(vec p) {pos = p;}
  void set_vel(vec p) {vel = p;}

  void increment(dyn *b) {
    nstar++;
    mass += b->get_mass();
  }

  void put(ostream &s = cerr) {
    s << "Cluster (id= " << h->get_index() 
      << "), N= " << nstar 
      << ", M= " << mass << endl; 
    s << "   pos= " << pos << "   vel= " << vel << endl;
  }

};


class hop {
 protected:

  int nn_search;
  vector<cluster> cl;
 public:

  hop() {}
  ~hop() {}

  void set_nn_search(int n) {nn_search = n;}

  dyn *densest_nth_nn(dyn *b);
  void find_clump_center(dyn *b);
  void find_primary_cluster(dyn *b);
  void add_cluster_center(dyn* bc);

  void put(ostream &s = cerr) {
    vector<cluster>::iterator ic;
    for (ic = cl.begin(); ic<cl.end(); ic++) {
      ic->put(s);
    }
  }

  void put_snap(ostream &s) {
    vector<cluster>::iterator ic;
    for (ic = cl.begin(); ic<cl.end(); ic++) {
      put_dyn(ic->get_h());
    }
  }
};

#endif    //STARLAB_HOP_H

