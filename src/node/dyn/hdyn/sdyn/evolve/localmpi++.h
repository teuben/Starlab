#ifndef STARLAB_LOCALMPI_H
#   define STARLAB_LOCALMPI_H

//#define HAS_MPI

#ifndef HAS_MPI
class MPI {
  public:
  class COMM_WORLD {};
  int Get_rank() {return 0;}
};


//void finalize_MPI() {}

struct MPI_Datatype {
  // empty
};

//void initialize_MPI(int &myid, int &nproc) {}
#endif

#endif







