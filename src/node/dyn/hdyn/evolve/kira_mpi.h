#include <hdyn.h>
extern int mpi_communicator;
extern int mpi_myrank;
extern int mpi_nprocs;
void mpi_sum_list(hdyn * next_nodes[], int n_next, MPI_Comm mpi_communicator);
