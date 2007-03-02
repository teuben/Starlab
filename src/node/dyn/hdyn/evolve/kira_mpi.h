#include <hdyn.h>
extern MPI_Comm mpi_communicator;
extern int mpi_myrank;
extern int mpi_nprocs;
extern long long mpi_allreduce_count;
extern long long mpi_allreduce_len;
extern long long mpi_allgather_count;
extern long long mpi_allgather_len;
void mpi_sum_list(hdyn * next_nodes[], int n_next, MPI_Comm mpi_communicator);
extern double *mpi_k_c_t_l_a_j_timings;
extern long long *hdyn_ev_counts;
extern long long hdyn_ev_count;

extern hdyn* my_start_daughter;
extern int my_daughter_count; 
extern int my_start_count; 
extern int my_end_count; 
extern int num_daughters;
