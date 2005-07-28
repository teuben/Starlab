
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


// Functions associated with force evaluations within kira.
//
// Externally visible functions:
//
//	void calculate_acc_and_jerk_for_list
//	void calculate_acc_and_jerk_on_all_top_level_nodes
//	void kira_synchronize_tree
//	void initialize_system_phase2
//	void clean_up_kira_ev

// Substantially rewrote code and removed timing and some other
// debugging #ifdefs -- too confusing, and no longer relevant
// after the rewrite.
//					      Steve 4/03, 8/03
//
// Major reorganization to remove compile-time GRAPE selection.
//					      Steve 6/04

#include "hdyn.h"
#include "kira_timing.h"

#ifdef USEMPI
#    include <mpi.h>
#    include "kira_mpi.h"
#    include "hdyn_inline.C"
#endif

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_kira_ev
#   undef T_DEBUG
#endif


//-------------------------------------------------------------------------
//
// NOTE:  Currently, kira_calculate_top_level_acc_and_jerk() is a simple
// compile-time switch between the various (GRAPE and non-GRAPE) methods
// available to compute the acc and jerk (due to internal forces only)
// on the particles listed in next_nodes[].
//
// Alternatively, it could be replaced by a function pointer leading
// directly to the functions involved, set at startup based on the
// configuration of the system.  Implemented and then removed, but
// most of the commented-out code can still be found...
//
// The function declaration is:
//
//	int kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
//						  int n_next,
//					  	  xreal time,
//					  	  bool restart_grape);
//
// where restart_grape may be ignored or used for other purposes in future
// applications.  It indicates a change in top-level tree structure that
// requires an internal reset in the GRAPE arrays.  However, the GRAPE-6
// code now sets an internal flag to convey the information, so only the
// GRAPE-4 code actually uses this variable.  The functions's return
// value is simply the number of top-level nodes in the system.
//
// Currently, the list of possible functions that can be used is
//
//	top_level_acc_and_jerk_for_list()	below
//	grape4_calculate_acc_and_jerk()		in hdyn_grape4.C
//	grape6_calculate_acc_and_jerk()		in hdyn_grape6.C
//
// Note the acc_function_ptr typedef in hdyn.h.
//
// Currently, the only call comes from calculate_acc_and_jerk_for_list(),
// and restart_grape is always set false there after the call.
//
//-------------------------------------------------------------------------


#ifdef USEMPI

int mpi_communicator;
int mpi_myrank;
int mpi_nprocs;

typedef struct { real val; int id; } real_int_type;

void hpsortreal_int(int n,  real_int_type ra[], int ia[])
{
	int  i,ir,j,l;
	real_int_type rra;
	int iia;

	ra--;  /* this is a numerical recipes function, */
	ia--;  /* so indexing starts with 1 (sic!)      */

	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
			iia = ia[l];
		} else {
			rra=ra[ir];
			iia=ia[ir];
			ra[ir]=ra[1];
			ia[ir]=ia[1];
			if (--ir == 1) {
				ra[1]=rra;
				ia[1]=iia;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j].id < ra[j+1].id) j++;
			if (rra.id < ra[j].id) {
				ra[i]=ra[j];
				ia[i]=ia[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
		ia[i]=iia;
	}
}

void mpi_sum_list(hdyn * next_nodes[], int n_next, 
                                   MPI_Comm mpi_communicator)
{
  // First a test if all n_next's are equal
  int n_max, n_min;
  MPI_Allreduce(&n_next,&n_max,1,MPI_INT,MPI_MAX,mpi_communicator);
  MPI_Allreduce(&n_next,&n_min,1,MPI_INT,MPI_MIN,mpi_communicator);
  if (n_max != n_min)
  {
    cout << ":"<<mpi_myrank<<": Catastrophic error in "<<
      __FILE__<<":"<<__LINE__<<" n_next's are not equal"
      ", my n_next is:"<<n_next<<endl;
    MPI_Finalize();
    exit(1);
  }
  {
    // s/rbuf for sending/receiving acc, jerk and pot
    real *sbuf = new real[7*n_next];
    real *rbuf = new real[7*n_next];
    int k=0;
    vec a;
    for (int i=0; i<n_next; i++)
      {
        if (next_nodes[i]){
  	  a = next_nodes[i] -> get_acc();
	  sbuf[k++] = a[0];
	  sbuf[k++] = a[1];
	  sbuf[k++] = a[2];
	  a = next_nodes[i] -> get_jerk();
	  sbuf[k++] = a[0];
	  sbuf[k++] = a[1];
	  sbuf[k++] = a[2];
	  sbuf[k++] = next_nodes[i] -> get_pot();
        }
        else
        {
 	  cerr << __FILE__ << ":"<< __LINE__ << 
		 ": NULL next_node " << i << endl;
        }
      }
    int mpi_real = MPI_DOUBLE;
    if (sizeof(real) != sizeof(double) )
         mpi_real = MPI_FLOAT;
    // add valuse of acc, jerk and pot in sbuf giving rbuf
    cout << ":"<<mpi_myrank<<": TD: allreducing n_next: "<<n_next<<endl;
    MPI_Allreduce(&sbuf[0],&rbuf[0],n_next*7,mpi_real,MPI_SUM,
            mpi_communicator);  
    cout << ":"<<mpi_myrank<<": TD: allreduced n_next: "<<n_next<<endl;
    // copy summed acc, jerk and pot values to the nodes
    k=0;
    for (int i=0; i<n_next; i++)
    {
      if (next_nodes[i]){
        next_nodes[i] -> set_acc_and_jerk_and_pot(
 	  vec(rbuf[k  ],rbuf[k+1],rbuf[k+2]),
	  vec(rbuf[k+3],rbuf[k+4],rbuf[k+5]),
	     rbuf[k+6]);
      }
      k += 7;
    }
    delete [] sbuf;
    delete [] rbuf;
  }
   // handling of nn and coll
   // 
   // Every process has it's own notion of who are the nn's and
   // the coll's. For each element of next_nodes we will
   // determine what is the minimum d_nn_sq and what is the
   // id of the nn.
   // 
   // MPI_All_reduce can be used for that purpose when the
   // datatype MPI_DOUBLE_INT is used. It consists of a double,
   // followed by an int. The allreduce finds the minimum of the doubles
   // and returns the double and the corrsponding int. 
   // For the double we use d_nn_sq, for the int we use nn_id.
   //
   // The same story for coll.  
   //
   // Eventually, the reduce for acc, jerk, pot and nn,coll can be done
   // in one call. TODO
  {
    // sbuf for sending nn and coll stuff,
    // rbuf for receiving
    real_int_type *sbuf = new real_int_type[2*n_next];
    real_int_type *rbuf = new real_int_type[2*n_next];

    for (int i=0; i<n_next; i++)
    {
      if (next_nodes[i]){
	sbuf[i+i].val = next_nodes[i] -> get_d_nn_sq();
	sbuf[i+i].id = next_nodes[i] -> get_nn_id();
	sbuf[i+i+1].val = next_nodes[i] -> get_d_coll_sq();
	sbuf[i+i+1].id = next_nodes[i] -> get_coll_id();
      }
      else
      {
	cout << ":" << mpi_myrank << ":" << __FILE__ << ":"<< __LINE__ << 
	       ": NULL next_node " << i << endl;
      }
    }
    for (int k=0; k<2*n_next; k+=2)
    {
	cout << ":" << mpi_myrank << ": TD:" << k << 
	  " val: " << sbuf[k  ].val << " id: " << sbuf[k  ].id << 
	  " val: " << sbuf[k+1].val <<  " id: " <<sbuf[k+1].id << endl;
    }

    int mpi_real_int = MPI_DOUBLE_INT;
    if (sizeof(real) != sizeof(double) )
         mpi_real_int = MPI_FLOAT_INT;

    // find minimum of distances, giving rbuf
    cout << ":"<<mpi_myrank<<": TD: allreducing coll nn n_next: "
               <<n_next<<endl; 
    MPI_Allreduce(&sbuf[0],&rbuf[0],n_next*2,mpi_real_int,MPI_MINLOC,
            mpi_communicator);  
    cout << ":"<<mpi_myrank<<": TD: ready allreducing coll nn n_next: "
               <<n_next<<endl; 
    //
    // Now, in rbuf[i].id we find the nn_id (i even)
    // and                            coll_id (i odd)
    //
    // Find now the addresses of the corresponding nodes and 
    // put these in the nn and coll fields of the next_nodes array
    //

    //
    // We loop over the whole universe, looking if the id
    // (i.e. the number of the daughter visited)
    // is in rbuf. For efficiency, it is best to sort rbuf
    // using the id field as key. 
    // We have to remember if it is a nn or a coll and
    // we have to keep track of original place in the
    // next_nodes array. 
    // We define an array, dimension 2*n_next which contains the
    // indexes into the next_nodes array:
    // 0 0 1 1 2 2 3 3 .....
    // Multiply this by 2 and add 1 for the odd places:
    // 0 1 2 3 4 5 6 7 .....
    //
    int * indices = new int[2*n_next];
    for (int i=0; i<2*n_next; i++)
    {
      indices[i] = i;
    }
    for (int k=0; k<2*n_next; k++)
    {
	printf(":%d: TD:k%d: rbuf val %f id %d index %d\n",mpi_myrank,k,rbuf[k].val,
	    rbuf[k].id,indices[k]);
    }

    printf(":%d: TD: calling hpsort\n",mpi_myrank);
    hpsortreal_int(n_next*2, &rbuf[0], &indices[0]);
    printf(":%d: TD: called hpsort\n",mpi_myrank);
    for (int k=0; k<2*n_next; k++)
    {
	printf(":%d: TD:k%d: rbuf sort val %f id %d index %d\n",mpi_myrank,k,rbuf[k].val,
	    rbuf[k].id,indices[k]);
    }
    int k=0;
    int l;
    int id = rbuf[k].id;
    int counter = 0;

    printf(":%d: TD: entering for_all_daughters loop\n",mpi_myrank);
    for_all_daughters(hdyn,next_nodes[0]->get_root(),bi)
    {
      int daughter_id = bi->get_node_id();
      // Just to test if numbering of the universe did not change
      if (daughter_id != counter)
      {
        cout << ":"<<mpi_myrank<<": Catastrophic error in "<<
          __FILE__<<":"<<__LINE__<<" daughter_id: "<<daughter_id<<
	  " counter: "<<counter<<endl;
	cout << ":"<<mpi_myrank<<" these should be equal"<<endl;
        MPI_Finalize();
        exit(1);
      }
      counter++;
      // printf (":%d: TD: examining id:%d daughter_id:%d \n",mpi_myrank,id,daughter_id);
      if (id != daughter_id) continue;
      // found an entry in rbuf with an id equal to the id
      // of the current daughter
      // Now, determine if this is a nn or a coll id, and act as appropriate

      while (id == daughter_id)
      {
	l = indices[k]/2;
	if (indices[k] & 1) // odd, so coll
	{
	  cout << ":"<<mpi_myrank<<": TD: setting coll "<<l<<" "<<rbuf[k].val<<endl;
	  next_nodes[l] -> set_coll(bi);
	  next_nodes[l] -> set_d_coll_sq(rbuf[k].val);
	}
	else // even, so nn
	{
	  cout << ":"<<mpi_myrank<<": TD: setting nn "<<l<<" "<<rbuf[k].val<<endl;
	  next_nodes[l] -> set_nn(bi);
	  next_nodes[l] -> set_d_nn_sq(rbuf[k].val);
	  // TODO:    also set_d_nn_sq(rbuf[k].val); ??
	}
	// TODO : take in consideration the possibility of next_node[l]
	// == NULL?
	k++;
	if (k >= 2*n_next)
	  break;
	id = rbuf[k].id;
      }
      if (k >= 2*n_next)
	break;
    }
    delete []sbuf;
    delete []rbuf;
    delete []indices;
  }
}
// USEMPI

#endif


local inline int _kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
							int n_next,
							xreal time,
							bool restart_grape)
{
    if (n_next <= 0) return 0;

    int n_top = 0;
    unsigned int config = next_nodes[0]->get_config();

    switch (config) {

	case 0:	n_top = top_level_acc_and_jerk_for_list(next_nodes, n_next,
							time);
		break;

	case 1:	n_top = grape4_calculate_acc_and_jerk(next_nodes, n_next,
						      time, restart_grape);
		break;

	case 2:	n_top = grape6_calculate_acc_and_jerk(next_nodes, n_next,
						      time, restart_grape);
		break;
    }

    return n_top;
}

// (Use the inline version within this file...)

int kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
					  int n_next,
					  xreal time,
					  bool restart_grape)
{
    return _kira_calculate_top_level_acc_and_jerk(next_nodes, n_next,
						  time, restart_grape);
}



// ************************************************************************
// ***  This function is a candidate for threading if n_next is large.  ***
// ***  Hmmm.  Simultaneous force calculations here may cause cache     ***
// ***  problems in top_level_node_real_force_calculation.              ***
// ***  Better to thread there? -- Probably too low-level??             ***
// ***                                                    Steve, 4/05   ***
// ************************************************************************

local inline int std_top_level_acc_and_jerk_for_list(hdyn **next_nodes,
						     int n_next,
						     xreal time)  // not used
{
    int n_top = 0;

    for (int i = 0; i < n_next; i++) {
	hdyn *bi = next_nodes[i];
	if (bi->is_top_level_node())
	    n_top += bi->top_level_node_real_force_calculation(); // hdyn_ev.C
    }

    return n_top;
}

#ifndef HAVE_LIBPTHREAD

int top_level_acc_and_jerk_for_list(hdyn **next_nodes,
				    int n_next,
				    xreal time)  	// not used
{
    return std_top_level_acc_and_jerk_for_list(next_nodes, n_next, time);
}


#else

// Threaded code:

#include <pthread.h>

typedef struct {hdyn **node_list;
		int n_list;
	        int thread;} arg_list;
typedef struct {int n_top;
		real cpu;} thread_stats;
static thread_stats *all_stats = NULL;

void *thread_top_level_acc_and_jerk_for_list(void *arg)
{
    arg_list *a	= (arg_list*)arg;
    hdyn *b = a->node_list[0];
    int n_threads = b->get_n_threads();

    int n_top = 0;
    real cpu = cpu_time();

    int t = a->thread;
    for (int i = t; i < a->n_list; i += abs(n_threads)) {
	hdyn *bi = a->node_list[i];
	if (bi->is_top_level_node()) {
	    if (n_threads > 0)
		n_top += bi->top_level_node_real_force_calculation();
	    else
		n_top += 1;
	}
    }

    all_stats[t].n_top = n_top;
    all_stats[t].cpu = cpu_time() - cpu;	// CPU time for this thread

    if (t == 0)
	return (void*)1;
    else
	pthread_exit((void*)1);			// anything non-NULL will do
}

int top_level_acc_and_jerk_for_list(hdyn **next_nodes,
				    int n_next,
				    xreal time)			  // not used
{
    if (n_next <= 0)

	return 0;

    else if (n_next <= 1 || next_nodes[0]->get_n_threads() == 0)

	return std_top_level_acc_and_jerk_for_list(next_nodes,
						   n_next, time);

    else {

	hdyn *b = next_nodes[0];
	int n_threads = b->get_n_threads();

	// Use threads to parallelize the calculation.

	if (!all_stats) {
	    all_stats = new thread_stats[abs(n_threads)];
	    if (!all_stats) err_exit("Can't create all_stats for threads.");
	}

	pthread_t threads[abs(n_threads)];
	arg_list a[abs(n_threads)];

	// Start up abs(n_threads)-1 threads.  Main program will be 0.

	for (int t = 1; t < abs(n_threads); t++) {

	    // Place argument in arg_list for transmittal to the thread.

	    a[t].node_list = next_nodes;
	    a[t].n_list = n_next;
	    a[t].thread = t;

	    int return_code
		= pthread_create(threads+t, NULL,
				 thread_top_level_acc_and_jerk_for_list,
				 (void *)(a+t));
	    if (return_code) {
		PRL(return_code);
		err_exit("thread creation error");
	    }
	}

	int n_top = 0;

	// We are thread 0.  This is ugly...

	a[0].node_list = next_nodes;
	a[0].n_list = n_next;
	a[0].thread = 0;
	void *stat = thread_top_level_acc_and_jerk_for_list((void*)a);
	if (stat)
	    n_top += all_stats[0].n_top;

	// Wait for the threads to finish and combine the results.

	for (int t = 1; t < abs(n_threads); t++) {
	    void *stat;
	    int rc = pthread_join(threads[t], &stat);

	    if (stat) {
		n_top += all_stats[t].n_top;
		b->set_thread_cpu(b->get_thread_cpu() + all_stats[t].cpu);
	    }
	}

	if (n_threads > 0)
	    return n_top;
	else
	    return std_top_level_acc_and_jerk_for_list(next_nodes,
						       n_next, time);
    }
}

#endif



int calculate_acc_and_jerk_for_list(hdyn **next_nodes,
				    int  n_next,
				    xreal time,
				    bool exact,
				    bool tree_changed,
				    bool &reset_force_correction, // no longer
								  // used...
				    bool &restart_grape)
{
    // Note that this function uses the same function calls as
    // calculate_acc_and_jerk_on_top_level_node, but in a different
    // order: all prologue calls are performed first, then all
    // top-level (i.e. GRAPE) force calculations, then all epilogue
    // calls.

    if (!next_nodes[0] || !next_nodes[0]->is_valid()) return 0;	  // unnecessary

    hdyn *b = next_nodes[0]->get_root();

    kira_counters *kc = b->get_kira_counters();
    kira_diag *kd = b->get_kira_diag();
    kira_options *ko = b->get_kira_options();

    bool ignore_internal = b->get_ignore_internal();

    xreal sys_t = b->get_system_time();

    // Note that time and system_time should be the same...

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 1 << endl;
	int p = cerr.precision(HIGH_PRECISION);
	PRI(7); PRC(n_next); PRC(exact); PRC(time); PRL(time-sys_t);
	cerr.precision(p);
    }
#endif

    // Explicitly split the calculation into top-level and low-level
    // nodes (Steve, 4/03).  Top-level nodes run from 0 to n_top-1;
    // low-level nodes from n_top to n_next-1.
    //
    // The first loop here does the splitting, at the same time
    // performing the first part of the force calculation.

    bool print = kd->kira_ev;
    int n_top = 0;

    for (int i = 0; i < n_next; i++) {

	hdyn *bi = next_nodes[i];
	bool top = false;

	if (bi->get_parent() == b) {

	    // This is a top-level node.

	    top = true;

	    predict_loworder_all(bi, sys_t);

	    // Prologue operations.  DO NOT clear the interaction here,
	    // as acc and jerk may be needed to update the GRAPE before
	    // the new force calculation.

	    // bi->clear_interaction();

	    // Note that top_level_node_prologue_for_force_calculation()
	    // does the *entire* calculation in the case exact = true.
	    // For exact = false it does almost nothing...

	    if (!ignore_internal)
		bi->top_level_node_prologue_for_force_calculation(exact);

	} else {

	    // This is a low-level node.
	    // Predict the entire clump containing node bi (Steve, 8/98):

	    if (!ignore_internal)
		predict_loworder_all(bi->get_top_level_node(), sys_t);
	}

	if (ignore_internal && (!top || !exact)) {

	    // Set flags in case of no internal forces...

	    bi->set_nn(bi);
	    bi->set_coll(bi);
	    bi->set_d_nn_sq(VERY_LARGE_NUMBER);
	}

	bi->inc_steps();

	// All initial operations on bi are complete.
	// Reorder the list if necessary.

	if (top) {
	    next_nodes[i] = next_nodes[n_top];
	    next_nodes[n_top++] = bi;
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 2
	     << endl << flush;
    }
#endif

    // Retain restart_grape for compatibility with other functions, but
    // in the GRAPE-6 version we actually transmit the information via
    // the static restart_grape flag.  The reason is that (as of 3/05)
    // there now are two functions which may need to reset the internal
    // data structures, and it is too complicated to carry the flag
    // through the many levels before we get to the new call.  The flag
    // will be unset as soon as it is tested and action taken.
    //
    // In general, the use of static flags is probably a better way to
    // handle this sort of communication (Steve, 3/05).

    // Synchronize the old and new flags.  Retain the old version until we
    // are sure the new one works, then we can remove all restart_grape
    // references from the GRAPE version of calculate_acc_and_jerk()...
    // 							(Steve 3/05)

    if (tree_changed) restart_grape = true;
    if (restart_grape) b->set_restart_grape_flag();

    // Complete the top-level force calculation (now top-level nodes
    // run from 0 to n_top-1):

    if (n_top > 0) {

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 3
		 << endl << flush;
	}
#endif

	int n_force = 1;

	// Calculate top-level forces.  Note new return value from
	// kira_calculate_top_level_acc_and_jerk() (Steve, 4/03).

	// Top_level force calculation actually uses
	//
	//	grape*_calculate_acc_and_jerk()
	//
	// or
	//
	//	top_level_node_real_force_calculation(), via
	//	top_level_acc_and_jerk_for list()
	//
	// as appropriate.

	if (!exact && !ignore_internal) {


// 	    acc_function_ptr get_acc_and_jerk
// 		= b->get_kira_calculate_top_level_acc_and_jerk();
//
// 	    n_force = get_acc_and_jerk(next_nodes,
// 				       n_top, time,
// 				       restart_grape);

#ifdef USEMPI
	  // fill in node_id in all nodes. We need it later
	  int node_id = 0;
          for_all_daughters(hdyn,next_nodes[0]->get_root(),bi)
          {
            bi->set_node_id(node_id);
	    node_id++;
          }
#endif

	    n_force = _kira_calculate_top_level_acc_and_jerk(next_nodes,
							     n_top, time,
							     restart_grape);

#ifdef USEMPI

	  // In the MPI case, next_nodes contain the values of acc, jerk and
	  // pot calculated using only a part of the universe b.
	  // mpi_sum_list knows how to add them up and give each process the
	  // same values.

            cout << ":" << mpi_myrank <<":"
		 << " TD: calling mpi_sum_list n_top:"<< n_top 
	         << " n_force:" << n_force << endl;

            mpi_sum_list(next_nodes, n_top, mpi_communicator);

	    // wwvv this is the place where the rest of the calculations
	    // are to be done: the calculations that were skipped in the
	    // for_all_daughters loop in hdyn_ev.C
	    // Does not work. A rethinking of the way the calculations
	    // are to be done is necessary.
	    //
	    // It would be nice, if the code leading to the for_all_daughters
	    // loop in hdyn_ev.C could be separated: one call that only 
	    // involves the calculations of the 'next_nodes' array, and
	    // thereafter a call that updates the universe.
	    //
	    // In the current code, in hdyn_ev.C, after the calculations on
	    // one element of 'next_nodes', the rest of the universe is updated,
	    // which is not efficient in the parallel version, maybe even
	    // not doable. In the parallel version, the sideeffects of the
	    // calculations on the 'next_nodes' should be taken care off
	    // after the MPI_Allreduce (... sum ...) of the 'next_nodes'.
	    //
#if 0
	    _kira_calculate_top_level_acc_and_jerk(next_nodes,
							     n_top, time,
							     restart_grape);
#endif

#endif

	    // Note that we now clear restart_grape explicitly here.

	    restart_grape = false;
	    b->clear_restart_grape_flag();	// should be unnecessary
	    
	}

//	if (sys_t >= 44.15329 && sys_t <= 44.1533) {
//	    cerr << "after top-level acc_and_jerk" << endl;
//	    pp3(next_nodes[0]->get_top_level_node());
//	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 4
		 << endl << flush;
	}
#endif

	// Complete calculation of top-level accs and jerks by correcting
	// for C.M. interactions.

	if (!exact) {

	    // Note: correct_acc_and_jerk() now checks list membership to
	    // determine if correction is needed.

	    // The new version of correct_acc_and_jerk appears to work, but
	    // retain the possibility of reverting to the old version until we
	    // are sure there are no problems.  Results of the two versions
	    // are similar, but *not* identical.

	    // On_integration_list flags are for use by correct_acc_and_jerk(),
	    // to ensure that we only correct interactions between objects on
	    // the list.  Should be unset immediately on return, but see the
	    // note below.

	    n_force--;

	    // Epilogue operations.

	    for (int i = 0; i < n_top; i++) {

		hdyn *bi = next_nodes[i];

		// Epilogue force calculation mostly performs CM corrections
		// and cleans up the perturber list.

		bi->inc_direct_force(n_force);		// direct force counter
		bi->top_level_node_epilogue_force_calculation();
		bi->set_on_integration_list();
	    }

	    if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list)

		correct_acc_and_jerk(b,			// old version
				     reset_force_correction);
	    else

		correct_acc_and_jerk(next_nodes,	// new version
				     n_top);

	    // Don't unset the "on_integration_list" flags here, because
	    // the loop through memory may cost more than it is worth.
	    // *** Must unset these flags in the calling function. ***

	    // for (int i = 0; i < n_top; i++)
	    //     next_nodes[i]->clear_on_integration_list();

	}

#ifdef T_DEBUG
	if (IN_DEBUG_RANGE(sys_t)) {
	    cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 5
		 << endl << flush;
	}
#endif

	// Add external forces, if any, to top-level nodes.

	if (b->get_external_field() > 0) {

	    for (int i = 0; i < n_top; i++) {

		hdyn *bi = next_nodes[i];

		real pot = bi->get_pot();
		vec acc = bi->get_acc(), jerk = bi->get_jerk();

		// (Really only need acc to be correctly set...)

		get_external_acc(bi, bi->get_nopred_pos(), bi->get_nopred_vel(),
				 pot, acc, jerk);
		bi->set_pot(pot);
		bi->set_acc(acc);
		bi->set_jerk(jerk);
	    }
	}
    }

    // Everything is done for top-level nodes.  Complete the low-level
    // force calculations.  Do this last, so any recomputation of
    // top-level perturber lists has already occurred.

    if (!ignore_internal) {

	for (int i = n_top; i < n_next; i++) {

	    hdyn *bi = next_nodes[i];
	    if (!bi->get_kepler()) {

		// Compute the forces (perturbed motion).

		if (print) {
		    cerr << "\nComputing force on low-level node "
			 << bi->format_label() << endl;
		    pp3(bi);
		    PRC(bi->get_system_time());
		    PRC(bi->get_time()); PRL(bi->get_t_pred());
		    PRL(bi->get_pred_pos());
		    PRL(bi->get_pred_vel());
		}

		hdyn *sister = bi->get_younger_sister();
		if (!sister) {
		    sister = bi->get_elder_sister();
		    if (!sister) continue;		  // really an error...
		}

		bi->clear_interaction();

		// Doing these steps here allows us to bypass
		// calculate_acc_and_jerk and go directly to
		// calculate_acc_and_jerk_on_low_level_node.

		bi->set_d_coll_sq(VERY_LARGE_NUMBER);
		bi->set_coll(NULL);
		if (sister) sister->set_d_coll_sq(VERY_LARGE_NUMBER);

		bi->calculate_acc_and_jerk_on_low_level_node();
		kc->pert_step++;

		if (print) {
		    cerr << "after..."<<endl;
		    pp3(bi);
		    PRL(bi->get_acc());
		    PRL(bi->get_binary_sister()->get_acc());
		    PRL(bi->get_pos());
		    PRL(bi->get_binary_sister()->get_pos());
		}
	    }
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: calculate_acc_and_jerk_for_list " << 6
	     << endl << flush;
    }
#endif

//   if (sys_t >= 44.15329 && sys_t <= 44.1533) {
//     cerr << "after low-level acc_and_jerk" << endl;
//     pp3(next_nodes[0]->get_top_level_node());
//   }

    return n_top;
}



void calculate_acc_and_jerk_on_all_top_level_nodes(hdyn * b)
{
    int n_top = b->n_daughters();
    hdynptr * list = new hdynptr[n_top];

    int i_top = 0;
    for_all_daughters(hdyn, b, bb)
	list[i_top++] = bb;

    bool reset_force_correction = true; 	// (no longer used)
    bool restart_grape = true;

    calculate_acc_and_jerk_for_list(list, i_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    false,	// called before tree changes
				    reset_force_correction, // obsolete
				    restart_grape);
    delete [] list;
}

void calculate_acc_and_jerk_on_top_level_binaries(hdyn * b)
{
    int n_top = b->n_daughters();
    hdynptr * list = new hdynptr[n_top];

    int i_top = 0;
    for_all_daughters(hdyn, b, bb)
	if (bb->is_parent()) list[i_top++] = bb;

    bool reset_force_correction = true; 	// (no longer used)
    bool restart_grape = true;

    calculate_acc_and_jerk_for_list(list, i_top,
				    b->get_system_time(),
				    false,	// usually what we want...
				    true,
				    reset_force_correction, // no longer used
				    restart_grape);
    delete [] list;
}



void kira_synchronize_tree(hdyn *b,
			   bool sync_low_level)		// default = false
{
    // GRAPE replacement for synchronize_tree().  Synchronize all
    // top-level nodes.  Called from integrate_list() in kira_ev.C and
    // hdyn::merge_nodes().  For now, at least, the entire algorithm
    // from hdyn_ev.C is repeated here.

    if (b->has_grape()) {

	// Code is similar to that in integrate_list(), but only top-level
	// nodes are considered and we don't check for errors in function
	// correct_and_update.  Possibly should merge this with (part of)
	// integrate_list() and drop synchronize_tree() completely.
	//						     (Steve, 1/02)

	// Make a list of top-level nodes in need of synchronization.
	// Generally interested in recomputation of acc and jerk, so
	// probably don't need to treat low-level nodes.  Default is
	// not to touch them.  Note that, even if sync_low_level is
	// true, we still won't synchronize unperturbed binaries.

	xreal sys_t = b->get_system_time();

#if 0
	cerr << endl
	     << "synchronizing tree using GRAPE at time " << sys_t
	     << endl << flush;
#endif

	// Note: no need to set time steps here, and in fact GRAPE will
	// complain if j-particle times and timesteps are not consistent.

	int n_next = 0, n_top = 0;
	for_all_daughters(hdyn, b, bi) {
	    n_top++;
	    if (bi->get_time() < sys_t) {
		// unnecessary and bad:
		// bi->set_timestep(sys_t - bi->get_time());
		n_next++;
	    }
	}

	hdyn **next_nodes = new hdynptr[n_next];
	n_next = 0;
	for_all_daughters(hdyn, b, bi)
	    if (bi->get_time() < sys_t) next_nodes[n_next++] = bi;

	// Integrate all particles on the list.  Start by computing forces.
	// (Assume exact = false and ignore_internal = false.)

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    predict_loworder_all(bi, sys_t);
	    bi->clear_interaction();
	    bi->top_level_node_prologue_for_force_calculation(false);
	}

	bool restart_grape = false;

//	b->get_kira_calculate_top_level_acc_and_jerk()(next_nodes, n_next,
//						       sys_t, restart_grape);

	kira_calculate_top_level_acc_and_jerk(next_nodes, n_next,
					      sys_t, restart_grape);

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    bi->inc_direct_force(n_top-1);
	    bi->top_level_node_epilogue_force_calculation();
	    bi->inc_steps();
	}

	// Complete calculation of accs and jerks by correcting for C.M.
	// interactions and applying external fields.

	kira_options *ko = b->get_kira_options();

	for (int i = 0; i < n_next; i++)
	    next_nodes[i]->set_on_integration_list();

	if (ko->use_old_correct_acc_and_jerk || !ko->use_perturbed_list) {
	    bool reset = false;
	    correct_acc_and_jerk(b, reset);		// old version
	} else
	    correct_acc_and_jerk(next_nodes, n_next);	// new version

	for (int i = 0; i < n_next; i++)
	    next_nodes[i]->clear_on_integration_list();

	if (b->get_external_field() > 0) {

	    // Add external forces.

	    for (int i = 0; i < n_next; i++) {
		hdyn *bi = next_nodes[i];
		real pot = bi->get_pot();
		vec acc = bi->get_acc(), jerk = bi->get_jerk();

		// (Really only need acc to be correctly set...)

		get_external_acc(bi, bi->get_pred_pos(), bi->get_pred_vel(),
				 pot, acc, jerk);
		bi->set_pot(pot);
		bi->set_acc(acc);
		bi->set_jerk(jerk);
	    }
	}

	// Apply corrector and redetermine timesteps.

	real st = sys_t;

	for (int i = 0; i < n_next; i++) {
	    hdyn *bi = next_nodes[i];
	    bi->correct_and_update();
	    bi->init_pred();
	    bi->store_old_force();

	    // As in synchronize_tree, make sure time step is consistent with
	    // system_time (= time).

	    real timestep = bi->get_timestep();

	    int iter = 0;
	    while (fmod(st, timestep) != 0) {
		if (iter++ > 30) break;
		timestep *= 0.5;
	    }

	    if (iter > 20) {
		cerr << "kira_synchronize_tree: " << bi->format_label() << " ";
		PRL(iter);
		PRI(4); PRC(fmod(st, timestep)); PRL(timestep);
	    }

	    bi->set_timestep(timestep);
	}

	delete [] next_nodes;

	// Top-level nodes are all synchronized.  Deal with low-level nodes.  

	if (sync_low_level) {
	    for_all_daughters(hdyn, b, bi) {
		hdyn *od = bi->get_oldest_daughter();
		if (od && !od->get_kepler()) {

		    // Synchronizing od will also synchronize its sister,
		    // but call synchronize_tree explicitly to ensure that
		    // substructure is properly handled.

		    synchronize_tree(od);

		    hdyn *yd = od->get_younger_sister();
		    if (yd) synchronize_tree(yd);		// else error?

		}
	    }
	}

	cerr << endl
	     << "end of synchronization"
	     << endl << flush;

    } else {

#if 0
	cerr << endl
	     << "synchronizing tree without GRAPE at time "
	     << b->get_system_time() << endl;
#endif

	synchronize_tree(b);
    }
}



// initialize_system_phase2:  Calculate acc, jerk, timestep, etc for all nodes.
//
// NOTE: System should be synchronized prior to calling this function.

// Static data:

static int work_size = 0;
static hdyn ** nodes = NULL;
static int nnodes = -1;

// Allow possibility of cleaning up if necessary:

void clean_up_kira_ev() {if (nodes) delete [] nodes;}

void initialize_system_phase2(hdyn *b,
			      int call_id,	// default = 0
			      int set_dt)	// 0 ==> set only if zero
						// 1 ==> set with limit (def)
						// 2 ==> always set
{
    // cerr << "initialize_system_phase2: "; PRC(call_id), PRL(set_dt);
    dbg_message("initialize_system_phase2", b);

    if (!b->is_root())
	err_exit("initialize_system_phase2 called with non-root");

    xreal time = b->get_system_time();

    int n = 0;
    for_all_nodes(hdyn, b, bb) n++;

    if  (work_size < n) {
	if (!nodes) delete [] nodes;
	work_size = n + 10;
	nodes = new hdynptr[work_size];
    }

    for_all_nodes(hdyn, b, bb) {
	if (!bb->get_kepler() && (bb->get_unperturbed_timestep() > 0)) {

	    // When is this necessary? (SLWM 3/98)

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl
		 << "initialize_system_phase2: "
		 << "creating kepler for unperturbed binary"
		 << endl
		 << "    " << bb->get_parent()->format_label()
		 << " at system time " << time
		 << "  call_id = " << call_id
		 << endl;
	    cerr.precision(p);

	    bb->update_kepler_from_hdyn();
	}
    }

    // Make a list of all nodes except unperturbed binary components.

    n = 0;
    for_all_nodes(hdyn, b, bb) {
	if ((bb != b) && (!bb->get_kepler())) {
	    nodes[n] = bb;
	    n++ ;
	}
    }

    bool tree_changed = true;
    bool reset_force_correction = true;	// no longer used
    bool restart_grape = true;
    bool exact = false;

    calculate_acc_and_jerk_for_list(nodes, n, time,
				    exact,
				    tree_changed,
				    reset_force_correction,  // no longer used
				    restart_grape);

    // (All perturber lists have been redetermined...)

    real min_dt = VERY_LARGE_NUMBER;

    for_all_nodes(hdyn, b, bb) {
	if ((bb != b) && (!bb->get_kepler())) {

	    bb->store_old_force();

	    if (set_dt || ((real)bb->get_time() <= 0
			   || bb->get_timestep() <= 0) ) {

		real dtlim = bb->get_timestep()/2;
		bb->set_first_timestep();
		if (set_dt == 1 && bb->get_timestep() < dtlim)
		    bb->set_timestep(dtlim);

		if (bb->get_time() + bb->get_timestep() < time) {

		    // Note from Steve and Simon, Feb 26, 1999:

		    // This bug can only occur if the function is
		    // improperly called.

		    cerr << endl << "warning: initialize_system_phase2: "
			 << "time will go backwards!" << endl;

		    int p = cerr.precision(HIGH_PRECISION);
		    PRL(bb->format_label());
		    PRC(time); PRL(bb->get_time() + bb->get_timestep());
		    cerr.precision(p);

		    cerr << "function should not be called with "
			 << "unsynchronized system." << endl << endl;
		    pp3(bb);

		    // Fix would be to force timestep to go past system
		    // time, but this is not in general possible while
		    // maintaining the block step structure.  Could
		    // simply terminate here, but for now we let the
		    // code die in kira.C.
		}
	    }

	    min_dt = Starlab::min(min_dt, bb->get_timestep());
	}
    }

    // Check for perturbed binaries in the input data...

    if (set_dt && (real)time <= 0) {
	for_all_nodes(hdyn, b, bb) {
	    if ((bb != b) && (!bb->get_kepler())) {
		if (bb->is_parent()
		    && bb->get_oldest_daughter()
			 ->get_perturbation_squared() > 1)
		    bb->set_timestep(min_dt);
	    }
	}
    }

    // cerr << "leaving initialize_system_phase2()\n\n";
}
