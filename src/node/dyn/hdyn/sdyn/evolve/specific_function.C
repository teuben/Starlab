#include "sigma.h"

#ifndef TOOLBOX

#ifdef USE_specific_function
#ifdef USE_MPI
void print_scatter_specific_information(sdyn *b,
					sigma_input input, 
					scatter_exp experiment) {

  char *prefered_string = "(p1,p2,(t1,t2))";
  //  if(experiment.get_scatter_discriptor()==single_ionization ||
  //     experiment.get_scatter_discriptor()==exchange_ionization) {
  if(!strcmp(experiment.get_final_form(), "((p1,t1),p2,t2)") ||
     !strcmp(experiment.get_final_form(), "(p1,(p2,t1),t2)") ||
     !strcmp(experiment.get_final_form(), "(p1,p2,(t1,t2))")) {
    //  if (true) {
    int myid = MPI::COMM_WORLD.Get_rank();
    MPI::Status status;  
    int request_for_data = myid;
    int request_granted = 1;
    int master = 0;
    int length = 1;
    //    cerr << "process " << myid << " requests for output"<<endl;
    MPI::COMM_WORLD.Send(&request_for_data, length, MPI::INT, master, 
			 WRITE_REQUEST_TAG);
    MPI::COMM_WORLD.Recv(&request_granted, length, MPI::INT, master, 
			 WRITE_REQUEST_TAG, status); 

    //    cerr << "slave " << myid << " is allowed to write" << endl;

    cerr << "Found prefered final state: "
	 << experiment.get_final_form() << endl;
    ppn(b, cerr);
    vector center = b->get_pos();
    print_structure_recursive(b, 0., center, true, true, 4);

    MPI::COMM_WORLD.Send(&request_for_data, length, MPI::INT,
			 master, WRITE_READY_TAG);
    //    cerr << "process " << myid << " finished output"<<endl;

  }

}

#else // No MPI
void print_scatter_specific_information(sdyn *b,
					sigma_input input, 
					scatter_exp experiment) {

  char *prefered_string = "(p1,p2,(t1,t2))";
  //  if(experiment.get_scatter_discriptor()==single_ionization ||
  //     experiment.get_scatter_discriptor()==exchange_ionization) {
  if(!strcmp(experiment.get_final_form(), "((p1,t1),p2,t2)") ||
     !strcmp(experiment.get_final_form(), "(p1,(p2,t1),t2)") ||
     !strcmp(experiment.get_final_form(), "(p1,p2,(t1,t2))")) {

    cerr << "Found prefered final state: "
	 << experiment.get_final_form() << endl;
    ppn(b, cerr);
    vector center = b->get_pos();
    print_structure_recursive(b, 0., center, true, true, 4);

  }
}
#endif

#else
void print_scatter_specific_information(sdyn *b,
					sigma_input input, 
					scatter_exp experiment) {
// do nothing
}
#endif
#endif
