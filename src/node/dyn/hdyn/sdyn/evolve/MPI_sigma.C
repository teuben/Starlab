////     MPI_sigma     toolbox with routines which call MPI
////
////
//
//                     Simon Portegies Zwart   MIT  Nov 2000

#include "sigma_MPI.h"

#ifndef TOOLBOX

#ifdef HAS_MPI

void initialize_MPI(int &myid, int &nproc) {

  int argc = 0;
  char **argv = NULL;
  MPI::Init(argc, argv);

  myid = MPI::COMM_WORLD.Get_rank();
  nproc = MPI::COMM_WORLD.Get_size();
}

void finalize_MPI() {

  MPI::Finalize();
}

MPI_Datatype sigma_input::initialize_data_structures_MPI() {

  // MPI definition of datatype input_type
  int inputblockcounts[4] = {255, 7, 7, 7};
  MPI::Aint displs[4];  
  MPI_Datatype inputtype;

  MPI_Address(&init_string, &displs[0]);
  MPI_Address(&delta_t, &displs[1]);
  MPI_Address(&n_experiments, &displs[2]);
  MPI_Address(&rho_sq_min, &displs[3]);
  MPI_Datatype ntypes[4] = {MPI::CHAR, MPI::DOUBLE, MPI::INT, MPI::DOUBLE};
  displs[1] -= displs[0];
  displs[2] -= displs[0];
  displs[3] -= displs[0];
  displs[0] = 0;

  MPI_Type_struct(4, inputblockcounts, displs, ntypes, &inputtype);
  MPI_Type_commit(&inputtype);

  return inputtype;
}

#else // No MPI

void initialize_MPI(int &myid, int &nproc) { 
  // do nothing
}

void finalize_MPI() {
  // do nothin
}

#endif

void execute_sigma_experiment(sigma_input &input) {

  int myid=0, nproc;
  initialize_MPI(myid, nproc);

  scatter_exp experiment; 

  MPI_Datatype inputtype = input.initialize_data_structures_MPI();
  MPI_Datatype scatter_exp_type = experiment.initialize_data_structures_MPI();

  int master = 0;
  if(myid==master) {
    get_sigma(input, inputtype, experiment, scatter_exp_type);
  }
  else {
    slave_process(master, 
		  input, inputtype, 
		  experiment, scatter_exp_type);
  }

  finalize_MPI();

}

local void slave_part_of_experiment(sigma_input input,
				    scatter_exp &experiment) {

  int random_seed = srandinter(input.seed, input.n_rand);

  if(input.verbose==1) 
    cerr << "Random seed = " << get_initial_seed()
         << "  n_rand = " << get_n_rand() << flush << "\n";
    
  sdyn *b = mkscat(input.init_string, input);
    
  b->flatten_node();
    
  real kin, pot;
  real etot_init = calculate_energy_from_scratch(b, kin, pot);
  if(input.verbose==1) 
    cerr << "Time = " << 0 << "  Etot = " << kin + pot 
         << " (" << kin << ", " << pot
         << ")  Tinf = " << kin/pot << endl;
  make_tree(b, !DYNAMICS, STABILITY, K_MAX, input.debug);
  b->set_name("root");
  
  // Initial output:
    
  if(input.verbose==2) {
    cerr << "*** Initial configuration (random seed = "
         << get_initial_seed() << "):\n";

    print_normal_form(b, cerr);
    vector center = b->get_pos();
    print_structure_recursive(b, 0., center, true, true, 4);
  }
    
  // Integrate the system to completion:
  int hit = single_scatter(b, (scatter_input)input, experiment);

  //  PRC(hit);
  //  PRC(experiment.get_nzone());
  //  PRL(experiment.get_nhits(experiment.get_nzone()));

  cerr.precision(6);
  if (input.debug) cerr << endl;

  if(input.verbose==2) {
    cerr << "*** Final system configuration (time = " 
	 << b->get_time() << "):\n";
    cerr << "Final Normal form:  ";
    print_normal_form(b, cerr);

    ppn(b, cerr);
    vector center = b->get_pos();
    print_structure_recursive(b, 0., center, true, true, 4);

    cerr << "Energy error = " << experiment.get_energy_error() << endl;
  }

  print_scatter_specific_information(b, input, experiment);

  delete b;

}

#ifdef HAS_MPI
void slave_process(int master, 
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  int myid = MPI::COMM_WORLD.Get_rank();
  MPI::Status status;  

  int length = 1;

  int request_for_data = 1;
  int data_available;
  int sending_data;

    do {

      request_for_data = myid;
      //      cerr << "Slave " << myid << " requests new job "<<endl;

      MPI::COMM_WORLD.Send(&request_for_data, length, MPI::INT, master, 
			   DATA_REQUEST_TAG);
      MPI::COMM_WORLD.Recv(&data_available, length, MPI::INT, master, 
			   DATA_REQUEST_TAG, status); 
      //      PRC(myid);PRL(data_available);
      if(data_available == SEND_HOLD_TAG) {
	//	cerr << "Slave " << myid << " is put on hold until furthre notice"<<endl;
	MPI::COMM_WORLD.Recv(&data_available, length, MPI::INT, master, 
			     STOP_HOLD_TAG, status); 
	//	cerr << "Slave " << myid << " continues his work"<<endl;
	MPI::COMM_WORLD.Send(&request_for_data, length, MPI::INT, master, 
			     DATA_REQUEST_TAG);
	MPI::COMM_WORLD.Recv(&data_available, length, MPI::INT, master, 
			     DATA_REQUEST_TAG, status); 
      }

      if(data_available<=0) {
	cerr << "Terminate slave "<< myid<<endl;
	return;
      }

      //      cerr << "Slave " << myid << " waiting for data "<<endl;

      MPI::COMM_WORLD.Recv(&input, length, inputtype, master, 
			   DATA_COMING_TAG);
      MPI::COMM_WORLD.Recv(&experiment, length, scatter_exp_type, master, 
      			   DATA_COMING_TAG);

      //      cerr << "Slave " << myid << " received data (now working) "<<endl;
      slave_part_of_experiment(input, experiment);
	
      //      cerr << "Slave " << myid << " requests to send data "<<endl;
      sending_data = myid;
      MPI::COMM_WORLD.Send(&sending_data, length, MPI::INT, master, 
			   DATA_COMING_TAG);
      MPI::COMM_WORLD.Send(&experiment, length, scatter_exp_type, 
			     master, DATA_COMING_TAG);
    }
    while(true);

}
#else
void slave_process(int master, 
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

      slave_part_of_experiment(input, experiment);
}
#endif

#ifdef HAS_MPI
void terminate_all_processors() {

  cerr << "Start terminating all processes"<<endl;

  int myid = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();

  int accept_all_senders = MPI::ANY_SOURCE;

  int slave;
  int request_for_data;
  int length = 1;
  int data_available = 0;
  for(int i=1; i<nproc; i++) {

      MPI::COMM_WORLD.Recv(&request_for_data, length, MPI::INT,
			   accept_all_senders, DATA_REQUEST_TAG);
    
      slave = request_for_data;
      MPI::COMM_WORLD.Send(&data_available, length, MPI::INT, slave, 
			   DATA_REQUEST_TAG);
  }

  cerr << "All processes terminated" << endl;
}
#else
void terminate_all_processors() {
  // dummy
}
#endif

#ifdef HAS_MPI
int master_process(sigma_out & out,
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  int myid = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();
  MPI::Status status;  

  // keep track of holding processes
  bool *hold      = new bool[nproc];
  real *hold_time = new real[nproc];
  for(int i=0; i<nproc; i++) {
    hold[i] = false;
    hold_time[i] = 0;
  }

  int accept_all_senders = MPI::ANY_SOURCE;
  int accept_all_tags    = MPI::ANY_TAG;

  int length = 1;
  int nsend = 0;

  int n_exp = 0;

  int slave, tag;
  int request_for_data;
  int data_available = 1;
  int slave_sends_data;

  do {

    real cpu_init = cpu_time();
    MPI::COMM_WORLD.Recv(&request_for_data, length, MPI::INT,
    			 accept_all_senders, accept_all_tags, status);

    slave = status.Get_source();
    tag = status.Get_tag();

    //    cerr << "Master received request from slave " << slave 
    //	 << " for " << tag << endl;

    switch(tag) {
    case DATA_REQUEST_TAG:
      if(nsend>=input.n_experiments) {
	//	cerr << "Put process " << slave << " on hold" << endl;
	hold[slave] = true;
	hold_time[slave] = cpu_time();
	int please_hold = SEND_HOLD_TAG;
	MPI::COMM_WORLD.Send(&please_hold, length, MPI::INT, 
			     slave, DATA_REQUEST_TAG);

      }
      else {
	//      cerr << "Master: Slave " << slave << " requests data"<<endl;
	MPI::COMM_WORLD.Send(&data_available, length, MPI::INT, slave, 
			     DATA_REQUEST_TAG);
	input.n_rand += input.n_rand_inc;
	
	MPI::COMM_WORLD.Send(&input, length, inputtype, slave, 
			     DATA_COMING_TAG);
	MPI::COMM_WORLD.Send(&experiment, length, scatter_exp_type, slave, 
			     DATA_COMING_TAG);
	nsend++;
      }
      break;
    case WRITE_REQUEST_TAG:
      //      cerr << "Master: Slave " << slave << " requests to write"<<endl;
      MPI::COMM_WORLD.Send(&data_available, length, MPI::INT, 
			   slave, WRITE_REQUEST_TAG);
      MPI::COMM_WORLD.Recv(&request_for_data, length, MPI::INT,
			   slave, WRITE_READY_TAG, status);
      break;
    case DATA_COMING_TAG:
      //      cerr << "Master: Slave " << slave << " requests send data"<<endl;
      MPI::COMM_WORLD.Recv(&experiment, length, scatter_exp_type, 
			   slave,  DATA_COMING_TAG, status);
      single_scatter_stats(&experiment, out);
      cerr << "Experiment (" << out.n_zone << ") #" << n_exp << ": ";
      cerr << " cpu= " << cpu_time() - cpu_init
	   << "   Time = " << experiment.get_time() 
	   << "   result: " << experiment.get_final_form() 
	   << " (" << experiment.get_scatter_discriptor() << ")"
	   << endl;
      n_exp++;
      break;
    };
    
    //    PRC(out.n_zone);PRC(nsend);PRC(n_exp);PRL(input.n_experiments);
  }
  while(n_exp < input.n_experiments);

  // resume holding processes.
  for(slave=1; slave<nproc; slave++) {
    if(hold[slave]) {
      //      cerr << "stop holding for process " << slave << endl;
	int dummy=1;
	MPI::COMM_WORLD.Send(&dummy, length, MPI::INT, 
			     slave, STOP_HOLD_TAG);
	hold[slave] = false;
	//	cerr << "Slave " << slave << " has hold for " 
	//	     << cpu_time() - hold_time[slave] << " seconds" << endl;

    }
  }
      
  cerr << "Master ("<<myid<<"): has done "<<n_exp << " experiments " <<endl;
  //  PRC(out.n_zone);PRC(n_hits);
  //  PRL(experiment.get_nhits(out.n_zone));

  return experiment.get_nhits(out.n_zone);

}
#else
int master_process(sigma_out & out,
		   sigma_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  input.n_rand_inc = get_n_rand() - input.n_rand; 
  int n_exp = 0;

  int master = 0;
  do {
    cerr << "Experiment (" << out.n_zone << ") #" << n_exp << ": ";

    input.n_rand += input.n_rand_inc;
    experiment.set_nzone(out.n_zone);

    real cpu_init = cpu_time();
    slave_process(master, 
		  input, inputtype, 
		  experiment, scatter_exp_type);

    cerr << " cpu= " << cpu_time() - cpu_init
	 << "   Time = " << experiment.get_time() 
	 << "   result: " << experiment.get_final_form() << endl;
    single_scatter_stats(&experiment, out);

    n_exp++;
  }
  while(n_exp < input.n_experiments);

  return experiment.get_nhits(out.n_zone);
}
#endif

#else

main(int argc, char **argv) {

  sigma_input input;

  // identical binary collision
  //char* default_init  
  //      = "-M 1 -rm 3 -v 1 -t -r1 0 -r2 0 -q 1 -p -a 1 -q 1 -r1 0 -r2 0";

    char* default_init  
      = "-M 0.66667 -rm 3 -S 30 -v 0 -t -r1 0 -r2 0 -e 0 -q 0.5 -p -a 1 -q 1 -r1 0 -r2 0";        // Iota Ori Probleem 


  // I-orionis problem with zero radii stars
  //char* default_init  
  // = "-M 0.5 -v 0 -r 1 -t -r1 0 -r2 0 -q 0.5 -p -a 1 -q 0.00001 -r1 0 -r2 0";        

  // Simplified I-orionis problem with zero radii stars
  //char* default_init  
  //  = "-M 0.66667 -v 0 -r 1 -t -r1 0 -r2 0 -q 0.50 -p -a 1 -q 1 -r1 0 -r2 0";

  // Simplified I-orionis problem with non-zero radii stars
  //  char* default_init  
  //      = "-M 0.879 -v 2 -r 1 -t -r1 0.0508 -r2 0.0348 -q 0.567 -p -a 1 -q 1 -r1 0.0394 -r2 0.0394";

  strcpy(&input.init_string[0], default_init);

  //    check_help();

    real  delta_t = VERY_LARGE_NUMBER;       // time span of the integration
    real  dt_out = VERY_LARGE_NUMBER;       // time output interval

    extern char *poptarg;
    int c;
    char* param_string = "A:c:C:d:D:e:g:Ii:m:M:N:pqQs:t:v:V:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'A': input.eta = atof(poptarg);
		      break;
	    case 'c': input.cpu_time_check = 3600*atof(poptarg);
	                                     // (Specify in hours)
		      break;
	    case 'd': input.max_trial_density = atof(poptarg);
		      break;
	    case 'D': input.dt_out = atof(poptarg);
	    case 'g': input.tidal_tol_factor = atof(poptarg);
		      break;
		      //case 'I': intermediate_sigma = 1 - intermediate_sigma;
		      //break;
	    case 'i': strcpy(input.init_string, poptarg);
		      break;
	    case 'M': input.pmass = atof(poptarg);
		      break;
	    case 'm': input.pmass = atof(poptarg);
		      break;
	    case 'N': input.n_rand = atoi(poptarg);
		      break;
		      // case 'p': print_counts = 1 - print_counts;
		      // break;
	    case 's': input.seed = atoi(poptarg);
		      break;
	    case 't': input.delta_t = atof(poptarg);
		      break;
	    case 'v': input.v_inf = atof(poptarg);
		      break;
   	    case 'V': //input.debug = atoi(poptarg);
	              input.verbose = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	      //		      get_help();
	}

    execute_sigma_experiment(input);
    
}

#endif
