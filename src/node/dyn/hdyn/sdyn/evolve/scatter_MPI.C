//// modified scatter program for multiple experiments
//// scatter:  simple N-body scattering program.
////
//// Options:      -A    specify accuracy parameter [0.02]
////               -c    specify CPU time check interval (s) [3600]
////               -C    specify cube size for snapshot output [10]
////               -d    specify log output interval [none]
////               -D    specify snapshot output interval [none]
////               -g    specify tidal tolerance [1.E-6]
////               -i    specify mkscat-style string to initialize
////                        the integration ["-M 1 -v 2 -t -p"]
////               -n    specify the number of scattering experiments [1]
////               -p    read initial state from cin [false]
////               -t    specify time span of integration [infty]
////               -v    verbose mode [false]
////               -s    seed [system clock]
////                     Use negative value to prevent mkscat from
////                     reinitializing the random number sequence.

#ifdef USE_MPI
#include "mpi++.h"
#else
#include "localmpi++.h"
#endif
//#include "sdyn.h"
#include "scatter.h"
#include "sigma.h"




#ifndef TOOLBOX


ostream& operator<<(ostream& s, scatter_input& i) {
 
  s << "Scatter_input: " << endl;
  s << "Input string: " << i.init_string << endl;
  s << " i.dt = " << i.delta_t << " " << i.dt_out << " " << i.dt_snap << endl;
  s << " eta= " << i.eta << endl;
  s << " tol= " << i.tidal_tol_factor  << endl;
  s << " Cs= " << i.snap_cube_size << " dt cpu= " << i.cpu_time_check << endl;
  s << " e_exp= " << i.n_experiments << endl;
  s << " rnd= " << i.seed << " "<< i.n_rand <<" "<<i.n_rand_inc<<" "
    << i.pipe  << " "<< i.debug;

  return s;

}

void ppn(sdyn* b, ostream & s, int level) {

  s.precision(3);

  for (int i = 0; i < 2*level; i++) s << " ";

  if(b != b->get_root()) {
    b->pretty_print_node(s);
    cerr << endl;
    s << "   "<< b->get_mass() << "  "
      << b->get_pos() << " (r= " << abs(b->get_pos()) << ")  " 
      << b->get_vel() << " (v= " << abs(b->get_vel()) << ")" << endl;
    //	<< "r= " << abs(b->get_pos()) << "  " 
    //	<< "v= " << abs(b->get_vel()) << endl;
  }

  for (sdyn * daughter = b->get_oldest_daughter();
       daughter != NULL;
       daughter = daughter->get_younger_sister())
    ppn(daughter, s, level + 1);	
}


#define DT_CHECK	20

// scatter: Take the system with root node b and integrate it forward
// 	    in time up to time delta_t.  Check the state of the system
//	    every DT_CHECK time units.

void scatter(sdyn* b, scatter_input input, 
	     scatter_exp &experiment) {

  scatter(b, input.eta,
	     input.delta_t, input.dt_out, input.cpu_time_check,
	     input.dt_snap, input.tidal_tol_factor, input.snap_cube_size,
	     input.debug,
	     experiment);

}

void scatter(sdyn* b, real eta,
	     real delta_t, real dt_out, real cpu_time_check,
	     real dt_snap, real ttf, real snap_cube_size,
	     int debug,
	     scatter_exp &experiment) {

  experiment.init_scatter_exp(b); 

  real de_merge = 0;
  int stop_at_cpu = 0;
  real cpu_save = cpu_time();
  bool terminate = false;

  real t_out = b->get_time_offset();
  real t_end = delta_t + b->get_time_offset();

  char previous_form[255];
  char current_form[255];
  strcpy(current_form, get_normal_form(b));

  real kin=0, pot=0;
  b->flatten_node();
  real etot_init = calculate_energy_from_scratch(b, kin, pot); 
  make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);

  // loop continues for ever or unit termination time is reached or
  // system is unbound (see extend_or_end.C)
  for (real t = 0; t <= VERY_LARGE_NUMBER; t += DT_CHECK) {

    terminate = tree_evolve(b,  DT_CHECK, dt_out, dt_snap, 
			    snap_cube_size, 
			    eta, cpu_time_check);

    calculate_energy(b, kin, pot);
    if (debug) {
      cerr.precision(6), cerr << "\nStatus at time " << b->get_time();
      cerr.precision(9), cerr << " (energy = " << kin + pot << "):\n";
      cerr.precision(6);
    }

    
    // Check to see if the scattering is over.
    
    if (cpu_time_check < 0
	&& cpu_time() - cpu_save > abs(cpu_time_check)) {
      cerr <<" tcpu>tcpu_check"<<endl;
      return;
    }
    
    if(t>=dt_out) {

      real ekin, epot;
      calculate_energy(b, ekin, epot);
      if(debug==1) {
	cerr << "Time = " << b->get_time() 
	     << "  Etot = " << ekin + epot 
	     << "  Tinf = " << ekin/epot << endl;
      }
      put_sdyn(cout, *b);
      t_out += dt_out;
    }
    
    if(t>=t_end) {
      cerr << "Early termination"<<endl;

      experiment.set_scatter_discriptor(stopped);
      experiment.set_stop(true);
      real ekin, epot;
      calculate_energy(b, ekin, epot);
      cerr << "Time = " << b->get_time() 
	   << "  Etot = " << ekin + epot << endl;
      break;
    }
    
    int coll_flag = 2;
    de_merge += merge_collisions(b, coll_flag);
    
    make_tree(b, DYNAMICS, STABILITY, K_MAX, false);
    
    int unbound = extend_or_end_scatter(b, ttf, false);
    if(unbound==2) {   // two body system is bound but we stop anyway
      experiment.set_final_bound(true);
      unbound = 1;   // but stop anyway
    }
    //bool unbound = tree_is_unbound(b, ttf, false);
    //make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);
    if (unbound || terminate) break;
    
    //	if(experiment.get_form_changes()>1) 
    //	  cerr << "Resonant encounter: " << get_normal_form(b) << endl;
    strcpy(previous_form, current_form);
    strcpy(current_form, get_normal_form(b));
    if(strcmp(previous_form, current_form))
      experiment.inc_form_changes();
    
    //	print_normal_form(b, cerr);
  }
  
  b->flatten_node();
  real etot_error = calculate_energy_from_scratch(b, kin, pot) 
    - etot_init - de_merge;
  experiment.set_energy_error(etot_error);
  //PRL(etot_error);

  make_tree(b, DYNAMICS, STABILITY, K_MAX, debug);
  
  // add final report
  experiment.final_scatter_exp(b);
  //PRL(experiment.get_final_form());

}


#else

void slave_part_of_experiment(scatter_input input,
			      scatter_exp &experiment) {

  int random_seed = srandinter(input.seed, input.n_rand);

  cerr << "Random seed = " << get_initial_seed()
       << "  n_rand = " << get_n_rand() << flush << "\n";
    
  sdyn *b = mkscat(input.init_string);
    
  b->flatten_node();
    
  real kin, pot;
  real etot_init = calculate_energy_from_scratch(b, kin, pot);
  cerr << "Time = " << 0 
       << "  Etot = " << kin + pot 
       << " (" << kin << ", " << pot
       << ")  Tinf = " << kin/pot << endl;
  make_tree(b, !DYNAMICS, STABILITY, K_MAX, input.debug);
  
  // Initial output:
  b->set_name("root");
  
#if 0  
  cerr << "*** Initial configuration (random seed = "
       << get_initial_seed() << "):\n";
  print_normal_form(b, cerr);
    
  vector center = b->get_pos();
  print_structure_recursive(b, 0., center, true, true, 4);
#endif
    
  // Integrate the system to completion:
  scatter(b, input, experiment);

  cerr.precision(6);
  if (input.debug) cerr << endl;

#if 0
  cerr << "*** Final system configuration (time = " 
       << b->get_time() << "):\n";
  print_normal_form(b, cerr);

  ppn(b, cerr);
  //	vector center = b->get_pos();
  print_structure_recursive(b, 0., center, true, true, 4);

  cerr << "Energy error = " << experiment.get_energy_error() << endl;

  cerr << "Final Normal form:  ";
  print_normal_form(b, cerr);
#endif

  delete b;

}

void slave_process(int master, 
		   scatter_input &input, MPI_Datatype inputtype,
		   scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  int myid = MPI::COMM_WORLD.Get_rank();

  int length = 1;

  int request_for_data = 1;
  int data_available;
  int sending_data;

    do {

      request_for_data = myid;
      MPI::COMM_WORLD.Send(&request_for_data, length, MPI::INT, master, 
			   DATA_REQUEST_TAG);
      MPI::COMM_WORLD.Recv(&data_available, length, MPI::INT, master, 
			   DATA_REQUEST_TAG); 

      if(data_available<=0) {
	cerr << "Terminate slave "<< myid<<endl;
	return;
      }

      MPI::COMM_WORLD.Recv(&input, length, inputtype, master, 
			   DATA_COMING_TAG);
      MPI::COMM_WORLD.Recv(&experiment, length, scatter_exp_type, master, 
      			   DATA_COMING_TAG);

      slave_part_of_experiment(input, experiment);
	
      cerr << "Send by slave "<<flush;
      PRL(experiment.get_final_form());
      
      sending_data = myid;
      MPI::COMM_WORLD.Send(&sending_data, length, MPI::INT, master, 
			   DATA_COMING_TAG);
      MPI::COMM_WORLD.Send(&experiment, length, scatter_exp_type, 
			     master, DATA_COMING_TAG);
    }
    while(true);

}

void master_process(scatter_input &input, MPI_Datatype inputtype,
		    scatter_exp &experiment, MPI_Datatype scatter_exp_type) {

  int myid = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();

  int accept_all_senders = MPI::ANY_SOURCE;

  // Initialize random number generator
  int random_seed = srandinter(input.seed, input.n_rand);
  input.n_rand_inc = get_n_rand() - input.n_rand; 
  
  // initialize scatter_hist
  sdyn *b = mkscat(input.init_string);
  scatter_hist *hi = initialize_scatter_hist(b);
  delete b;
  b = NULL;

  input.n_rand_inc = get_n_rand() - input.n_rand; 

  int length = 1;
  int nsend = 0;
  int i;

  int n_exp = 0;

  int slave;
  int request_for_data;
  int data_available = 1;
  int slave_sends_data;

  do {

    MPI::COMM_WORLD.Recv(&request_for_data, length, MPI::INT,
    			 accept_all_senders, DATA_REQUEST_TAG);
    slave = request_for_data;
    MPI::COMM_WORLD.Send(&data_available, length, MPI::INT, slave, 
			 DATA_REQUEST_TAG);

    if (input.n_experiments > 1) cerr << "Experiment #" << n_exp
				      << ": " << "\n";

    input.n_rand += input.n_rand_inc;
    PRC(input.n_rand);PRL(input.n_rand_inc);

    MPI::COMM_WORLD.Send(&input, length, inputtype, slave, 
			 DATA_COMING_TAG);
    MPI::COMM_WORLD.Send(&experiment, length, scatter_exp_type, slave, 
			 DATA_COMING_TAG);
    nsend++;

    MPI::COMM_WORLD.Recv(&slave_sends_data, length, MPI::INT,
    			 accept_all_senders,  DATA_COMING_TAG);
    slave = slave_sends_data;
    MPI::COMM_WORLD.Recv(&experiment, length, scatter_exp_type, 
    			 slave,  DATA_COMING_TAG);

    //    PRL(experiment.get_final_form());
    //    cerr << experiment << endl;
    hi->add_scatter_hist(experiment, 0);
    //cerr << "scatter history" << flush << endl;
    //hi->put_scatter_hist(cerr, false);
    //cerr << " N done = " << n_exp+1 << endl;

    n_exp++;
  }
  while(n_exp < input.n_experiments);

  cerr << "Start terminating all processes"<<endl;

  data_available = 0;
  for(i=1; i<nproc; i++) {

      nsend++;

      MPI::COMM_WORLD.Recv(&request_for_data, length, MPI::INT,
			   accept_all_senders, DATA_REQUEST_TAG);
    
      slave = request_for_data;
      cerr << "Send termination signal to slave " << slave << endl;
      MPI::COMM_WORLD.Send(&data_available, length, MPI::INT, slave, 
			   DATA_REQUEST_TAG);
  }

  cerr << "scatter history" << flush << endl;
  hi->put_scatter_hist(cerr, false);
  cerr << "Master ends"<<endl;
}

void execute_scatter_experiment(scatter_input input) {

  int argc = 0;
  char ** argv = NULL;
  MPI::Init(argc, argv);
  int myid = MPI::COMM_WORLD.Get_rank();
  int nproc = MPI::COMM_WORLD.Get_size();

  // MPI definition of datatype input_type
  int inputblockcounts[3] = {255, 7, 6};
  MPI_Datatype ntypes[3];
  MPI::Aint displs[3];
  MPI_Datatype inputtype;

  MPI_Address(&input.init_string, &displs[0]);
  MPI_Address(&input.delta_t, &displs[1]);
  MPI_Address(&input.n_experiments, &displs[2]);
  ntypes[0] = MPI::CHAR;
  ntypes[1] = MPI::DOUBLE;
  ntypes[2] = MPI::INT;
  displs[1] -= displs[0];
  displs[2] -= displs[0];
  displs[0] = 0;

  MPI_Type_struct(3, inputblockcounts, displs, ntypes, &inputtype);
  MPI_Type_commit(&inputtype);

  scatter_exp experiment; 

  // MPI definition of datatype scatter_exp
  int charlength = 255 + 255;
  int reallength = 4;
  int intlength = 10 + 2*N_RHO_ZONE_MAX + N_STEP_BIN + N_OSC_BIN;
  int boolength = 3 + intlength;
  int blockcounts[3] = {charlength, reallength, boolength};
  MPI::Aint exp_displs[3];
  MPI_Datatype scatter_exp_type;

  MPI_Address(experiment.get_initial_form(), &exp_displs[0]);
  MPI_Address(experiment.get_time_ptr(), &exp_displs[1]);
  MPI_Address(experiment.get_scatter_discriptor_ptr(), &exp_displs[2]);
  MPI_Datatype nexp_types[3] = {MPI::CHAR, MPI::DOUBLE, MPI::INT};
  exp_displs[1] -= exp_displs[0];

  exp_displs[2] -= exp_displs[0];
  exp_displs[0] = 0;

  MPI_Type_struct(3, blockcounts, exp_displs, nexp_types, &scatter_exp_type);
  MPI_Type_commit(&scatter_exp_type);

  int master = 0;
  if(myid==master) {
    master_process(input, inputtype, 
		   experiment, scatter_exp_type);
  }
  else {
    slave_process(master, 
		  input, inputtype, 
		  experiment, scatter_exp_type);
  }

  MPI::Finalize();

}

main(int argc, char **argv) {

  scatter_input input;
  char* default_init  
    = "-M 0.879 -rm 3 -v 2 -r 1 -t -r1 0.0508 -r2 0.0348 -e 0 -q 0.567 -p -a 1 -q 1 -r1 0.0394 -r2 0.0394";        // Iota Ori Probleem 

  strcpy(&input.init_string[0], default_init);
//  input.init_string = default_init;	

  bool D_flag = false;
  check_help();

  extern char *poptarg;
  int c;
  char* param_string = "A:c:C:d:D:g:i:n:N:ps:t:v";

  while ((c = pgetopt(argc, argv, param_string)) != -1)
    switch(c) {

    case 'A': input.eta = atof(poptarg);
      break;
    case 'c': input.cpu_time_check = atof(poptarg);
      break;
    case 'C': input.snap_cube_size = atof(poptarg);
      break;
    case 'd': input.dt_out = atof(poptarg);
      break;
    case 'D': D_flag = TRUE;
      input.dt_snap = atof(poptarg);
      break;
    case 'g': input.tidal_tol_factor = atof(poptarg);
      break;
    case 'i': strcpy(&input.init_string[0], poptarg);
      break;
    case 'n': input.n_experiments = atoi(poptarg);
      break;   
    case 'N': input.n_rand = atoi(poptarg);
      break;
    case 'p': input.pipe = 1 - input.pipe;
      break;
    case 't': input.delta_t = atof(poptarg);
      break;
    case 'v': input.debug = 1 - input.debug;
      break;
    case 's': input.seed = atoi(poptarg);
      break;
    case '?': params_to_usage(cerr, argv[0], param_string);
      get_help();
    }            

  cerr << input << endl;

  if (!D_flag) input.dt_snap = input.delta_t; // Guarantee output at end

  execute_scatter_experiment(input);

}
#endif








