
#ifdef HAS_MPI
#include <mpi++.h>
#else
#include "localmpi++.h"
#endif

#include "stdio.h"
#include "scatter_exp.h"

void scatter_exp::initialize_to_zero() {

  min_min_ssd = VERY_LARGE_NUMBER;
  min_min_time = 0;
  min_nn_dr2 = NULL;
  min_nn_label = NULL;
  min_nn_of_label = NULL;
        
  //  initial_form = new char[255];
  //  final_form = new char[255];

  int i;
  for(i=0; i< N_RHO_ZONE_MAX; i++) 
    n_hits[i] = n_hit[i] = 0;
  for(i=0; i<N_STEP_BIN; i++)
    step_counter[i] = 0;
  for(i=0; i<N_OSC_BIN; i++)
    osc_counter[i] = 0;

  n_coll = 0;
  n_star = 0;

  energy_error = 0;

  final_bound = false;
  resonance = false;
  stop = false;
  sd = not_specified; 
  n_zone = 0;
  form_changes = 0;
  n_initial = 0;
  n_final = 0;
  id_scenario = 0;
  n_found = 0;
  sigma = 0;
  sigma_err_sq = 0; 
}

real scatter_exp::get_specific_sigma_err_sq(scatter_discriptor d) {

    if(d == sd) 
      return sigma_err_sq;
    else 
      return 0;
}

real scatter_exp::get_specific_sigma(scatter_discriptor d) {

    if(d == sd) 
      return sigma;
    else 
      return 0;
}

int scatter_exp::get_specific_counts(scatter_discriptor d) {

    if(d == sd) 
      return n_found;
    else 
      return 0;
}

int scatter_exp::get_specific_counts(scatter_discriptor d, int i) {

    if(d == sd) 
      return n_hits[i];
    else 
      return 0;
}

bool scatter_exp::identical_to(scatter_exp *ha) {

  bool the_same = false;

  if (!strcmp(&final_form[0], ha->get_final_form()) && 
      resonance == ha->get_resonance() &&
      sd == ha->get_scatter_discriptor()) {

    the_same = true;
    set_nzone(ha->get_nzone());
  }

  return the_same;
}

bool scatter_exp::operator == (scatter_exp& ha) const {

  bool the_same = false;
  if(strcmp(&initial_form[0], ha.get_initial_form())==0 && 
     strcmp(&final_form[0], ha.get_final_form())==0 &&
     resonance == ha.get_resonance() &&
    sd == ha.get_scatter_discriptor()) {
      
    the_same = true;
  }

  return the_same;
}

bool scatter_exp::operator != (scatter_exp& ha) const {

  bool different = true;
  if (*this == ha)
    different = false;

  return different;
}
 
local char * type_string(scatter_discriptor sd) {
  
  switch(sd) {
  case preservation:                     return "pres";
       break;
  case exchange:                         return "exch";
       break;
  case exchange_ionization:              return "ex_ion";
       break;
  case single_ionization:                return "s_ion";
       break;
  case multiple_ionization:              return "m_ion";
       break;
  case total_ionization:                 return "t_ion";
       break;
  case collision_binary:                 return "bcoll";
       break;
  case collision_ionization:             return "c_ion";
       break;
  case stable_triple:                    return "trip";
       break;
  case stable_higher_order:              return "ntup";
       break;
  case two_body_collision:               return "2_coll";
       break;
  case multiple_body_collision:          return "n_coll";
       break;
  case unknown:                          return "unknown";
       break;
  case error:                            return "error";
       break;
  case stopped:                          return "stopped";
       break;
  case not_specified:                    return "NAS";
       break;
  case number_of_scatter_discriptors:    return "N";
       break;
  default:                               return "???";
  };
}

//char dummy_string[74];

ostream& operator<<(ostream& s, scatter_exp& hi) {

#if 0
  s << hi.id_scenario << " " << hi.n_found << " "
    << hi.sd << " " <<hi.resonance << " "
    << hi.get_final_form() << " \t";

  for(int i=0; i<N_RHO_ZONE_MAX; i++) 
    s << hi.get_nhits(i) << " ";
  s << endl;
#endif

  int i;
  s << hi.id_scenario;
  for(i=0; i<3-floor(log10((real)hi.id_scenario)+1); i++)
    s << " ";
  s << hi.n_found;
  for(i=0; i<2-floor(log10((real)hi.n_found+1)); i++)
    s << " ";
  s << type_string(hi.sd);
  for(i=0; i<8-strlen(type_string(hi.sd)); i++)
    s << " ";
  if(hi.resonance)
    s << "r  ";
  else
    s << "   ";
  int istrl = strlen(hi.get_initial_form());
  for(i=0; i<istrl-strlen(hi.get_final_form()); i++)
    s << " ";
  s << hi.get_final_form();
  //       for(i=0; i<N_RHO_ZONE_MAX; i++) 
  //	 s << hi.get_nhits(i);

  
  int n_space = 0;
  for(i=0; i<N_RHO_ZONE_MAX; i++) 
    n_space += (int)floor(log10((real)hi.get_nhits(i)+1));
  for(i=0; i<Starlab::max(0,N_RHO_ZONE_MAX-n_space-15); i++) 
    s << ".";

  for(i=0; i<N_RHO_ZONE_MAX; i++) {
    for(int j=0; j<1-floor(log10((real)hi.get_nhits(i)+1)); j++)
      s << " ";
    s << hi.get_nhits(i);
  }

  return s;

}

scatter_exp* read_scatter_exp(istream& s) {

  scatter_exp* hi = new scatter_exp();
    if(s.eof() || hi == NULL) 
      return NULL;
    
    int id;
    char iform[255];
    char fform[255];
    s >> id >> &iform[0] >> &fform[0]; 

    hi->set_id_scenario(id);
    hi->set_final_form(fform);
    hi->set_initial_form(iform);

    return hi;
}


void scatter_exp::init_scatter_exp(sdyn* b) {

  time = 0;
  form_changes = 0;
  resonance = false;
  sd = not_specified; 
  //  initial_form = new char[255];
  strcpy(&initial_form[0], get_normal_form(b));
  n_initial = b->n_leaves();
  //  final_form = NULL;

  // allocate space to keep track of nn's and their names
  min_nn_dr2 = new real[n_initial];
  min_nn_label = new int[n_initial];
  min_nn_of_label = new int[n_initial];
}

void scatter_exp::final_scatter_exp(sdyn* b) {

  if(form_changes>1)
    resonance = true;

  time = b->get_time();
  n_star = b->n_leaves();
  strcpy(&final_form[0], get_normal_form(b));
  // is the final system a binary?
  n_final = b->n_leaves();
  sd = classify_scatter();

  //cerr << "Update min_nn_dr2 for experiment in final_scatter_exp()" << endl;
  
  int i=0;
  for_all_leaves(sdyn, b, bi) {
    min_nn_label[i] = bi->get_min_nn_label();
    min_nn_dr2[i] = bi->get_min_nn_dr2();
    min_nn_of_label[i] = bi->get_index();

    i++;
  }
  
}


int scatter_exp::count_character_in_string(char string[], char search) {


  int counts = 0;
  for(int i=0; i<strlen(string); i++) {
    if(strncmp(&string[i], &search, 1)==0)
      counts++;
  }

  return counts;

}

int scatter_exp::count_character_in_string(char string[], char search[],
					   int n_char) {


  int counts = 0;
  for(int i=0; i<strlen(string); i++) {
    if(strncmp(&string[i], &search[0], n_char)==0)
      counts++;
  }

  return counts;

}

// returns 3 if triple, 4 if quadruple etc.
int scatter_exp::count_multiplicity() {

  char *final = get_final_form();
  int n_fstr = strlen(final);

  int b=0, bn=0;
  for(int i=0; i<n_fstr; i++) {
    if (strncmp(&final[i], "(", 1)==0)
      b++;
    else if(strncmp(&final[i], ")", 1)==0)
      b--;
    bn = Starlab::max(bn, b);
  }
  if(b!=0)
    cerr << "Not enough braces found in scatter_exp::count_mulitplicity()"
	 << endl;

  return bn;
}

// At the moment quite rudimentary
bool scatter_exp::check_for_exchange() {

  bool exchange = false;

  char *init = get_initial_form();
  char *final = get_final_form();
  int n_istr = strlen(init);
  int n_fstr = strlen(final);
  char inames[64];
  char fnames[64];

  int icnt=0;
  for(int i=0; i<n_istr; i++) {
    if (!(strncmp(&init[i], "(", 1)==0 ||
	  strncmp(&init[i], ")", 1)==0 ||
	  strncmp(&init[i], "+", 1)==0 ||
	  strncmp(&init[i], ",", 1)==0)) {
      sprintf(&inames[icnt], "%c", init[i]);
      icnt++;
    }
  }

  int fcnt=0;
  for(int f=0; f<n_fstr; f++) {
    if (!(strncmp(&final[f], "(", 1)==0 ||
	  strncmp(&final[f], ")", 1)==0 ||
	  strncmp(&final[f], "+", 1)==0 ||
	  strncmp(&final[f], ",", 1)==0)) {
      sprintf(&fnames[fcnt], "%c", final[f]);
      fcnt++;
    }
  }

  if(icnt!=fcnt) {
    cerr << "ERROR: icnt (" << icnt << ") != fcnt (" << fcnt << ")" << endl;
    cerr << "       in   bool scatter_exp::check_for_exchange()" << endl;
    icnt = Starlab::min(icnt, fcnt);
  }
  if(strncmp(inames, fnames, icnt)!=0)
    exchange = true;

  return exchange;
}


scatter_discriptor scatter_exp::classify_scatter() {

  scatter_discriptor discriptor = preservation;

  if(sd == stopped) {
    discriptor = stopped;
    stop = true;
  }
  else if(abs(energy_error)>1.e-5)
    discriptor = error; 
  else if(strcmp(get_initial_form(), get_final_form())==0) 
    discriptor = preservation;                           // ((p1,p2),(t1,t2))
  else if(count_character_in_string(get_final_form(), '+')>=1 &&
	  count_character_in_string(get_final_form(), '(')==1) {
    if (get_n_final()==2 && get_final_bound())
      discriptor = collision_binary;                    // (p1+p2,t1+t2)
    else
      discriptor = collision_ionization;                // (p1+p2,t1,t2)
  }
  else if(count_character_in_string(get_final_form(), '+')>=1 &&
	  count_character_in_string(get_final_form(), '(')>=2)
    discriptor = collision_binary;                      // (p1+p2,(t1,t2))
  else if(count_multiplicity()==3) 
    discriptor = stable_triple;                         // (p1,(p2,(t1,t2)))
  else if(count_multiplicity()>3) 
    discriptor = stable_higher_order;               // (p11,(p12,(p2,(t1,t2))))
  else if(count_character_in_string(get_final_form(), '(')==1)
    discriptor = total_ionization;                      // (p1,p2,t1,t2)
  else if(check_for_exchange()) {
    if (count_character_in_string(get_initial_form(), '(') >
	count_character_in_string(get_final_form(), '('))
      discriptor = exchange_ionization;                 // (p1,t1,(p2,t2))
    else
      discriptor = exchange;                            // ((p1,t1),(p2,t2))
  }
  else if(count_character_in_string(get_initial_form(), '(') >
	  count_character_in_string(get_final_form(), '(')+1)
    discriptor = multiple_ionization;                   // (p11,p12,p2,t1,t2)
  else if(count_character_in_string(get_initial_form(), '(') >
	  count_character_in_string(get_final_form(), '('))
    discriptor = single_ionization;                     // (p1,p2,(t1,t2))
  else
    discriptor = unknown;

  n_coll = count_character_in_string(get_final_form(), '+');

  return discriptor;
}

#ifdef HAS_MPI

MPI_Datatype scatter_exp::initialize_data_structures_MPI() {

  // MPI definition of datatype scatter_exp
  int charlength = 2*255;
  int reallength = 4;
  int intlength = 10 + 2*N_RHO_ZONE_MAX + N_STEP_BIN + N_OSC_BIN;
  int boolength = 3;
  int blockcounts[3] = {charlength, reallength, intlength+boolength};
  MPI::Aint exp_displs[3];
  MPI_Datatype scatter_exp_type;

  MPI_Address(get_initial_form(), &exp_displs[0]);
  MPI_Address(get_time_ptr(), &exp_displs[1]);
  MPI_Address(get_scatter_discriptor_ptr(), &exp_displs[2]);
  MPI_Datatype nexp_types[3] = {MPI::CHAR, MPI::DOUBLE, MPI::INT};
  exp_displs[1] -= exp_displs[0];
  exp_displs[2] -= exp_displs[0];
  exp_displs[0] = 0;

  MPI_Type_struct(3, blockcounts, exp_displs, nexp_types, &scatter_exp_type);
  MPI_Type_commit(&scatter_exp_type);

  return scatter_exp_type;
}

#else

MPI_Datatype scatter_exp::initialize_data_structures_MPI() {
  MPI_Datatype dummy;
  return dummy;
}

#endif

#ifdef TOOLBOX

void main(int argc, char **argv) {

    check_help();
    extern char *poptarg;
    int c;
    char* param_string = "c:";

    real cpu_time_check;
    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c': cpu_time_check = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
	}            

    do {
 
      //      scatter_exp *hi = read_scatter_exp(cin);
      scatter_exp *hi = new scatter_exp(cin);
      cerr << *hi << endl; 
      delete hi; 
    }
    while(!cin.eof());
}

#endif
