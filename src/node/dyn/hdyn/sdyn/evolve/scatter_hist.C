#include "scatter_hist.h"

scatter_hist* initialize_scatter_hist(sdyn* b) {

  if(b == NULL)
    return NULL;
  
  scatter_hist* hi = new scatter_hist();
  
  hi->init_scatter_exp(b);
  hi->final_scatter_exp(b);
  hi->set_id_scenario(1);
  hi->set_n_found(0);
  hi->set_resonance(false);
  
  return hi;
}

void scatter_hist::inc_id_scenario() {

  id_scenario = past->get_id_scenario()+1;

}

scatter_hist* scatter_hist::get_identical_scatter_hist() {

  scatter_hist *ha = NULL;
  for_all_scatter_hist(scatter_hist, this, ho) {
    if(this == ha) {
      ha = ho;
      break;
    }
  }

  return ha;
}

void scatter_hist::put_scatter_hist(ostream& s, bool verbose) {

     if (verbose) {
       s << "put_history of scatter" << flush << endl;
     }

     //     int i;
     for_all_scatter_hist(scatter_hist, this->get_first(), hi) {
       // s << *hi << " \t";
       //       for (int i = 0; i < Starlab::max(0, 7 - strlen(dummy_string)); i++)
       //	 cerr << " ";

       s << *hi << endl;

       //       for(i=0; i<N_RHO_ZONE_MAX; i++) 
       //	 s << hi->get_nhits(i);
       //       for(i=0; i<2-floor(log10(hi->get_nhits(i)+0.1)); i++)
       //	 s << " ";
       //       cerr << endl;
     }
       
}
scatter_hist* scatter_hist::contains(scatter_hist *hi) {

  for_all_scatter_hist(scatter_hist, get_first(), ha) {
    if(hi->identical_to(ha)) { 
      return ha;
    }
  }

  return NULL;
}

void scatter_hist::add_scatter_hist(scatter_exp he, int zone) {

  scatter_hist *hi = new scatter_hist(he);
  if (!hi) {
    cerr << "hi==NULL in add_scatter_hist(istream& s)" << endl;
    exit(-1);
  }

    if(hi->get_initial_form() == NULL || hi->get_final_form() == NULL)
      return;

    scatter_hist *identical_scenario = this->contains(hi);
    
    if (identical_scenario == NULL) {
      int n_scenario = 0;
      for_all_scatter_hist(scatter_hist, get_first(), ho) {
	n_scenario = Starlab::max(n_scenario, ho->get_id_scenario());
      }

      hi->set_id_scenario(n_scenario+1);
      hi->set_n_found(1);
      for(int i=0; i< N_RHO_ZONE_MAX; i++) 
	hi->n_hits[i] = 0;
      hi->n_hits[zone] = 1;

      scatter_hist *ha = get_last();
      ha->set_future(hi);
      hi->set_past(ha);
    }
    else {
      identical_scenario->inc_n_found(); 
      identical_scenario->n_hits[zone]++;
    }

}

void scatter_hist::add_scatter_hist(istream& s) {

    scatter_hist *hi = new scatter_hist(cin);
    if (!hi) {
      cerr << "hi==NULL in add_scatter_hist(istream& s)" << endl;
      exit(-1);
    }

    scatter_hist *identical_scenario = this->contains(hi);
    
    if (identical_scenario == NULL) {
      int n_scenario = 0;
      for_all_scatter_hist(scatter_hist, get_first(), ho) {
	n_scenario = Starlab::max(n_scenario, ho->get_id_scenario());
      }
      hi->set_id_scenario(n_scenario+1);
      hi->set_n_found(1);
      hi->n_hits[n_zone] = 1;

      scatter_hist *ha = get_last();
      ha->set_future(hi);
      hi->set_past(ha);
    }
    else {
      //      delete hi;
      identical_scenario->inc_n_found(); 
      identical_scenario->n_hits[n_zone]++;
    }
}

void scatter_hist::put_state(ostream & s) {
 
  s << "put_state" << endl;
  s << id_scenario << " " << initial_form << " " << final_form << endl;
}

void put_state(scatter_hist * hi, ostream & s) {

  for_all_scatter_hist(scatter_hist, hi, ha) {
    ha->put_state(s);
  }

  s << endl;
}

#ifdef TOOLBOX

void main(int argc, char **argv) {

    check_help();
    extern char *poptarg;
    int c;
    char* param_string = "c:";

    real cpu_time_check;
    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'c': cpu_time_check = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
	}            

    scatter_hist *hi = new scatter_hist(cin);
    if (!hi) {
      cerr << "hi==NULL in main" << endl;
      exit(-1);
    }

    do {

      hi->add_scatter_hist(cin);
    }
    while (!cin.eof());

    hi->put_scatter_hist(cerr);
}


#endif
