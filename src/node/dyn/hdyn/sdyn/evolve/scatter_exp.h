
//// scatter_hist: keeps track of scatter experiments
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Aug 2000   Simon Portegies Zwart   spz@mit.edu
//                                                  mit
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "stdinc.h"
#include "sdyn.h"
#include "string.h"
#ifdef HAS_MPI
#include  <mpi++.h>
#else
#include "localmpi++.h"
#endif
//
// scatter_exp.h
//

#ifndef    _SCATTER_EXP
#  define  _SCATTER_EXP

#define N_RHO_ZONE_MAX 20                 // Number of scatter annuli
#define N_STEP_BIN 10                    // Number of time-step bins
#define N_OSC_BIN 15                     // Number of oscillation bins

#if 0
enum scatter_discriptor {preservation=0, 
			 exchange, ionization, collision, 
			 single, binary, triple,
                         unknown, error, stopped,
			 number_of_scatter_discriptors};
#endif

enum scatter_discriptor {preservation=0, 
			 exchange, exchange_ionization,
			 single_ionization, multiple_ionization, 
			                    total_ionization,
			 collision_binary, collision_ionization,
			 stable_triple, stable_higher_order, 
			 two_body_collision, multiple_body_collision,
			 unknown, error, stopped,
			 not_specified,
			 number_of_scatter_discriptors};

#if 0
struct scatter_experiment_summary {

  char final_form[255];
  real energy_error;
  scatter_discriptor sd;
  int stop;
  int final_bound;
  int n_changes;
};
#endif

class scatter_exp {
    protected:

    char initial_form[255];
    char final_form[255];

    real time;
    real min_min_ssd;              // minimum distance squared
    real min_min_time;             // time at mimimun distance


    // Add information about minimal distance to next stellar surface
    // and information about which star is near to what other.
    // That information should be initialized in scatter.C 
    
    real *min_nn_dr2;     // minimum-ever nearest neighbor distance squared
    int *min_nn_label;    // identity of nearest-ever neighbor
    int *min_nn_of_label; // we should also keep track of the one who'se nn this is..

    real energy_error;
    real angular_momentum_error;
    real sigma;                   // breakdown of cross-section by type
    real sigma_err_sq;            // breakdown of squared errors by type

    scatter_discriptor sd;        // Is this considered an MPI::INT ?

    int n_star;
    int n_steps;
    int form_changes;

    int n_found;
    int  id_scenario;

    int n_zone;
    int n_hits[N_RHO_ZONE_MAX];
    int n_hit[N_RHO_ZONE_MAX];

    int  step_counter[N_STEP_BIN];
    int  osc_counter[N_OSC_BIN];

    int  n_initial;               //  5
    int  n_final;                 //  5
    int  n_coll;

    int final_bound;
    int resonance;
    int stop;

#if 0
    bool final_bound;
    bool resonance;
    bool stop;
#endif

    public:
       scatter_exp() {
	 initialize_to_zero();
       }
       scatter_exp(const scatter_exp &exp) {

	 time = exp.time;
	 min_min_ssd = exp.min_min_ssd;
	 min_min_time = exp.min_min_time;
	 final_bound = exp.final_bound;
	 n_coll = exp.n_coll;
	 sd = exp.sd;
	 n_star = exp.n_star;
	 n_zone = exp.n_zone;
	 energy_error = exp.energy_error;
	 angular_momentum_error = exp.angular_momentum_error;
	 form_changes = exp.form_changes;
	 resonance = exp.resonance;
	 stop = exp.stop;
	 int i;
	 for(i=0; i< N_OSC_BIN; i++) 
	   osc_counter[i] = exp.osc_counter[i];
	 for(i=0; i< N_STEP_BIN; i++) 
	   step_counter[i] = exp.step_counter[i];

	 for(i=0; i< N_RHO_ZONE_MAX; i++) {
	   n_hits[i] = exp.n_hits[i];
	   n_hit[i] = exp.n_hit[i];
	 }

	 //	 initial_form = new char[255];
	 strcpy(&initial_form[0], &exp.initial_form[0]);
	 //	 final_form = new char[255];
	 strcpy(&final_form[0], &exp.final_form[0]);
	 n_initial = exp.n_initial;
	 n_final = exp.n_final;
	 id_scenario = exp.id_scenario;
	 n_found = exp.n_found;
	 sigma_err_sq = exp.sigma_err_sq;
	 sigma = exp.sigma;

	 // copy nn's
	 // Note that n_star can be smaller than n_initial and larger than n_final
	 min_nn_dr2 = new real[n_initial];
	 min_nn_label = new int[n_initial];
	 min_nn_of_label = new int[n_initial];
	 for(i=0; i<n_initial; i++) {
	   min_nn_dr2[i] = exp.min_nn_dr2[i];
	   min_nn_label[i] = exp.min_nn_label[i];
	   min_nn_of_label[i] = exp.min_nn_of_label[i];
	 }
       }

       scatter_exp(istream &s) {
	 initialize_to_zero();

	 int id;

	 //	 char *iform = new char[255];
	 //	 char *fform = new char[255];
	 char iform[255];
	 char fform[255];
	 s >> id >> iform >> fform; 
	 //    cerr << iform << " " << fform << endl;
	 set_id_scenario(id);
	 set_final_form(fform);
	 set_initial_form(iform);
       }

       ~scatter_exp() {
	 //	 if(initial_form) delete []initial_form;
	 //	 if(final_form) delete []final_form;
	 if(min_nn_dr2)
	   delete []min_nn_dr2;
	 if(min_nn_label)
	   delete []min_nn_label;
	 if(min_nn_of_label)
	   delete []min_nn_of_label;
       }

  void initialize_to_zero();

  void set_min_min_ssd(real m) {min_min_ssd = m;}
  real get_min_min_ssd() {return min_min_ssd;}
  real get_min_min_time() {return min_min_time;}

  real *get_min_nn_dr2() {return min_nn_dr2;}
  int *get_min_nn_label() {return min_nn_label;}
  int *get_min_nn_of_label() {return min_nn_of_label;}

  void set_nstar(int n) {n_star = n;}
  int get_nstar() {return n_star;}
  int get_n_found() {return n_found;}
  void set_n_found(int nf) {n_found = nf;}
  void inc_n_found() {n_found++;}
  int get_id_scenario() {return id_scenario;}
  void set_id_scenario(int id) { id_scenario = id;}

  void  set_n_coll(int n) {n_coll = n;}
  int get_n_coll() {return n_coll;}

  void  set_final_bound(bool fb) {final_bound = fb;}
  bool get_final_bound() {return final_bound;}

  char *get_initial_form() {return &initial_form[0];}
  void set_initial_form(char *nf) {strcpy(&initial_form[0], nf);}
  char *get_final_form() {return final_form;}
  void set_final_form(char *nf) {strcpy(&final_form[0], nf);}
  int  get_n_initial() {return n_initial;}
  void set_n_initial(int n) {n_initial = n;}
  int  get_n_final() {return n_final;}
  void set_n_final(int n) {n_final = n;}
  void set_nzone(int z) {n_zone = z;}
  int  get_nzone() {return n_zone;}

  real get_time() {return time;}
  real *get_time_ptr() {return &time;}
  void set_energy_error(real e) {energy_error = e;}
  real get_energy_error() {return energy_error;}
  void set_angular_momentum_error(real l) {angular_momentum_error = l;}
  real get_angular_momentum_error() {return angular_momentum_error;}

  void set_n_steps(int n) {n_steps = n;}
  int  get_n_steps() {return n_steps;}
  int get_form_changes() {return form_changes;}
  void inc_form_changes() {form_changes++;}
  //  void set_scatter_discriptor(int s, scatter_discriptor h,
  //                                     scatter_discriptor v) {sd[h][v] = s;}
  void set_scatter_discriptor(scatter_discriptor s) {sd = s;}
  scatter_discriptor get_scatter_discriptor()  {return sd;}
  scatter_discriptor *get_scatter_discriptor_ptr()  {return &sd;}
  bool identical_to(scatter_exp *ha);
  void init_scatter_exp(sdyn* b);
  void final_scatter_exp(sdyn* b);
  void set_resonance(bool r) {resonance = r;}
  bool get_resonance() {return resonance;}

  void set_stop(bool s) {stop = s;}
  bool get_stop() {return stop;}

  real get_specific_sigma_err_sq(scatter_discriptor d);
  real get_specific_sigma(scatter_discriptor d);
  int get_specific_counts(scatter_discriptor d);
  int get_specific_counts(scatter_discriptor d, int i);


  void inc_sigma(real s) {sigma += s;}
  void inc_sigma_err_sq(real e2) {sigma_err_sq += e2;}

  void inc_step_counter(int index) {step_counter[index]++;}
  int get_step_counter(int index) {return step_counter[index];}
  void inc_osc_counter(int index) {osc_counter[index]++;}
  int get_osc_counter(int index) {return osc_counter[index];}

  void  set_nhits(int i, int n) {n_hits[i] = n;}
  int  get_nhits(int i) {return n_hits[i];}
  void inc_n_hits(int i) {if(1) n_hits[i]++;}
  void inc_n_hit(int i) {if(1) n_hit[i]++;}

  scatter_discriptor classify_scatter();
  //  scatter_discriptor classify_scatter(char*);
  int count_character_in_string(char string[], char search);
  int count_character_in_string(char string[], char search[], int n_char);
  int count_multiplicity();
 
  bool exch_or_exchion(char string[], char search1, char search2); 
  bool check_for_exchange_ionization();
 

  friend ostream& operator<<(ostream& s, scatter_exp&);
  bool operator == (scatter_exp& ha) const;
  bool operator != (scatter_exp& ha) const;

  MPI_Datatype initialize_data_structures_MPI();
};

scatter_exp* read_scatter_exp(istream& s);

#endif
