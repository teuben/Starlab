
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  hdyn.h: derived class for nbody systems, when using a Hermite integrator
//.............................................................................
//    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
//    version 2:
//.............................................................................
//     This file includes:
//		1) definition of class hdyn
//		2) declaration of member functions
//		3) declaration of related functions
//.............................................................................

#ifndef  STARLAB_HDYN_H
#  define  STARLAB_HDYN_H

#define MAX_PERTURBERS				256
#define ALLOW_LOW_LEVEL_PERTURBERS		true
#define RESOLVE_UNPERTURBED_PERTURBERS		false

// For full_dump binary ouptut:

#define SMALL_TDYN_PERT_SQ 1.e-4

//#include "sint.h"

#include  "_dyn_.h"

class kira_counters;	// (defined later)
class kira_options;
class kira_diag;

// class hdyn;
// typedef int (*acc_function_ptr)(hdyn **, int, xreal, bool);

//-----------------------------------------------------------------------------
//  hdyn  --  a derived class of dynamical particles, with enough information
//             to integrate the equations of motions with a 4th-order Hermite
//             polynomial integrator.
//-----------------------------------------------------------------------------

class  hdyn : public _dyn_ {
    protected:

        // ---------------------- Global variables! ----------------------
	// --------------- (initialized in util/hdyn_io.C) ---------------

	// Counters and diagnostics:

        static kira_counters	*kc;
	static kira_options	*options;
	static kira_diag	*diag;

	// Binary evolution:

        static bool use_dstar;  // activate binary evolution if true

	// Stellar encounters and mergers:

        static real stellar_encounter_criterion_sq;
        static real stellar_merger_criterion_sq;
        static real stellar_capture_criterion_sq;

	// Perturbed top-level binaries:

	static hdyn **perturbed_list;	// list of binaries
	static int n_perturbed;		// number of binaries on the list

	// Run-time integration parameters:

	static real eta;	// time step parameter
	static real eps;	// softening length
	static real eps2;	// softening length squared

        static real d_min_fac;	// scale term governing tree adjustment
        static real d_min_sq;	// length scale governing tree adjustment
	static real lag_factor;	// squared hysteresis factor
        static real mbar;	// mass scale

        static real gamma2;	// squared threshhold for unperturbed motoin
        static real gamma23;	// gamma^{-2/3}

	static real initial_step_limit;	// limit on first time step
	static real step_limit;		// limit on all time steps
	static real unpert_step_limit;	// limit on unperturbed time steps
					// to allow binary evolution to be
					// regularly updated

	// Escaper removal:

	static real scaled_stripping_radius; // stripping radius for unit mass
  
	// Slow binary motion:

	static int  max_slow_factor;
	static real max_slow_perturbation;
	static real max_slow_perturbation_sq;

	// Hardware/software configuration:

	static unsigned int config;	// 0: normal, 1: GRAPE-4, 2: GRAPE-6

//	static acc_function_ptr kira_calculate_top_level_acc_and_jerk;
//					// function for force computation

	// ---------------------------------------------------------------

	// Run-time integration flag.  Used in correct_acc_and_jerk to
	// determine whether or not a node is on the current integration
	// list and hence should be considered for correction.  Adopted
	// in order to avoid checking get_next_time(), which can be slow,
	// especially if time is an xreal.

	bool on_integration_list;

        // Variables for unperturbed kepler motion:

        real  perturbation_squared;  // relative perturbation squared
        bool  fully_unperturbed;     // true if orbit is fully unperturbed
        real  unperturbed_timestep;  // timestep for the unperturbed motion

        // Perturber information:

        int  n_perturbers;	// number of perturbers of this node
        hdyn** perturber_list;  // pointer to perturber array
        bool valid_perturbers;  // true if any particle is within the
			        // perturbation radius and the perturber
				// has not overflowed
        real  perturbation_radius_factor;

	int n_perturbers_low;	// number of perturbers of lower-level nodes
				// -- for use in multiples when a top-level
				// binary has too many perturbers, but still
				// need a list of perturbers for daughter nodes
				//				(Steve, 3/03)
        bool valid_perturbers_low;  // true low-level list is usable

	real posvel;		// convenient (?) to save pos*vel, at least
				// for low-level nodes
	real prev_posvel;

        // Other (colliding) neighbor information:

        hdyn* nn;               // pointer to the nearest neighbor
        real d_nn_sq;           // distance squared to the nearest neighbor
    
        hdyn* coll;		// pointer to the neighbor whose surface is
                                // closest to 'this' surface (i.e. as
                                // measured by  distance - r1 - r2 )
        real d_coll_sq;		// distance squared to this neighbor
  
        // Variables for GRAPE-4/6:

        int grape_index;	// address of this particle in GRAPE memory

        real grape_rnb_sq;      // the neighbor sphere radius of the particle
			        // used on the GRAPE
	int grape_nb_count;	// counter to limit computation of full
				// neighbor information (coll and perturbers)

	// Bookkeeping (keeping track of costs):

	real steps;
	real direct_force;
	real indirect_force;

        // Protected functions - criteria for a multiple system:

        bool is_weakly_perturbed(int& status);
        bool is_stable(int& status, bool top_level = true);

    public:

        hdyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase)
	    : _dyn_(the_hbpfp, the_sbpfp) {

	    on_integration_list = false;

	    perturbation_squared = -1;
            unperturbed_timestep  =  -VERY_LARGE_NUMBER;
            fully_unperturbed = false;

            n_perturbers = n_perturbers_low = 0;
	    perturber_list = NULL;
            valid_perturbers = valid_perturbers_low = false;
	    perturbation_radius_factor = -1;

	    posvel = VERY_LARGE_NUMBER;
	    prev_posvel = -VERY_LARGE_NUMBER;	// avoid spurious apocenter
						// in a hyperbolic encounter

            nn = NULL;
            d_nn_sq = 0;

            coll = NULL;
	    d_coll_sq = 0;

	    grape_index = 0;
	    grape_rnb_sq = 0;
	    grape_nb_count = 0;

	    steps = direct_force = indirect_force = 0;

        }

	virtual ~hdyn() {
	    if (perturber_list)
		delete [] perturber_list;	// added 7/29/98 (SLWM & JM)
	}

        // ------------------ Global variable accessors ------------------

	// Options, counters, and diagnostics:

        inline kira_counters *get_kira_counters()  {return kc;}
        void set_kira_counters(kira_counters *k){kc = k;}

        inline kira_options *get_kira_options()	{return options;}
        void set_kira_options(kira_options *o)	{options = o;}

        inline kira_diag *get_kira_diag()	{return diag;}
        void set_kira_diag(kira_diag *d)	{diag = d;}

	// Binary evolution:

	inline bool get_use_dstar()		{return use_dstar;}
	void set_use_dstar(bool u)		{use_dstar = u;}

	// Stellar encounters and mergers:

        real get_stellar_encounter_criterion_sq()
		{return stellar_encounter_criterion_sq;} 
        void set_stellar_encounter_criterion_sq(real d_sq)
		{stellar_encounter_criterion_sq = d_sq;}
        inline real get_stellar_capture_criterion_sq()
		{return stellar_capture_criterion_sq;} 
        void set_stellar_capture_criterion_sq(real d_sq)
		{stellar_capture_criterion_sq = d_sq;}

        inline real get_stellar_merger_criterion_sq()
		{return stellar_merger_criterion_sq;} 
        void set_stellar_merger_criterion_sq(real d_sq)
		{stellar_merger_criterion_sq = d_sq;}

	// Perturbed top-level binaries:

	inline int get_n_perturbed() {return n_perturbed;};
	inline hdyn** get_perturbed_list() {return perturbed_list;};

	// Run-time integration parameters:

	void set_eta(real e)			{eta = e;}
	inline real get_eta()			{return eta;}

	void set_eps(real e)			{eps = e;
					         eps2 = e*e;}
	real get_eps()				{return eps;}
	inline real get_eps2()			{return eps2;}

        void set_d_min_fac(real d)		{d_min_fac = d;}
        inline real get_d_min_fac()		{return d_min_fac;}

        void set_d_min_sq(real d)		{d_min_sq = d;}
        inline real get_d_min_sq()		{return d_min_sq;}

	void set_lag_factor(real f)		{lag_factor = f;}
	inline real get_lag_factor()		{return lag_factor;}

        void set_mbar(real m)			{mbar = m;}
        inline real get_mbar()			{return mbar;}

        void set_gamma2(real g2)		{gamma2 = g2;
					         gamma23 = pow(g2, -1.0/3.0);}
        inline real get_gamma2()		{return gamma2;}
        inline real get_gamma23()		{return gamma23;}

        void set_initial_step_limit(real s)	{initial_step_limit = s;}
        inline real get_initial_step_limit()	{return initial_step_limit;}

        void set_step_limit(real s)		{step_limit = s;}
        inline real get_step_limit()		{return step_limit;}

        void set_unpert_step_limit(real s)	{unpert_step_limit = s;}
        inline real get_unpert_step_limit()	{return unpert_step_limit;}

	// Escaper removal:

	void set_scaled_stripping_radius(real r)
	    {scaled_stripping_radius = r;}
	inline real get_scaled_stripping_radius()
	    {return scaled_stripping_radius;}

	// Slow binary motion:

	void set_max_slow_factor(int f = 1)	{max_slow_factor = f;}
	inline int get_max_slow_factor()	{return max_slow_factor;}

	void set_max_slow_perturbation(real p) {
	    max_slow_perturbation = p;
	    max_slow_perturbation_sq = p*p;
	}
	void set_max_slow_perturbation_sq(real p) {
	    max_slow_perturbation_sq = p;
	    max_slow_perturbation = sqrt(p);
	}
	inline real get_max_slow_perturbation() {
	    return max_slow_perturbation;
	}
	inline real get_max_slow_perturbation_sq() {
	    return max_slow_perturbation_sq;
	}

	// Hardware/software configuration:

	inline unsigned int get_config()	{return config;}
	inline void set_config(unsigned int c)	{config = c;}

// 	inline acc_function_ptr get_kira_calculate_top_level_acc_and_jerk()
// 	    {return kira_calculate_top_level_acc_and_jerk;}
//
// 	inline void
// 	       set_kira_calculate_top_level_acc_and_jerk(acc_function_ptr f)
// 		    {kira_calculate_top_level_acc_and_jerk = f;}

        // ----------------- Other member data accessors ------------------

	// Integration flag:

	inline bool is_on_integration_list()	{return on_integration_list;}
	void set_on_integration_list(bool on = true) {on_integration_list = on;}
	void clear_on_integration_list()	{on_integration_list = false;}

        // Unperturbed kepler motion:

        inline real get_perturbation_squared()	{return perturbation_squared;}
        void set_perturbation_squared(real p)	{perturbation_squared = p;}

        inline bool get_fully_unperturbed()	{return fully_unperturbed;}
        void set_fully_unperturbed(bool f)	{fully_unperturbed = f;}
	
        inline real get_unperturbed_timestep()   {return unperturbed_timestep;}

	void set_unperturbed_timestep(real step) {unperturbed_timestep
							     = step;}

        // Perturber information:

        void zero_perturber_list()		{n_perturbers
						     = n_perturbers_low = 0;}
        void set_n_perturbers(int n) {
	    if (is_leaf())
		cerr << "warning: setting n_perturbers for leaf "
		     << format_label() << endl;
	    n_perturbers = n;
	}
        void set_n_perturbers_low(int n) {
	    if (is_leaf())
		cerr << "warning: setting n_perturbers_low for leaf "
		     << format_label() << endl;
	    n_perturbers_low = n;
	}
        inline int get_n_perturbers()		{return n_perturbers;};
        inline int get_n_perturbers_low()	{return n_perturbers_low;};
        inline hdyn ** get_perturber_list()	{return perturber_list;}

        void remove_perturber_list()		{n_perturbers
						     = n_perturbers_low = 0;
					         valid_perturbers
						     = valid_perturbers_low
						     = false;
					         if (perturber_list)
						     delete [] perturber_list;
						 perturber_list = NULL;}

        void print_perturber_list(ostream & s = cerr, char* pre = "");
        void find_print_perturber_list(ostream & s = cerr, char* pre = "");

        void new_perturber_list() {
	    if (perturber_list == NULL) {
		perturber_list = new hdyn *[MAX_PERTURBERS];
	    }
	    n_perturbers = n_perturbers_low = 0;
	    valid_perturbers = valid_perturbers_low = false;
	}

        inline bool  get_valid_perturbers()       {return valid_perturbers;}
        inline bool  get_valid_perturbers_low()   {return valid_perturbers_low;}
        void  set_valid_perturbers(const bool vp) {valid_perturbers = vp;}
        void  set_valid_perturbers_low(const bool vp) {valid_perturbers_low
							   = vp;}

        inline real get_perturbation_radius_factor()
		{return perturbation_radius_factor;}
        void set_perturbation_radius_factor(real f) 
		{perturbation_radius_factor = f;}

	inline bool passed_apo()		{return (posvel < 0
							 && prev_posvel >= 0);}

	inline real get_posvel()		{return posvel;}

	// NOTE: This may run into problems with nearly circular,
	//       lightly perturbed orbits...

        // Neighbor information:

	inline hdyn * get_nn()			{return nn;}
	void set_nn(hdyn * new_nn)		{nn = new_nn;}
	inline real get_d_nn_sq()		{return d_nn_sq;}
	void set_d_nn_sq(real d)		{d_nn_sq = d;}
  
        inline hdyn * get_coll()		{return coll;}
        void  set_coll(hdyn * new_coll)		{coll = new_coll;}
        inline real get_d_coll_sq()		{return d_coll_sq;}
        void  set_d_coll_sq(real d)		{d_coll_sq = d;}

        // Variables for GRAPE-4/6

        inline int get_grape_index()		{return grape_index;}
        void set_grape_index(int hindex)	{grape_index = hindex;}
        inline real get_grape_rnb_sq()		{return grape_rnb_sq;}
        void set_grape_rnb_sq(real rnb_sq)	{grape_rnb_sq = rnb_sq;}
        inline int get_grape_nb_count()		{return grape_nb_count;}
        void set_grape_nb_count(int n)		{grape_nb_count = n;}

	bool has_grape()			{return (config > 0);}
	bool has_grape4()			{return (config > 0
							 && config%2 != 0);}
	bool has_grape6()			{return (config > 0
							 && config%2 == 0);}

	// Experimental:

	inline real *get_pos_addr()		{return &pos[0];}
	inline real *get_vel_addr()		{return &vel[0];}
	inline real *get_acc_addr()		{return &acc[0];}
	inline real *get_jerk_addr()		{return &jerk[0];}
	inline real *get_old_acc_addr()		{return &old_acc[0];}
	inline real *get_old_jerk_addr()	{return &old_jerk[0];}
	inline real *get_k18_addr()		{return &k_over_18[0];}

	// Bookkeeping:

	inline void inc_steps(int i = 1)	{steps += i;}
	inline void inc_steps(real i)		{steps += i;}
	inline real get_steps()			{return steps;}
	inline void inc_direct_force(int i = 1)	{direct_force += i;}
	inline void inc_direct_force(real i)	{direct_force += i;}
	inline real get_direct_force()		{return direct_force;}
	inline void inc_indirect_force(int i = 1) {indirect_force += i;}
	inline void inc_indirect_force(real i)	{indirect_force += i;}
	inline real get_indirect_force()	{return indirect_force;}

	// Convenient restatements of inherited functions:

        inline hdyn * get_parent()
		{return (hdyn*) node::get_parent();}
        inline hdyn * get_oldest_daughter()
		{return (hdyn*)node::get_oldest_daughter();}
        inline hdyn * get_younger_sister()
		{return (hdyn*) node::get_younger_sister();}
        inline hdyn * get_elder_sister()
		{return (hdyn*) node::get_elder_sister();}

        inline hdyn * get_root()
		{return (hdyn*) node::get_root();}
        inline hdyn * get_top_level_node()
		{return (hdyn*) node::get_top_level_node();}
        inline hdyn * get_binary_sister()
		{return (hdyn*) node::get_binary_sister();}

	virtual void null_pointers();
	virtual void print_static(ostream &s = cerr);

        virtual istream& scan_dyn_story(istream&);

	virtual bool check_and_correct_node(bool verbose = true);

        virtual ostream& print_dyn_story(ostream &s,
					 bool print_xreal = true,
					 int short_output = 0);
	
        inline real get_true_timestep() {
            if (kep)
                return unperturbed_timestep;
            else
                return timestep;
        }

	// Potentially slow function -- avoid where possible!

        inline xreal get_next_time() {
            if (kep)
                return get_time() + (xreal)unperturbed_timestep;
            else
                return get_time() + (xreal)timestep;
        }

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// Member functions (declarations as expressed in named files):

	// ----- In util/hdyn_init.C: -----

	void initialize_unperturbed();
	void initialize_slow();

	// ----- In util/hdyn_kepler.C: -----

        void update_kepler_from_hdyn();
        void reinitialize_kepler_from_hdyn();

	// ----- In util/hdyn_tt: -----

	hdyn* next_node(hdyn*, bool resolve_unperturbed = true);
	hdyn* find_perturber_node();
        bool nn_stats(real energy_cutoff, real kT,
		      vec center, bool verbose,
		      bool long_binary_output = true,
		      int which = 0);
        real print_pert(bool long_binary_output = true,
			int indent = BIN_INDENT);

        void setup_binary_node();

	// ----- In hdyn_ev.C: -----

        void  synchronize_node(bool exact = true);

        void set_first_timestep(real additional_step_limit = 0);

        // bool correct();		// changed to bool by Steve, 9/15/98

        void update(vec& bt2,	// called only from correct_and_update;
		   vec& at3);	// member function for convenience

	bool correct_and_update();	// combined 8/99 by Steve

        int flat_calculate_acc_and_jerk(hdyn* root_node,
					bool make_perturber_list);

        void perturber_acc_and_jerk_on_leaf(vec & a,
					    vec & j,
					    real & p,
					    real & distance_squared,
					    hdyn*& nn,
					    hdyn* pnode,
					    hdyn* step_node);

	void tree_walk_for_partial_acc_and_jerk_on_leaf(hdyn * b,
							hdyn * mask,
							vec& offset_pos,
							vec& offset_vel,
							vec& a,
							vec& j,
							real & p,
							real & distance_squared,
							hdyn * &p_nn,
							bool point_mass_flag,
							hdyn * step_node);

	void calculate_partial_acc_and_jerk_on_leaf(hdyn * top,
						    hdyn * common,
						    hdyn * mask,
						    vec & a,
						    vec & j,
						    real & p,
						    real & distance_squared,
						    hdyn * &nn,
						    bool point_mass_flag,
						    hdyn * pnode,
						    hdyn * step_node);

	void calculate_partial_acc_and_jerk(hdyn * top,
					    hdyn * common,
					    hdyn * mask,
					    vec & a,
					    vec & j,
					    real & p,
					    real & distance_squared,
					    hdyn * &p_nn,
					    bool point_mass_flag,
					    hdyn* pnode,
					    hdyn* step_node);

	void check_add_perturber(hdyn* p, vec& this_pos);
	void create_low_level_perturber_list(hdyn* pnode);
	void create_low_level_perturber_lists(bool only_if_null = true);

	void calculate_acc_and_jerk_on_low_level_node();
        void calculate_acc_and_jerk_on_top_level_node(bool exact);

        void top_level_node_prologue_for_force_calculation(bool exact);
        int top_level_node_real_force_calculation();
        void top_level_node_epilogue_force_calculation();

        void calculate_acc_and_jerk(bool exact);

        void integrate_node(bool exact = true,
			    bool integrate_unperturbed = true,
			    bool force_unperturbed_time = false);

	// ----- In hdyn_tree.C: -----

        real distance_to_sister_squared();
        hdyn* new_sister_node(bool& top_level_combine);
        int adjust_tree_structure(int full_dump = 0);

	hdyn* check_periapo_node();

	hdyn* check_merge_node();
        hdyn* merge_nodes(hdyn * bcoll, int full_dump = 0);
	void  merge_logs_after_collision(hdyn *bi, hdyn* bj);


	// ----- In hdyn_unpert.C: -----

        void print_unperturbed_binary_params();

        void update_dyn_from_kepler(bool need_acc_and_jerk = true);
        bool is_close_pair();

        bool is_unperturbed_and_approaching();
        void startup_unperturbed_motion();

	// Accidentally overloaded (see real argument version above...)!

        real set_unperturbed_timestep(bool check_phase);
        // int set_unperturbed_timestep(bool check_phase);
        // sint set_unperturbed_timestep(bool check_phase);
        // unsigned long set_unperturbed_timestep(bool check_phase);

        real get_unperturbed_steps(bool to_apo = true,
				   bool predict = false);
        // int get_unperturbed_steps(bool to_apo = true,
	//			     bool predict = false);
        // unsigned long get_unperturbed_steps(bool to_apo = true,
	//				       bool predict = false);

	void recompute_unperturbed_step();
	void recompute_unperturbed_steps();

        int integrate_unperturbed_motion(bool& reinitialize,
					 bool force_time = false);

	// ----- In hdyn_slow.C: -----

        void startup_slow_motion();
        void extend_or_end_slow_motion(real P = 0);

	// ----- In perturbed_list.C: -----

	bool is_perturbed_cpt();
	bool is_perturbed_cm();
	void reconstruct_perturbed_list();
	void dump_perturbed_list();
	int  get_index_on_perturbed_list(bool debug = false);
	bool on_perturbed_list();
	void check_perturbed_list();
	void add_to_perturbed_list(int id = 0);
	void remove_from_perturbed_list(int id = 0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
};

//--------------------------------------------------------------------------

// The following files need hdyn to have been defined:

#include  "hdyn_kepler.h"	// not clear that it is really necessary
				// to separate out these declarations...
#include  "kira.h"

//--------------------------------------------------------------------------

#define  for_all_leaves_or_unperturbed(dyntype, base, node_name)              \
         for (dyntype* node_name = base;                                      \
	      node_name != NULL;                                              \
	      node_name = (dyntype*) node_name->next_node(base, false))       \
	      if (node_name->get_oldest_daughter() == NULL                    \
		  || (node_name->get_oldest_daughter()->get_kepler()	      \
		      && node_name->get_oldest_daughter()		      \
		      		  ->get_fully_unperturbed()))


// Non-member functions:

typedef hdyn * hdynptr; // to enable dynamic array declarations such as
			//    hdyn** hdyn_ptr_array = new hdynptr[n];
			// (note that the following expression is illegal:
			//    hdyn** hdyn_ptr_array = new (hdyn *)[n];)

// Third argument below is present because get_node expects it...

inline node * new_hdyn(hbpfp the_hbpfp,
		       sbpfp the_sbpfp,
		       bool use_stories)
    {return (node *) new hdyn(the_hbpfp, the_sbpfp);}	 // ignore 3rd argument

// See kira_init.C...

// get_hdyn() modified to follow Ernie's new get_dyn() (Steve, 5/03).
//
//  inline hdyn * get_hdyn(istream & s = cin,
//  		       hbpfp the_hbpfp = new_hydrobase,
//  		       sbpfp the_sbpfp = new_starbase)
//  {
//      hdyn *root = (hdyn *) get_node(s, new_hdyn, the_hbpfp, the_sbpfp);
//
//      // Complete the initialization of special configurations:
//
//      root->initialize_unperturbed();
//      root->initialize_slow();
//
//      return root;
//  }

inline  hdyn * get_hdyn(istream & s = cin,
			hbpfp the_hbpfp = new_hydrobase,
			sbpfp the_sbpfp = new_starbase)
{
    // Nasty kludge because of the C/C++ I/O problem...

  if (dyn::get_ifp() ?
      ungetc(getc(dyn::get_ifp()), dyn::get_ifp()) : s.peek() == '(') {
    // Old code:
    hdyn *root = (hdyn*)get_node(s, new_hdyn, the_hbpfp, the_sbpfp);
    // Complete the initialization of special configurations:
    root->initialize_unperturbed();
    root->initialize_slow();
    return root;
  }
  else return (hdyn*)get_col(s, new_hdyn, the_hbpfp, the_sbpfp, true);
}

#define  put_hdyn  put_node

inline hdyn * common_ancestor(hdyn * bi, hdyn * bj)
		{return (hdyn*) common_ancestor((node*)bi, (node*)bj);}

// ----- In util/hdyn_pp3.C: -----

void pp3_maximal(hdyn *, ostream & s = cerr, int i = 0);
void pp3(hdyn *, ostream & s = cerr, int i = 0);
void pp3_minimal(hdyn *, ostream & s = cerr, int i = 0);
void pp3_tree(hdyn *, ostream & s = cerr, int i = 0);

void pp3_maximal(char *, ostream & s = cerr, int i = 0);
void pp3(char *, ostream & s = cerr, int i = 0);
void pp3_minimal(char *, ostream & s = cerr, int i = 0);
void pp3_tree(char *, ostream & s = cerr, int i = 0);

// ----- In util/Qt_pp3.C: -----

class QApplication;
int Qt_pp3(hdyn *b,
	   QApplication *app);
int Qt_pp3(hdyn *b,
	   const char * disp = NULL);

// ----- In util/hdyn_tt.C: -----

void print_nn(hdyn* b, int level = 0, ostream& s = cerr);
void print_coll(hdyn* b, int level = 0, ostream& s = cerr);

// The second definition below, of hdyn_MF_ptr, should probably be local to
// the file hdyn_tt.C  -- we should look at this later.

typedef vec (hdyn::*hdyn_VMF_ptr)(void);   // vec member function pointer
typedef void (hdyn::*hdyn_MF_ptr)(const vec &);   // member function pointer

//  NOTE:  all the "hdyn" declarations below should be changed to "node"
//	   here (just as was done with get_node()).

void check_consistency_of_node(hdyn * node,
			       hdyn_VMF_ptr get_something,
			       char *id);
void check_consistency_of_nodes(hdyn * node);

void create_binary_from_toplevel_nodes(hdyn * bi, hdyn * bj);

void remove_node_and_correct_upto_ancestor(hdyn * ancestor, hdyn * node);

vec something_relative_to_ancestor(hdyn * bj,
				      hdyn * bi,
				      hdyn_VMF_ptr get_something);

// This is the same function as something_relative_to_root, but
// g++ versions after 2.7.0 complain about the ambiguity for pos,
// vel, and acc if we name it something_relative_to_root!
// This version is actually needed (and used) only for get_jerk.

// DEC C++ also complains about all uses of something_relative_to_root
// when hdyn member functions are involved...

vec hdyn_something_relative_to_root(hdyn * bi,
				       hdyn_VMF_ptr get_something);

void correct_leaf_for_change_of_mass(hdyn * node, real dm);

void correct_leaf_for_change_of_vector(hdyn * node,
				       vec d_something,
				       hdyn_VMF_ptr get_something,
				       hdyn_MF_ptr inc_something);

void move_node(hdyn * node_to_move, hdyn * place_to_insert);
void move_node_by_index(int i, int j, hdyn * root);

void calculate_energies(hdyn * root, real eps2,
			real & epot, real & ekin, real & etot,
			bool cm = false);

// ---- In util/reset_counters.C: -----

void reset_counters(hdyn* bi);

// ----- In correct_perturbers.C: -----

void correct_perturber_lists(hdyn * b,
			     hdyn ** list,
			     int n,
			     hdyn * cm = NULL);

void expand_perturber_lists(hdyn * b,
			    hdyn * bb,
			    bool verbose = false);

// ----- In hdyn_ev.C: -----

real kepler_step(hdyn *b, real correction_factor = 1);
real timestep_correction_factor(hdyn *b);

void update_binary_sister(hdyn* bi);

hdyn* node_to_move(hdyn* b, real& tmin);
void get_nodes_to_move(hdyn * b,
		       hdyn * list[],
		       int &nlist,
		       real & tmin);

void initialize_system_phase1(hdyn* b, xreal t);

void correct_acc_and_jerk(hdyn * root,		// OLD!
			  bool& reset);

void correct_acc_and_jerk(hdyn** next_nodes,	// NEW
			  int n_next);

// ----- In hdyn_grape4/6.C: -----

int get_grape4_chip(hdyn *b);
void grape6_set_dma(bool jp_dma = true);

void check_release_grape(unsigned int config,
			 kira_options *ko, xreal time, bool verbose = true);
void check_release_grape4(kira_options *ko, xreal time, bool verbose = true);
void check_release_grape6(kira_options *ko, xreal time, bool verbose = true);

void grape4_calculate_energies(hdyn *b,
			       real &epot,
			       real &ekin,
			       real &etot,
			       bool cm = false);
void grape6_calculate_energies(hdyn *b,
			       real &epot,
			       real &ekin,
			       real &etot,
			       bool cm = false);

int grape4_calculate_acc_and_jerk(hdyn **next_nodes,
				  int n_next,
				  xreal time,
				  bool restart);
int grape6_calculate_acc_and_jerk(hdyn **next_nodes,
				  int n_next,
				  xreal time,
				  bool restart);

bool grape4_calculate_densities(hdyn *root,
				real h2_crit = 4);
bool grape6_calculate_densities(hdyn *root,
				real h2_crit = 4);

// ----- In util/hdyn_io.C: -----

void set_write_unformatted(bool u = true);
bool get_write_unformatted();
void set_complete_system_dump(bool d = true);
void set_hdyn_check_timestep(bool check = true);


// ----- In hdyn_schedule.C: -----

void print_sort_counters();

void fast_get_nodes_to_move(hdyn * b,
			    hdyn * list[],
			    int  & nlist,
			    xreal & tmin,
			    bool & reset);

void dump_node_list(int n = 1000000000);
void dump_node_list_for(char *s);

// ----- In hdyn_slow.C: -----

bool has_binary_perturbers(hdyn* b);
void clear_perturbers_slow_perturbed(hdyn *b);
void list_slow_perturbed(hdyn* b);
void check_slow_consistency(hdyn* b);

// ----- In hdyn_tree.C: -----

void synchronize_tree(hdyn * b);

void split_top_level_node(hdyn * bi, int full_dump = 0);
void combine_top_level_nodes(hdyn * bj, hdyn * bi, int full_dump = 0);

// ----- In hdyn_unpert.C: -----

void set_allow_unperturbed(hdyn *b, bool value = true);
void set_allow_multiples(hdyn *b, bool value = true);
void toggle_unperturbed(hdyn *b, int level);
void print_unperturbed_options(hdyn *b);
void enable_isolated_multiples(bool value = true);

// ----- In kira.C: -----

void set_n_top_level(hdyn *b);
int  get_n_top_level();

// ----- In kira_approx.C: -----

real pow_approx(real x);

// ----- In kira_config.C: -----

unsigned int kira_config(hdyn *b,
			 int force_config = -1);

void kira_print_config(unsigned int config);

// ----- In kira_counters.C: -----

void initialize_counters_from_log(hdyn* b);
void write_counters_to_log(hdyn* b);
void print_counters(kira_counters* kc,
		    bool print_all = true,
		    kira_counters* kc_prev = NULL);

// ----- In kira_density.C: -----

bool kira_calculate_densities(hdyn* b, vec& cod_pos, vec& cod_vel);

// ----- In kira_energy.C: -----

void kira_calculate_internal_energies(hdyn* b,
				      real& epot, real& ekin, real& etot,
				      bool cm = false,
				      bool use_grape = true);

void kira_calculate_energies(dyn* b, real eps2, 
			     real &potential, real &kinetic, real &total,
			     bool cm);

void kira_top_level_energies(dyn *b, real eps2,
			     real& potential_energy,
			     real& kinetic_energy);

void print_recalculated_energies(hdyn* b,
				 bool print_dde = false,
				 bool save_story = false);

void calculate_energies_with_external(hdyn* b,
				      real& epot, real& ekin, real& etot,
				      bool cm = false,
				      bool use_grape = true);
// ----- In kira_escape.C: -----

void check_and_remove_escapers(hdyn* b,
			       xreal& ttmp, hdyn** next_nodes,
			       int& n_next, bool& tree_changed);

// ----- In kira_ev.C: -----

int kira_calculate_top_level_acc_and_jerk(hdyn **next_nodes,
				    int n_next,
				    xreal time,
				    bool restart_grape);

int top_level_acc_and_jerk_for_list(hdyn **next_nodes,
				    int n_next,
				    xreal time,
				    bool restart_grape);

int calculate_acc_and_jerk_for_list(hdyn ** next_nodes,
				    int n_next,
				    xreal time,
				    bool exact,
				    bool tree_changed,
				    bool & reset_force_correction,
				    bool & restart_grape);

void calculate_acc_and_jerk_on_all_top_level_nodes(hdyn * b);
void calculate_acc_and_jerk_on_top_level_binaries(hdyn * b);

void kira_synchronize_tree(hdyn *b, bool sync_low_level = false);

void initialize_system_phase2(hdyn * b,
			      int id = 0,
			      int set_dt = 1);

// ----- In kira_check.C: -----

bool check_kira_flag(hdyn* b, char* kira_flag);
bool check_allowed(bool allow_kira_override,
			 char * what_is_allowed,
			 bool verbose, bool& need_skip);

// ----- In kira_init.C: -----

bool kira_initialize(int argc, char** argv,
		     hdynptr& b,	// hdyn root node
		     real& delta_t,	// time span of the integration
		     real& dt_log,	// output interval--power of 2
		     int&  long_binary_out, // frequency of full binary info
		     real& dt_snap,	// snap output interval
		     real& dt_sstar,	// stellar evolution timestep
		     real& dt_esc,	// escaper removal
		     real& dt_reinit,	// reinitialization interval
		     real& dt_fulldump,	// full dump interval
		     bool& exact,	// no perturber list if true
		     real& cpu_time_limit,
		     bool& verbose,
		     bool& snap_init,
		     bool& save_last_snap,
		     char* snap_save_file,
		     int& n_stop,	// n to terminate simulation
		     bool& alt_flag);	// enable alternative output

// In kira.C:

void kira(hdyn * b,	       // hdyn array
	  real delta_t,	       // time span of the integration
	  real dt_log,	       // time step of the integration
	  int long_binary_out,
	  real dt_snap,	       // snapshot output interval
	  real dt_sstar,       // timestep for single stellar evolution
	  real dt_esc,	       // escaper removal
	  real dt_reinit,      // reinitialization interval
	  real dt_fulldump,    // full dump interval
	  bool exact,	       // exact force calculation
	  real cpu_time_limit,
	  bool verbose,
	  bool snap_init,
	  bool save_snap_at_log, // save snap at log output
	  char* snap_save_file,  // filename to save in
	  int n_stop,	         // when to stop
	  bool alt_flag);	 // enable alternative output

void kira_finalize(hdyn *b);

// In kira_int.C:

int integrate_list(hdyn * b,
		   hdyn ** next_nodes, int n_next,
		   bool exact, bool & tree_changed,
		   int& n_list_top_level,
		   int full_dump,
		   real r_reflect);

// ----- In kira_log.C: -----

int get_effective_block(real dt);
void print_dstar_params(dyn* b);
bool print_dstar_stats(dyn* b, bool mass_function,
		       vec center, bool verbose);
void get_energies_with_external(dyn* b, real eps2, 
				real &kinetic, real &potential, real &total,
				bool cm = false);
void print_statistics(hdyn* b, int long_binary_output = 2);
void set_block_check(bool b = true);

void update_cpu_counters(hdyn * b);
void log_output(hdyn * b,
		real count, real steps,
		real count_top_level, real steps_top_level,
		kira_counters* kc_prev,
		int long_binary_output = 2);
void force_nodensity();

void snap_output(hdyn * b, real steps, int& snaps,
		 bool reg_snap, bool last_snap,
		 char * snap_save_file,
		 xreal t, xreal ttmp, real t_end,
		 real& t_snap, real dt_snap, bool verbose);

// ----- In kira_runtime.C: -----

void clean_up_files();
bool check_file(char* name,
		bool del = true);

void check_kira_init(hdyn* b);

bool check_kira_runtime(hdyn* b,
			real& t_end, real& new_dt_log, real& new_dt_snap,
			int& long_binary_output, char* new_snap_save_file,
			bool& tree_changed);

// ----- In kira_stellar.C: -----

bool evolve_stars(hdyn* b, int full_dump = 0);

// ----- In kira_encounter.C: -----

real get_sum_of_radii(hdyn* bi, hdyn* bj, bool check_story = false);
real print_encounter_elements(hdyn* bi, hdyn* bj,
			      char* s = "Collision",
			      bool verbose = true);

void check_print_close_encounter(hdyn *bi);

void print_close_encounter(hdyn* bi,
			   hdyn* bj);

// ----- In kira_smallN.C -----

void set_smallN_eta(real a);
void set_smallN_gamma(real g);
void set_smallN_niter(int n);
void set_smallN_dtcrit(real dt);
void set_smallN_rcrit(real r);

real get_smallN_eta();
real get_smallN_gamma();
int  get_smallN_niter();
real get_smallN_dtcrit();
real get_smallN_rcrit();

real smallN_evolve(hdyn *b,
		   real t_end = VERY_LARGE_NUMBER,
		   real break_r2 = VERY_LARGE_NUMBER,
		   real end_on_unpert = false,
		   real dt_log = 0,
		   real dt_energy = 0,
		   real dt_snap = 0);

real integrate_multiple(hdyn *b);

// ----- Handy clean-up functions (in files of ~same name): -----

void clean_up_hdyn_schedule();
void clean_up_hdyn_ev();
void clean_up_hdyn_grape(unsigned int config);
void clean_up_hdyn_grape4();
void clean_up_hdyn_grape6();
void clean_up_kira_ev();

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/hdyn.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
