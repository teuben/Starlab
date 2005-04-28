
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  _dyn_.h:  Derived class for nbody systems using a Hermite integrator
//	      Base class for hdyn, sdyn3, and sdyn.
//.............................................................................
//    version 1:  Aug 1998   Steve McMillan, Simon Portegies Zwart
//    version 2:
//.............................................................................
//     This file includes:
//  1) definition of class _dyn_
//.............................................................................

#ifndef  STARLAB__DYN__H
#  define  STARLAB__DYN__H

#include "dyn.h"
#include "slow_binary.h"

//-----------------------------------------------------------------------------
//  _dyn_  --  a derived class of dynamical particles, with enough
//             information to integrate the equations of motion using
//             a 4th-order Hermite integrator.
//-----------------------------------------------------------------------------

class  _dyn_ : public dyn
{
    protected:

	xreal time;
	real  timestep;

	real pot;		// potential
	vec jerk;		// (d/dt) acc
	vec old_acc;
	vec old_jerk;

	vec k_over_18;		// for 5th-order prediction (GRAPE-6)

	// Relocated radius to try to reduce cache misses in
	// critical functions.  Would like to make it private to force
	// use of the accessor function, since the actual value may need
	// processing before use (e.g. black holes).  However, this
	// apparently also changes the data layout of the class and
	// worsens the caching problem.  For now, discourage use of the
	// raw class data by general functions, and return something
	// usable (e.g. not negative) in the accessor function.
	//						(Steve, 1/05)

	real radius;		// effective (or actual) radius of a node.
	vec pred_pos;		// current predicted pos
	vec pred_vel;		// current predicted vel
	xreal t_pred;		// time corresponding to pred_pos and pred_vel

	slow_binary * slow;	// indicator of "slow" binary motion
				// -- affects all time, prediction, and
				//    correction functions
				// -- actually used in kira, but all relevant
				//    member functions are defined here...

	// Pointer to linked list of slow binary CMs perturbed by this
	// node (also used in kira only):

	slow_perturbed * sp;

	// Note: "time" is always time, and "timestep" is always timestep,
	// regardless of the "slow" settings.  However, dtau is actually
	// used during the integration of slow binary motion.

    public:

        _dyn_(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase)
	   : dyn(the_hbpfp, the_sbpfp, true) {

	    time = timestep = pot = 0;
	    k_over_18 = jerk = old_acc = old_jerk = pred_pos = pred_vel = 0;
	    t_pred = -VERY_LARGE_NUMBER;
	    radius = 0;
	    slow = NULL;
	    sp = NULL;
	}

	virtual ~_dyn_() {

	    // Some care is required when deleting the slow structure,
	    // as binary components share it.  See the note on keplers
	    // in the dyn destructor.

	    if (slow) {
		_dyn_ *s = get_binary_sister();
		if (s) {
		    if (s->slow == slow) s->slow = NULL;
		}
		delete slow;
		slow = NULL;
	    }

	    if (sp) {
		delete sp;	// deletes the entire chain
		sp = NULL;
	    }

	}

	void set_timestep(real dt) {		// dt is always timestep
	    timestep = dt;
	    if (slow) slow->set_dtau(dt/slow->get_kappa());
	}

	inline void set_time(xreal t) {		// t is always time
	    time = t;
	    if (slow) slow->set_tau(slow->time_to_tau(time));
	}

	// The "get_time" accessors should be OK given the above "set"
	// functions.

	inline real get_timestep()	const	{return timestep;}
	inline xreal get_time()		const	{return time;}

	void  clear_interaction()         {acc = jerk = 0;
					   pot = 0;}

	void  set_acc_and_jerk_and_pot(const vec& a,
				       const vec& j, real p)
	   {acc = a; jerk = j; pot = p;}

	inline vec  get_pred_pos()     {predict_loworder(get_system_time());
					   return pred_pos;}
	inline vec  get_pred_vel()     {predict_loworder(get_system_time());
					   return pred_vel;}
	inline vec  get_nopred_pos()	const	{return pred_pos;}
	inline vec  get_nopred_vel()	const	{return pred_vel;}

	inline void set_pred_pos(vec p)	  {pred_pos = p;}
	inline void set_pred_vel(vec v)	  {pred_vel = v;}

	void  clear_t_pred()		  {t_pred = -VERY_LARGE_NUMBER;
					   if (slow) slow->clear_tau_pred();}
	inline xreal get_t_pred()         {return t_pred;}
	inline void  set_t_pred(xreal t)  {t_pred = t;}
	
	real  get_pot()		const	  {return pot;}
	void  set_pot(const real p)	  {pot = p;}
	void  clear_pot()                 {pot = 0;}
	void  inc_pot(const real& d_pot)  {pot += d_pot;}

        inline void set_jerk(const vec& new_jerk)     {jerk = new_jerk;}
        inline void set_old_acc(const vec& new_acc)   {old_acc = new_acc;}
        inline void set_old_jerk(const vec& new_jerk) {old_jerk = new_jerk;}

	inline void clear_jerk()                         {jerk = 0;}

	inline void scale_jerk(const real scale_factor)
	    {jerk *= scale_factor;}
	inline void scale_old_acc(const real scale_factor)
	    {old_acc *= scale_factor;}
	inline void scale_old_jerk(const real scale_factor)
	    {old_jerk *= scale_factor;}

	inline void inc_jerk(const vec& d_jerk)       {jerk += d_jerk; }
	inline void inc_old_acc(const vec& d_acc)     {old_acc += d_acc; }
	inline void inc_old_jerk(const vec& d_jerk)   {old_jerk += d_jerk; }

	inline vec get_jerk()		const	{return jerk;}
	inline vec get_old_acc()	const	{return old_acc;}
	inline vec get_old_jerk()       const	{return old_jerk;}
	inline vec get_k_over_18()	const	{return k_over_18;}

	inline real get_radius()	const	{return abs(radius);}
	void set_radius(real r)			{radius = r;}

	// Slow-binary manipulation functions defined in _dyn_slow.C:

	void create_slow(int k = 1);
	void delete_slow();
	void extend_slow(int k);

	slow_binary* get_slow()		const	{return slow;}
	slow_perturbed* get_sp()	const	{return sp;}

	inline int get_kappa()		const	{if (slow)
						     return slow->get_kappa();
						 else
						     return 1;
					        }

	inline _dyn_ * get_parent() const
	    {return (_dyn_*) node::get_parent();}
	inline _dyn_ * get_oldest_daughter() const
	    {return (_dyn_*)node::get_oldest_daughter();}
	inline _dyn_ * get_younger_sister() const
	    {return (_dyn_*) node::get_younger_sister();}
	inline _dyn_ * get_elder_sister() const
	    {return (_dyn_*) node::get_elder_sister();}

        inline _dyn_ * get_root() const
            {return (_dyn_*) node::get_root();}
        inline _dyn_ * get_top_level_node() const
            {return (_dyn_*) node::get_top_level_node();}
        inline _dyn_ * get_binary_sister()
            {return (_dyn_*) node::get_binary_sister();}

	void  init_pred() {
	    pred_pos = pos;
	    pred_vel = vel;
	    t_pred = get_time();
	    if (slow) slow->init_tau_pred();
	}

	void  store_old_force() {
	    old_acc = acc;
	    old_jerk = jerk;

	    // Save the perturbative forces separately for slow binaries
	    // and their perturbers.  Note that old_acc and old_jerk must
	    // *include* the perturbative pieces.

	    _dyn_ *od = get_oldest_daughter();

	    if (od && od->slow)
		od->slow->store_old_force();

	    if (sp)
		sp->store_old_force();			// store entire chain
	}

	inline void predict_loworder(xreal t) {		// t is always time

	    if (kep) {

		init_pred();

	    } else {

		if (t_pred < t) {

		    real dt;

		    // Note that prediction actually occurs in tau if
		    // slow is set.

		    if (slow)
			dt = slow->time_to_tau(t) - slow->get_tau();
		    else
			dt = t - get_time();

		    // Note that this prediction makes *no* assumptions
		    // about the prediction of binary sisters -- this
		    // may be quite inefficient if no special precautions
		    // are taken when many nodes are predicted.

		    real dt3 = dt * ONE_THIRD;
		    real dt2 = dt * 0.5;
		    pred_pos = ((old_jerk * dt3
				 	+ old_acc) * dt2 + vel) * dt + pos;
		    pred_vel = (old_jerk * dt2
					+ old_acc) * dt + vel;

		    t_pred = t;
		    if (slow)
			slow->set_tau_pred(slow->time_to_tau(t_pred));
		}
	    }
	}

	inline void predict_loworder5(xreal t) {	// t is always time

	    // 5th-order prediction of a top-level node (not checked)

	    if (t_pred < t) {

		real dt = t - get_time();

		real dt3 = dt * ONE_THIRD;
		real dt2 = dt * 0.5;
		pred_pos = (((4.5 * k_over_18 * dt + old_jerk) * dt3
			     	+ old_acc) * dt2 + vel) * dt + pos;
		pred_vel = ((6 * k_over_18 + old_jerk) * dt2
			    	+ old_acc) * dt + vel;

		t_pred = t;
	    }
	}

	inline void predict_from_elder_sister() {

	    // For low-level binary nodes only -- NOT checked here.

	    _dyn_ *s = get_elder_sister();

	    real factor = -s->mass / mass;
	    pred_pos = factor * s->pred_pos;
	    pred_vel = factor * s->pred_vel;
	    t_pred = s->t_pred;
	}

	// Handling the slow_perturbed lists (_dyn_slow.C):

	void clear_slow_perturbed() {if (sp) delete sp; sp = NULL;}
	void check_slow_perturbed(bool verbose = false);
	slow_perturbed *find_slow_perturbed(_dyn_ *n, bool verbose = false);
	slow_perturbed *add_slow_perturbed(_dyn_ *n, bool verbose = false);
	bool copy_slow_perturbed(_dyn_ *to,
				 bool overwrite = false,
				 bool verbose = false);
	void remove_slow_perturbed(_dyn_ *n, bool verbose = false);
	void dump_slow_perturbed(char *string = "");
	void print_slow_perturbed();
	int count_slow_perturbed();
	    
	virtual void null_pointers();
	virtual void print_static(ostream &s = cerr);

        virtual  istream& scan_dyn_story(istream&);
        virtual  ostream& print_dyn_story(ostream &s,
					  bool print_xreal = true,
					  int short_output = 0);

	// Potentially slow function -- avoid wherever possible!

        inline xreal get_next_time() const {return get_time()
						 + (xreal)timestep;}

	// Flat-specific member functions defined in _dyn_/evolve/_dyn_ev.C:
	
	void flat_accumulate_acc_and_jerk(_dyn_ *, real);
	void flat_calculate_acc_and_jerk(_dyn_ *, real);
        void flat_set_first_timestep(real, real);
	void flat_update(const real, const real);
        bool flat_correct();
};

// Also defined in _dyn_ev.C:

void predict_loworder_all(_dyn_* b, xreal t);

// In _dyn_slow.C:

bool is_valid_slow(_dyn_ *pert_node);

typedef _dyn_ * _dyn_ptr;  // to enable dynamic array declarations such as
                           //    _dyn_** _dyn__ptr_array = new _dyn_ptr[n];
                           // (note that the following expression is illegal:
                           //    _dyn_** _dyn__ptr_array = new (_dyn_ *)[n];)

// Third argument below is present because get_node expects it...

inline  node * new__dyn_(hbpfp the_hbpfp, sbpfp the_sbpfp,
			 bool use_stories = true)
    {return (node *) new _dyn_(the_hbpfp, the_sbpfp);}	 // ignore 3rd argument

// "True" below means "use stories".

inline  _dyn_ * get__dyn_(istream & s = cin,
			  hbpfp the_hbpfp = new_hydrobase,
			  sbpfp the_sbpfp = new_starbase)
    {return  (_dyn_ *) get_node(s, new__dyn_, the_hbpfp, the_sbpfp, true);}

#define  put__dyn_  put_node

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/_dyn_.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================
