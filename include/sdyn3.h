
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  sdyn3.h: derived class for nbody systems, for scattering experiments
 *.............................................................................
 *    version 1:  Dec 1993   Piet Hut, Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class sdyn3
 *.............................................................................
 */

#ifndef  STARLAB_SDYN3_H
#  define  STARLAB_SDYN3_H

#define DEFAULT_TIDAL_TOL_FACTOR 1e-6

#define SSD_HYSTERESIS_FACTOR 2     // multiplicative safety margin:
                                    // a minimum in ssd only counts as a 
                                    // certified minimum if it is surrounded at
                                    // both sides by maxima that are at least
                                    // this factor larger than the minimum;
                                    // two successive certified minima must
                                    // have at least one intermediate point in
                                    // time at which ssd is larger than either
                                    // minimum by this factor.

#define MAX_N_STEPS	100000000

#include  "_dyn_.h"

/*-----------------------------------------------------------------------------
 *  sdyn3  --  a derived class of dynamical particles, with enough information
 *            to integrate the equations of motions with a 4th-order
 *            generalized Hermite polynomial integrator.
 *-----------------------------------------------------------------------------
 */
class  sdyn3 : public _dyn_
{
    protected:

	xreal  time_offset;
    
	real  energy_dissipation;	   // dissipation associated with this
					   // star due to its merger history

	real  min_encounter_time_sq;       // pair-wise
	real  min_free_fall_time_sq;       // pair-wise

	// In the current implementation of sdyn3_ev, the "nn" quantities
	// are only set if no_ssd_flag = 0 (i.e. on the final iteration).

        sdyn3* nn;               // pointer to the nearest neighbor
        real d_nn_sq;           // distance squared to the nearest neighbor
						
	real  nn_dr2;         // current nearest neighbor distance squared
	int   nn_label;       // identity of current nearest neighbor
	real  min_nn_dr2;     // minimum-ever nearest neighbor distance squared
	int   min_nn_label;   // identity of nearest-ever neighbor

	int   init_nn_label;  // Identity of initial nearest neighbor
	int   nn_change_flag; // While zero, nearest neighbor remains unchanged

	real  new_pot;
	vec  new_pos;
	vec  new_vel;
	vec  new_acc;
	vec  new_jerk;

	void  accumulate_new_acc_and_jerk_from_new(sdyn3 *, real, int, int &);

// the following nine (?) are only used in the root node, but given for all
// particles, to keep the data structure homogeneous.

	real  ssd; // current sum of squares of inter-particle distances

	real  min_ssd;              // minimum sum of squares of all
	                            // inter-particle distances, after the
	                            // certified maximum
	real  max_ssd;              // maximum sum of squares of all
	                            // inter-particle distances, after the
	                            // certified maximum

	real  min_min_ssd;          // minimum-ever sum of squares of all
	                            // inter-particle distances
	int   n_ssd_osc;            // number of oscillations in ssd, so far

	int  ssd_ingoing_flag;      // true iff the last certified extremum was
	                            // a maximum; false iff the last certified
	                            // extremum was a minimum.

// the following three quantities could be put in the dyn_story.

	real n_steps;       // number of integration steps, to make it easy
	                   // to store a partially completed experiment, and 
                           // complete the run later.
	real e_tot_init;       // initial total energy of the system
	real de_tot_abs_max;   // maximum error at any previous step in the
	                       // total energy.

    public:

        sdyn3(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase)
	   : _dyn_(the_hbpfp, the_sbpfp)
	    {
	    energy_dissipation = 0;
	    time_offset = 0;
	    nn = NULL;
	    d_nn_sq = 0;
	    min_nn_dr2 = VERY_LARGE_NUMBER;
	    min_nn_label = -1;   // illegal value; no nearest-ever neighbor set
            init_nn_label = -1;  // illegal value; no initial neighbor set
	    nn_change_flag = 0;
	    n_steps = 0;
	    de_tot_abs_max = 0;
	    max_ssd = min_ssd = min_min_ssd = VERY_LARGE_NUMBER;
	    n_ssd_osc = 0;
	    ssd_ingoing_flag = 1;
	    ssd = VERY_LARGE_NUMBER;
	    }
 
	real  get_pot()                         {return pot;}
	xreal  get_time_offset()                {return time_offset;}

	real  get_new_pot()                     {return new_pot;}
	vec  get_new_pos()                   {return new_pos;}
	vec  get_new_vel()                   {return new_vel;}
	vec  get_new_acc()                   {return new_acc;}
	vec  get_new_jerk()                  {return new_jerk;}

	real  get_min_encounter_time_sq()       {return min_encounter_time_sq;}
	real  get_min_free_fall_time_sq()       {return min_free_fall_time_sq;}

	void  set_min_encounter_time_sq(real t) {min_encounter_time_sq = t;}
	void  set_min_free_fall_time_sq(real t) {min_free_fall_time_sq = t;}

	real  get_n_steps()                     {return n_steps;}
	void  set_n_steps(real n)               {n_steps = n;}

	real  get_e_tot_init()           	{return e_tot_init;}
	void  set_e_tot_init(real en)           {e_tot_init = en;}

	void  set_min_nn_dr2(real r)            {min_nn_dr2 = r;}
	void  set_min_nn_label(int label)       {min_nn_label = label;}
	void  set_nn_dr2(real r)                {nn_dr2 = r;}
	void  set_init_nn_label(int label)      {init_nn_label = label;}
	void  set_nn_change_flag(int flag)	{nn_change_flag = flag;}

	sdyn3 * get_nn(){return nn;}
	void set_nn(sdyn3 * new_nn){nn = new_nn;}
	real get_d_nn_sq(){return d_nn_sq;}
	void set_d_nn_sq(real d){d_nn_sq = d;}
  
	real  get_nn_dr2()                      {return nn_dr2;}
	real  get_min_nn_dr2()                  {return min_nn_dr2;}
	int   get_min_nn_label()                {return min_nn_label;}
	int   get_nn_label()                    {return nn_label;}
	int   get_init_nn_label()               {return init_nn_label;}
	int   get_nn_change_flag()		{return nn_change_flag;}
	real  get_min_min_ssd()                 {return min_min_ssd;}
	int   get_n_ssd_osc()                   {return n_ssd_osc;}

	real  get_ssd()                         {return ssd;}

	real  get_energy_dissipation()		{return energy_dissipation;}
	void  set_energy_dissipation(real d)	{energy_dissipation = d;}

	void  clear_new_interaction()           {new_acc=new_jerk=new_pot = 0;}
	void  clear_de_tot_abs_max()            {de_tot_abs_max = 0;}

        void  prepare_root()                  
	    {
	    min_nn_dr2 = -1;                   // not a valid number
	    }

        void  prepare_branch()
	    {
	    min_min_ssd = -1;                  // not a valid number
	    ssd = -1;                          // not a valid number
	    min_ssd = -1;                      // not a valid number
	    max_ssd = -1;                      // not a valid number
	    }

	void  inc_time(real dt)                 {time += dt;}

	void  begin_offset_time(xreal t_off)	{time -= t_off;
						 time_offset += t_off;}

	void  end_offset_time()			{time += time_offset;
						 time_offset = 0;}

	void  inc_new_pot(real dp)              {new_pot += dp;}
	void  inc_new_acc(vec da)            {new_acc += da;}
	void  inc_new_jerk(vec dj)           {new_jerk += dj;}

	void  inc_ssd(real ds)                  {ssd += ds;}
	    
	void  calculate_new_acc_and_jerk_from_new(sdyn3 *, real, int, int &);

        void  taylor_pred_new_pos_and_vel(const real);

        void  taylor_pred_new_pos(const real dt)
	    {
	    new_pos = pos + dt * ( vel + 0.5*dt * acc );
	    new_pos += (1./6) * dt*dt*dt * jerk;
	    }

        void  taylor_pred_new_vel(const real dt)
	    {
	    new_vel = vel + dt * acc;
	    new_vel += (1./2) * dt*dt * jerk;
	    }

        void  correct_new_acc_and_jerk(const real, const real);
        void  correct_new_pos_and_vel(const real);

        void  store_old_into_new()
	    {
	    new_pot = pot;
	    new_pos = pos;
	    new_vel = vel;
	    new_acc = acc;
	    new_jerk = jerk;
	    }

        void  store_new_into_old()
	    {
	    pot = new_pot;
	    pos = new_pos;
	    vel = new_vel;
	    acc = new_acc;
	    jerk = new_jerk;
	    }

	sdyn3 * get_parent()
	    {return (sdyn3*) node::get_parent();}
	sdyn3 * get_oldest_daughter()
	    {return (sdyn3*)node::get_oldest_daughter();}
	sdyn3 * get_younger_sister()
	    {return (sdyn3*) node::get_younger_sister();}
	sdyn3 * get_elder_sister()
	    {return (sdyn3*) node::get_elder_sister();}

	virtual  istream& scan_dyn_story(istream&);
	virtual  ostream& print_dyn_story(ostream&s,
					  bool print_xreal = true,
					  int short_output = 0);
};

typedef sdyn3 * sdyn3ptr;  // to enable dynamic array declarations such as
                         //    sdyn3** sdyn3_ptr_array = new sdyn3ptr[n];
                         // (note that the following expression is illegal:
                         //    sdyn3** sdyn3_ptr_array = new (sdyn3 *)[n];)

// Third argument below is present because get_node expects it...

inline  node * new_sdyn3(hbpfp the_hbpfp, sbpfp the_sbpfp,
			 bool use_stories = true)
    {return (node *) new sdyn3(the_hbpfp, the_sbpfp);}	 // ignore 3rd argument

inline sdyn3 * mksdyn3(int n, hbpfp the_hbpfp = new_hydrobase,
	                    sbpfp the_sbpfp = new_starbase)
    {return  (sdyn3 *) mk_flat_tree(n, new_sdyn3, the_hbpfp, the_sbpfp);}

inline  sdyn3 * get_sdyn3(istream & s = cin,
			  hbpfp the_hbpfp = new_hydrobase,
			  sbpfp the_sbpfp = new_starbase)
    {return  (sdyn3 *) get_node(s, new_sdyn3, the_hbpfp, the_sbpfp);}
    
#define  put_sdyn3  put_node

typedef  void (*sdyn3_print_fp)(sdyn3ptr);	// to extract information
						// from the integrator

// Useful functions:

void set_kepler_from_sdyn3(kepler&, sdyn3*, sdyn3*);
void kepler_pair_to_triple(kepler&, kepler&, sdyn3*, sdyn3*, sdyn3*);

int system_in_cube(sdyn3*, real);

// Standard integrators for scattering problems:

bool  tree3_evolve(sdyn3*, real, real, real, real, real,
		   real cpu_time_check = 3600,
		   real dt = VERY_LARGE_NUMBER,
		   sdyn3_print_fp p = NULL,
		   int snap_limit = -1);
bool  low_n3_evolve(sdyn3*, real, real, real, real, real, real,
		    int, char *, int, int, real,
		    real cpu_time_check = 3600,
		    real dt = VERY_LARGE_NUMBER,
		    sdyn3_print_fp p = NULL,
		    int snap_limit = -1);

real energy(sdyn3 * b);
real potential_energy(sdyn3 * b);
real total_angular_momentum(sdyn3 * b);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/sdyn3.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

