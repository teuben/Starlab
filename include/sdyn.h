
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  sdyn.h: derived class for nbody systems, for scattering experiments
 *.............................................................................
 *    version 1:  Mar 1994   Piet Hut, Steve McMillan
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class sdyn
 *.............................................................................
 */

#ifndef  STARLAB_SDYN_H
#  define  STARLAB_SDYN_H

#define MAX_N_STEPS	VERY_LARGE_NUMBER
//#define MAX_N_STEPS	1000000

#include  "_dyn_.h"

/*-----------------------------------------------------------------------------
 *  sdyn  --  a derived class of dynamical particles, with enough information
 *            to integrate the equations of motions with a 4th-order
 *            generalized Hermite polynomial integrator.
 *-----------------------------------------------------------------------------
 */
class  sdyn : public _dyn_
{
    protected:

	xreal  time_offset;
    
	real  energy_dissipation;	   // dissipation associated with this
					   // star due to its merger history

	real  min_encounter_time_sq;       // pair-wise
	real  min_free_fall_time_sq;       // pair-wise

        sdyn* nn;               // pointer to the nearest neighbor
        real d_nn_sq;           // distance squared to the nearest neighbor
						
	real  nn_dr2;         // current nearest neighbor distance squared
	int   nn_label;       // identity of current nearest neighbor
	sdyn* nn_ptr;         // pointer to current nearest neighbor
	sdyn* nnn_ptr;        // pointer to current next-nearest neighbor
	real  min_nn_dr2;     // minimum-ever nearest neighbor distance squared
	int   min_nn_label;   // identity of nearest-ever neighbor

	// For use in determination of substructure:

	real  nnn_dr2;	      // distance squared to next-nearest neighbor
	int   nnn_label;      // identity of next nearest neighbor

	int   init_nn_label;  // Identity of initial nearest neighbor
	int   nn_change_flag; // While zero, nearest neighbor remains unchanged

	real  new_pot;
	vec  new_pos;
	vec  new_vel;
	vec  new_acc;
	vec  new_jerk;

	void  accumulate_new_acc_and_jerk_from_new(sdyn *, real, int, int &, real&);

// the following three quantities could be put in the dyn_story.

	real n_steps;       // number of integration steps, to make it easy
	                   // to store a partially completed experiment, and 
                           // complete the run later.
	real e_tot_init;       // initial total energy of the system
	real de_tot_abs_max;   // maximum error at any previous step in the
	                       // total energy.

// for use in determining quarantine status:

	int  temp_quarantine_flag;
	int  quarantine_flag;
	real quarantine_time;
	real quarantine_sma;
	real quarantine_ecc;

    public:

        sdyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase)
	   : _dyn_(the_hbpfp, the_sbpfp)
	    {
	      nn = NULL;
	      d_nn_sq = 0;
	      //	    radius = 0;
	    energy_dissipation = 0;
	    time_offset = 0;
	    min_nn_dr2 = VERY_LARGE_NUMBER;
	    min_nn_label = -1;   // illegal value; no nearest-ever neighbor set
            init_nn_label = -1;  // illegal value; no initial neighbor set
	    nn_change_flag = 0;
	    n_steps = 0;
	    de_tot_abs_max = 0;
	    }
 
	real  get_pot()                         {return pot;}
	xreal  get_time_offset()                 {return time_offset;}

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

	void  set_e_tot_init(real en)           {e_tot_init = en;}

	void  set_min_nn_dr2(real r)            {min_nn_dr2 = r;}
	void  set_min_nn_label(int label)       {min_nn_label = label;}
	void  set_init_nn_label(int label)      {init_nn_label = label;}
	void  set_nn_change_flag(int flag)	{nn_change_flag = flag;}

	sdyn * get_nn(){return nn;}
	void set_nn(sdyn * new_nn){nn = new_nn;}
	real get_d_nn_sq(){return d_nn_sq;}
	void set_d_nn_sq(real d){d_nn_sq = d;}
  
	void  set_nn_dr2(real r)                {nn_dr2 = r;}
	void  set_nnn_dr2(real r)               {nnn_dr2 = r;}
	void  set_nn_label(int label)           {nn_label = label;}
	void  set_nnn_label(int label)          {nnn_label = label;}
	void  set_nn_ptr(sdyn* ptr)             {nn_ptr = ptr;}
	void  set_nnn_ptr(sdyn* ptr)            {nnn_ptr = ptr;}

	real  get_nn_dr2()                      {return nn_dr2;}
	real  get_nnn_dr2()                     {return nnn_dr2;}
	int   get_nn_label()                    {return nn_label;}
	int   get_nnn_label()                   {return nnn_label;}
	sdyn* get_nn_ptr()                      {return nn_ptr;}
	sdyn* get_nnn_ptr()                     {return nnn_ptr;}

	real  get_min_nn_dr2()                  {return min_nn_dr2;}
	int   get_min_nn_label()                {return min_nn_label;}
	int   get_init_nn_label()               {return init_nn_label;}
	int   get_nn_change_flag()		{return nn_change_flag;}

	real  get_radius()			{return radius;}
	void  set_radius(real r)		{radius = r;}

	real  get_energy_dissipation()		{return energy_dissipation;}
	void  set_energy_dissipation(real d)	{energy_dissipation = d;}

        void  set_temp_quarantine_flag(int f)	{temp_quarantine_flag = f;}
        void  set_quarantine_flag(int f)	{quarantine_flag = f;}
        void  set_quarantine_time(real t)	{quarantine_time = t;}
        void  set_quarantine_sma(real a)	{quarantine_sma = a;}
        void  set_quarantine_ecc(real e)	{quarantine_ecc = e;}

        int   get_temp_quarantine_flag()	{return temp_quarantine_flag;}
        int   get_quarantine_flag()		{return quarantine_flag;}
        real  get_quarantine_time()		{return quarantine_time;}
        real  get_quarantine_sma()		{return quarantine_sma;}
        real  get_quarantine_ecc()		{return quarantine_ecc;}

	void  clear_new_interaction()           {new_acc=new_jerk=new_pot = 0;}
	void  clear_de_tot_abs_max()            {de_tot_abs_max = 0;}

        void  prepare_root()                  
	    {
	    min_nn_dr2 = -1;                   // not a valid number
	    }

        void  prepare_branch()
	    {
	    // nothing in here for now (there was ssd info in the sdyn3 case)
	    }

	void  inc_time(real dt)                 {time += dt;}

	void  begin_offset_time(real t_off)	{time -= t_off;
						 time_offset += t_off;}

	void  end_offset_time()			{time += time_offset;
						 time_offset = 0;}

	void  inc_new_pot(real dp)              {new_pot += dp;}
	void  inc_new_acc(vec da)            {new_acc += da;}
	void  inc_new_jerk(vec dj)           {new_jerk += dj;}

	void  calculate_new_acc_and_jerk_from_new(sdyn *, real, int, int &, real&);

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

        void  store_new_into_old()
	    {
	    pot = new_pot;
	    pos = new_pos;
	    vel = new_vel;
	    acc = new_acc;
	    jerk = new_jerk;
	    }

	sdyn * get_parent()
	    {return (sdyn*) node::get_parent();}
	sdyn * get_oldest_daughter()
	    {return (sdyn*)node::get_oldest_daughter();}
	sdyn * get_younger_sister()
	    {return (sdyn*) node::get_younger_sister();}
	sdyn * get_elder_sister()
	    {return (sdyn*) node::get_elder_sister();}

	virtual  istream& scan_dyn_story(istream&);
	virtual  ostream& print_dyn_story(ostream &s,
					  bool print_xreal = true,
					  int short_output = 0);
};

typedef sdyn * sdynptr;  // to enable dynamic array declarations such as
                         //    sdyn** sdyn_ptr_array = new sdynptr[n];
                         // (note that the following expression is illegal:
                         //    sdyn** sdyn_ptr_array = new (sdyn *)[n];)

inline  node * new_sdyn(hbpfp the_hbpfp, sbpfp the_sbpfp,
			bool use_stories = true)
    {return (node *) new sdyn(the_hbpfp, the_sbpfp);}

inline sdyn * mksdyn(int n, hbpfp the_hbpfp = new_hydrobase,
	                    sbpfp the_sbpfp = new_starbase)
    {return  (sdyn *) mk_flat_tree(n, new_sdyn, the_hbpfp, the_sbpfp);}

inline  sdyn * get_sdyn(istream & s = cin,
			hbpfp the_hbpfp = new_hydrobase,
			sbpfp the_sbpfp = new_starbase)
    {return  (sdyn *) get_node(s, new_sdyn, the_hbpfp, the_sbpfp);}
    
#define  put_sdyn  put_node

typedef  void (*sdyn_print_fp)(sdynptr);	// to extract information
						// from the integrator
void make_tree(sdyn*, bool, bool, int, int);

real calculate_energy(sdyn*, real&, real&);
real calculate_energy_from_scratch(sdyn*, real&, real&);

sdyn* first_leaf(sdyn*);
sdyn* next_leaf(sdyn*);

char* id(sdyn*);

// Standard integrators for scattering problems:

bool  tree_evolve(sdyn*, real, real, real, real, real,
		  real &,
		  real cpu_time_check = 3600,
		  real dt = VERY_LARGE_NUMBER,
		  sdyn_print_fp p = NULL);
bool  low_n_evolve(sdyn*, real, real, real, real, real, real,
		   int, char *, int, int, real, 
		   real&,
		   real cpu_time_check = 3600,
		   real dt = VERY_LARGE_NUMBER,
		   sdyn_print_fp p = NULL);

real merge_collisions(sdyn * b, int ci);
void merge(sdyn * bi, sdyn * bj);
bool tree_is_unbound(sdyn* root, real ttf, int debug);
#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/sdyn.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

