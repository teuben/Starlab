
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  dyn.h: base class for nbody systems
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class dyn
 *.............................................................................
 */

#ifndef  STARLAB_DYN_H
#  define  STARLAB_DYN_H

#include <vector>
#include  "util_math.h"
#include  "node.h"
#include  "kepler.h"

#define BIN_INDENT	21	// for print_pert() and print_binary_params()

static const real rnull = 0.0;

//#include  "dyn_kepler.h"   NOTE: this file is included at the end instead,
//                           because it needs to know about the dyn class

/*-----------------------------------------------------------------------------
 *  dyn  --  the simplest class of dynamical particles
 *-----------------------------------------------------------------------------
 */
class  dyn : public node
{
    protected:

        // Global variables:

    	static xreal system_time;	// extended-precision time
    	static real real_system_time;	// convenient real copy

    	static bool use_sstar;  	// activate single star evolution

	// Flag to neglect all internal interactions:

	static bool ignore_internal;

	//------------------------------------------------------------

	// External field specification (moved from hdyn to dyn
	// by Steve, 7/01):

	static unsigned int external_field;
					// 0 ==> none
					// 1 ==> tidal (bit 0)
					// 2 ==> Plummer (bit 1)
					// 3 ==> Plummer+tidal (bits 0+1)
					// 4 ==> power-law (bit 2)
					// etc.

	// Tidal field specifications (changed by Steve, 6/99):
	// Note that alpha and omega are in general independent.

	static int  tidal_type;		// 0 ==> no tidal field
					// 1 ==> point-mass field
					// 2 ==> isothermal halo field
					// 3 ==> disk field (Oort constants)

        static real omega;		// system (circular) angular speed
        static real omega_sq;		// omega squared (probably not needed)

	static real alpha1;		// tidal field is conventionally taken
	static real alpha3;		// to be (-alpha1*x, 0, -alpha3*z)

	static vec tidal_center;	// (fixed) point relative to which
					// tidal forces are computed (7/01)

	// Confining Plummer model parameters (Steve, 7/01):

	static real p_mass;		// total mass
	static real p_scale_sq;		// scale radius squared

	static vec p_center;		// (fixed) point relative to which
					// Plummer forces are computed
	static bool p_friction;		// flag to include dynamical friction
					// (currently for Plummer only)

	// Confining power-law model parameters (Steve, 7/01).
	// New implementation (Steve, 2/04) means that Plummer is
	// no longer a special case (with exponent = 0).

	static real pl_coeff;		// overall coefficient
	static real pl_scale;		// scale/cutoff radius
	static real pl_exponent;	// power-law exponent

	static vec pl_center;		// (fixed) point relative to which
					// power-law forces are computed

	static FILE *ifp, *ofp;  // if not NULL, use this instead of C++ streams

	static bool col_output;  // if true, output in col format

	//------------------------------------------------------------

	vec  pos;	         // position (3-D Cartesian vector)
	vec  vel;	         // velocity: (d/dt) pos
	vec  acc;	         // acceleration: (d/dt) vel
	kepler * kep;		 // actually a pointer to a kepler orbit object

    public:

	static FILE* set_ifp(FILE* p)			{return ifp = p;}
	static FILE* get_ifp()				{return ifp;}
	void clear_ifp()				{ifp = NULL;}
	static FILE* set_ofp(FILE* p)			{return ofp = p;}
	static FILE* get_ofp()				{return ofp;}
	void clear_ofp()				{ofp = NULL;}

	static bool set_col_output(bool b = true)	{return col_output = b;}
	static bool get_col_output()			{return col_output;}

	inline void dyn_init() {
	    pos = vel = acc = 0.0;
	    kep = NULL;
	}

        dyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase,
	    bool use_stories = true)
	   : node(the_hbpfp, the_sbpfp, use_stories)	{dyn_init();}
 
	inline void rmkepler() {
	    if (kep) {
		dyn* s = get_binary_sister();

		if (s && s->kep == kep)
		    s->kep = NULL;

		// changed (unsigned int) to unsigned long
		if ((unsigned long) kep != 1 && (unsigned long) kep != 2)
		    delete kep;		// in case kep is used as a flag

		kep = NULL;
	    }
	}

	virtual ~dyn() {

	    // Note added by Steve 7/7/98:
	    //
	    // Some care is required when deleting the kepler structure,
	    // as binary components in kira share a kepler.  In that case,
	    // when the kepler of the first component is deleted, the
	    // kepler pointer of the other component must also be set
	    // NULL.  Otherwise, an attempt will be made to delete a
	    // nonexistent structure.  Do this via another function so
	    // that the same procedure can be used from other classes.

	    rmkepler();
	}

	inline xreal get_system_time()		{return system_time;}
	inline real get_real_system_time()	{return real_system_time;}
	void set_system_time(xreal t)		{system_time = t;
					         real_system_time = t;}

	void set_ignore_internal(bool i = true)	{ignore_internal = i;}
	bool get_ignore_internal()		{return ignore_internal;}

	//-----------------------------------------------------------------
	// External field:
  
	inline unsigned int get_external_field()       {return external_field;}
	inline void set_external_field(unsigned int e) {external_field = e;}

	// Tidal field (external field #1, bit 0):

	inline void set_tidal_field(int t)	{tidal_type = t;
					         if (t > 0)
						     SETBIT(external_field, 0);
					         else
					     	     CLRBIT(external_field, 0);}
	inline int get_tidal_field()		{return tidal_type;}
	inline void unset_tidal_field(int t)	{set_tidal_field(0);}

        void set_omega(real o)			{omega = o;
						 omega_sq = o*o;}
        inline real get_omega()			{return omega;}
        inline real get_omega_sq()		{return omega_sq;}

	void set_alpha(real a1, real a3)	{alpha1 = a1;
						 alpha3 = a3;}
	inline real get_alpha1()		{return alpha1;}
	inline real get_alpha3()		{return alpha3;}
	inline void set_tidal_center(vec c)	{tidal_center = c;}
	inline vec get_tidal_center()	{return tidal_center;}

	void set_tidal_field(real alpha1, real alpha3,
			     int type = 0, real omega = 0, vec c = 0.0) {
	    if (type > 0) set_tidal_field(type);
	    set_alpha(alpha1, alpha3);
	    set_tidal_center(c);
	}

	// Plummer field (external field #2, bit 1):

	inline void set_plummer()		{SETBIT(external_field, 1);}
	inline bool get_plummer()		{return (GETBIT(external_field,
								1) > 0);}
	inline void unset_plummer()		{CLRBIT(external_field, 1);}

	inline real get_p_mass()		{return p_mass;}
	inline void set_p_mass(real m)		{p_mass = m;}
	inline real get_p_scale_sq()		{return p_scale_sq;}
	inline void set_p_scale_sq(real r2)	{p_scale_sq = r2;}
	inline void set_p_center(vec c) 	{p_center = c;}
	inline vec  get_p_center()		{return p_center;}
	inline bool get_p_friction()		{return p_friction;}
	inline void set_p_friction(bool f)	{p_friction = f;}

	inline void set_plummer(real m, real r2,
				vec c = 0.0, bool f = false) {  // all in one
	    set_plummer();
	    p_mass = m;
	    p_scale_sq = r2;
	    p_center = c;
	    p_friction = f;
	}

	// Power-law field (external field #3, bit 2):

	inline void set_pl()			{SETBIT(external_field, 2);}
	inline bool get_pl()			{return (GETBIT(external_field,
								2) > 0);}
	inline void unset_pl()			{CLRBIT(external_field, 2);}

	inline real get_pl_coeff()		{return pl_coeff;}
	inline void set_pl_coeff(real c)	{pl_coeff = c;}
	inline real get_pl_scale()		{return pl_scale;}
	inline void set_pl_scale(real r)	{pl_scale = r;}
	inline real get_pl_exponent()		{return pl_exponent;}
	inline void set_pl_exponent(real e)	{pl_exponent = e;}
	inline void set_pl_center(vec c) 	{pl_center = c;}
	inline vec get_pl_center()		{return pl_center;}

	inline void set_pl(real A, real a,
			   real e, vec c = 0.0) {	  // all in one

	    set_pl();

	    pl_coeff = A;
	    pl_scale = a;
	    pl_exponent = e;
	    pl_center = c;
	}

	// Generic accessors for (single) external fields (see dyn_external.C):

	real get_external_scale_sq();
	vec get_external_center();

	//-----------------------------------------------------------------

        void  set_pos(const vec& new_pos)      {pos = new_pos;}
        void  set_vel(const vec& new_vel)      {vel = new_vel;}
        void  set_acc(const vec& new_acc)      {acc = new_acc;}

	void  clear_pos()                         {pos = 0.0;}
	void  clear_vel()                         {vel = 0.0;}
	void  clear_acc()                         {acc = 0.0;}

	inline void  inc_pos(const vec& d_pos) {pos += d_pos;}
	inline void  inc_vel(const vec& d_vel) {vel += d_vel;}
	inline void  inc_acc(const vec& d_acc) {acc += d_acc;}

	inline void  scale_pos(const real scale_factor)  {pos *= scale_factor;}
	inline void  scale_vel(const real scale_factor)  {vel *= scale_factor;}
	inline void  scale_acc(const real scale_factor)  {acc *= scale_factor;}

	inline vec  get_pos()                  {return pos;}
	inline vec  get_vel()                  {return vel;}
	inline vec  get_acc()                  {return acc;}

	inline dyn * get_parent()
	    {return (dyn*) node::get_parent();}
	inline dyn * get_oldest_daughter()
	    {return (dyn*)node::get_oldest_daughter();}
	inline dyn * get_younger_sister()
	    {return (dyn*) node::get_younger_sister();}
	inline dyn * get_elder_sister()
	    {return (dyn*) node::get_elder_sister();}

        inline dyn * get_root()
            {return (dyn*) node::get_root();}
        inline dyn * get_top_level_node()
            {return (dyn*) node::get_top_level_node();}
        inline dyn * get_binary_sister()
            {return (dyn*) node::get_binary_sister();}

	void  calculate_acceleration(dyn *, real);

	inline kepler * get_kepler()		    {return kep;}
	void  set_kepler(kepler * new_kep)	    {kep = new_kep;}

	virtual void null_pointers();
	virtual void print_static(ostream &s = cerr);

        virtual istream& scan_dyn_story(istream&);
	virtual bool check_and_correct_node(bool verbose = true);

	virtual ostream& print_dyn_story(ostream &s,
					 bool print_xreal = true,
					 int short_output = 0);

        void  to_com();      // transformation to center-of-mass coordinates
	void set_com(vec r = 0.0, vec v = 0.0);
	void offset_com();
	void reset_com();
        int flatten_node();

	// Declaration of member function defined in dyn_tt.C

	virtual real get_radius();
	virtual void set_radius(real) {};

        bool get_use_sstar()			{return use_sstar;}
	void set_use_sstar(bool u)		{use_sstar = u;}

	// Placeholders (~null in dyn, realized in hdyn)

	virtual bool nn_stats(real energy_cutoff, real kT,
			      vec center, bool verbose,
			      bool long_binary_output = true,
			      int which = 0);
	virtual real print_pert(bool long_binary_output = true,
				int indent = BIN_INDENT);
};

typedef dyn * dynptr;  // to enable dynamic array declarations such as
                       //    dyn** dyn_ptr_array = new dynptr[n];
                       // (note that the following expression is illegal:
                       //    dyn** dyn_ptr_array = new (dyn *)[n];)

typedef vec (dyn::*dyn_VMF_ptr)(void);     // vec member function pointer
typedef void (dyn::*dyn_MF_ptr)(const vec &);     // member function pointer

inline  node * new_dyn(hbpfp the_hbpfp,
		       sbpfp the_sbpfp,
		       bool use_stories)
    {return (node *) new dyn(the_hbpfp, the_sbpfp, use_stories);}

inline  dyn * mkdyn(int n, hbpfp the_hbpfp = new_hydrobase,
	                   sbpfp the_sbpfp = new_starbase)
    {return  (dyn *) mk_flat_tree(n, new_dyn, the_hbpfp, the_sbpfp);}

dyn *get_col(istream& s = cin,
	     npfp the_npfp = new_dyn,
	     hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true);

void put_col(dyn*, ostream& s = cout);

dyn *get_dyn(istream & s = cin,
	     hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true);

inline void put_dyn(dyn *b,			// note: not put_node, now!!
		    ostream &s = cout,
		    bool print_xreal = true,
		    int short_output = 0) {
    return dyn::get_col_output() ? put_col(b, s) :
	put_node(b, s, print_xreal, short_output);
}

dyn* fget_dyn(FILE* fp = dyn::get_ifp()?:stdin);

inline dyn * common_ancestor(dyn * bi, dyn * bj)
    {return (dyn*) common_ancestor((node*)bi, (node*)bj);}

void dbg_message(char*, dyn*);
void dbg_message(char*, dyn*, dyn*);

vec something_relative_to_root(dyn*, dyn_VMF_ptr);

// From dyn/init (used in sdyn3/evolve/bound3.C):

void  makesphere(dyn * root, int n, real R = 1, int u_flag = 0);

// From dyn/util:

real pot_on_general_node(dyn * bj, dyn * bi, real eps2, bool cm);
void calculate_energies(dyn * root, real eps2,
			real & epot, real & ekin, real & etot,
			bool cm = false);
void print_recalculated_energies(dyn *, int, real, real e_corr = 0);

void compute_density(dyn* b,
		     int k = 12,
		     dyn** list = NULL,
		     int n_list = 0);

void merge_low_level_nodes(dyn* b, real frac = 1, int option = 1);

// Special-case functions (superceded):

vector<real>& get_radial_densities(dyn *, vec, vector<real>&,
				   bool (*)(dyn*) = 0);

vector<real>& get_radial_numdensities(dyn *, vec, vector<real>&,
				      bool (*)(dyn*) = 0);

int get_radial_vdisp(dyn *b, vec cpos, vec cvel,
		     int n_zones, real r[], real v2[]);

int get_density_profile(dyn *b, vec cpos,
			int n_zones, real r[], real rho[],
			bool (*select)(dyn*) = NULL);
int get_profile(dyn *b, vec cpos,
		int n_zones, real r[], real q[],
		real (*Q)(dyn*));

// Overloaded functions:

void compute_com(dyn*, vec&, vec&);
void compute_com(dyn*);

void compute_mcom(dyn*, vec&, vec&, real f = 0.9, int n_iter = 2);
void compute_mcom(dyn*, real f = 0.9, int n_iter = 2);

int  get_std_center(dyn*, vec&, vec&, bool verbose = false);
int  get_std_center(dyn*, bool verbose = false);

void compute_max_cod(dyn*, vec&, vec&);
void compute_max_cod(dyn*);

void compute_mean_cod(dyn*, vec&, vec&);
void compute_mean_cod(dyn*);

// From lagrad.C:

void compute_mass_radii(dyn*);
void compute_mass_radii_percentiles(dyn*);
void compute_mass_radii_quartiles(dyn*);

typedef bool boolfn(dyn*);
void compute_general_mass_radii(dyn*, int,
				bool nonlin = false,
				boolfn *bf = NULL);

// From sys_stats.C:

bool parse_sys_stats_main(int argc, char *argv[],
			  int  &which_lagr,
			  bool &binaries,
			  bool &short_output,
			  bool &B_flag,
			  bool &calc_e,
			  bool &n_sq,
			  bool &out,
			  bool &verbose);
void check_addstar(dyn* b);
void sys_stats(dyn* root,
	       real energy_cutoff = 1,
	       bool verbose = true,
	       bool binaries = true,
	       bool long_binary_output = false,
	       int which_lagr = 0,
	       bool print_time = false,
	       bool compute_energy = false,
	       bool allow_n_sq_ops = false,
	       void (*compute_energies)(dyn*, real, real&, real&, real&, bool)
			= calculate_energies,
	       void (*dstar_params)(dyn*) = NULL,
	       bool (*print_dstar_stats)(dyn*, bool, vec, bool) = NULL);

void refine_cluster_mass(dyn *b, int verbose = 0);
void refine_cluster_mass2(dyn *b, int verbose = 0);

// From dyn_stats.C:

real print_binary_params(kepler* k, real m1, real kT,
			 real dist_from_cod,
			 bool verbose = true,
			 bool long_output = true,
			 int init_indent = 0,
			 int indent = BIN_INDENT);

real get_total_energy(dyn* bi, dyn* bj);
real get_period(dyn* bi, dyn* bj);
void get_total_energy_and_period(dyn* bi, dyn* bj, real& E, real& P);

void initialize_kepler_from_dyn_pair(kepler& k, dyn* bi, dyn* bj,
				     bool minimal = false);

void print_binary_from_dyn_pair(dyn* bi, dyn* bj,
				real kT = 0,
				vec center = vec(0,0,0),
				bool verbose = true,
				bool short_output = false);

real print_structure_recursive(dyn* b,
			       void (*dstar_params)(dyn*),
			       int& n_unp, real& e_unp,
			       real kT = 0.0,
			       vec center = vec(0,0,0),
			       bool verbose = true,
			       bool short_output = false,
			       int indent = 0);

real print_structure_recursive(dyn* b,
			       real kT = 0.0,
			       vec center = vec(0,0,0),
			       bool verbose = true,
			       bool short_output = false,
			       int indent = 0);

void compute_core_parameters(dyn*, int, bool, vec&, real&, int&, real&);

// From plot_stars.C:

void plot_stars(dyn * bi,
		int n = 5,
		int k = 3);

// From scale.C:

real get_mass(dyn *b);
void scale_mass(dyn *b, real mscale);
void scale_pos(dyn *b, real rscale, vec com_pos = 0.0);
void scale_vel(dyn *b, real vscale, vec com_vel = 0.0);
real get_top_level_kinetic_energy(dyn *b);
real get_kinetic_energy(dyn *b);
void get_top_level_energies(dyn *b, real eps2, real &potential, real &kinetic);
void scale_virial(dyn *b, real q, real potential, real& kinetic,
		  vec com_vel = 0.0);
real scale_energy(dyn *b, real e, real& energy,
		  vec com_pos = 0.0,
		  vec com_vel = 0.0);
bool parse_scale_main(int argc, char *argv[],
		      real& eps, bool& c_flag,
		      bool& e_flag, real& e,
		      bool& m_flag, real& m,
		      bool& q_flag, real& q,
		      bool& r_flag, real& r,
		      bool& debug);
void scale(dyn *b, real eps,
	   bool c_flag,
	   bool e_flag, real e,
	   bool m_flag, real m,
	   bool q_flag, real q,
	   bool r_flag, real r,
	   bool debug = false,
	   void (*top_level_energies)(dyn*, real, real&, real&)
			= get_top_level_energies);

// From dyn_external.C:

void get_external_acc(dyn * b,
		      vec pos,
		      vec vel,
		      real& pot,
		      vec& acc,
		      vec& jerk,
		      bool pot_only = false);
real get_external_pot(dyn * b,
		      void (*pot_func)(dyn *, real) = NULL);

real vcirc(dyn *b, vec r);

real get_tidal_pot(dyn *b);
real get_plummer_pot(dyn *b);
real get_power_law_pot(dyn *b);
real get_external_virial(dyn * b);
void print_external(unsigned int ext, bool shortbits = false);

void set_friction_beta(real b);
void set_friction_mass(real m);
void set_friction_vel(vec v);
void set_friction_acc(dyn *b, real r);

// From add_plummer.C and add_power_law.C:

bool get_physical_scales(dyn *b, real& mass, real& length, real& time);

void add_plummer(dyn *b,
		 real coeff, real scale,
		 vec center = 0.0,
		 bool n_flag = false,
		 bool verbose = false,
		 bool fric_flag = false);

void toggle_plummer_friction(dyn *b);

void add_power_law(dyn *b,
		   real coeff, real exponent, real scale,
		   vec center = 0.0,
		   bool n_flag = false,
		   bool verbose = false,
		   bool G_flag = false);

// From dyn_story.C:

real get_initial_mass(dyn* b,
		      bool verbose = false,
		      bool mass_set = false,
		      real input_mass = 0);
real get_initial_virial_radius(dyn* b,
			       bool verbose = false,
			       bool r_virial_set = false,
			       real input_r_virial = 0);
real get_initial_jacobi_radius(dyn* b,
			       real r_virial,
			       bool verbose = false,
			       bool r_jacobi_set = false,
			       real input_r_jacobi = 0);
void set_tidal_params(dyn* b,
		      bool verbose,
		      real initial_r_jacobi,
		      real initial_mass,
		      int& tidal_field_type);
void test_tidal_params(dyn* b,
		       bool verbose,
		       real initial_r_jacobi,
		       real initial_r_virial,
		       real initial_mass);
int check_set_tidal(dyn *b, bool verbose = false);
void check_set_plummer(dyn *b, bool verbose = false);
void check_set_power_law(dyn *b, bool verbose = false);
void check_set_external(dyn *b, bool verbose = false, int fric_int = -1);
void check_set_ignore_internal(dyn *b, bool verbose = false);

//----------------------------------------------------------------------
//
// Standard wrappers (for starcluster).
// See dyn/util/wrappers.C for functions not specified inline.

// "System" quantities (b is root node):

int  bound_number(dyn *b);
real bound_mass(dyn *b);
real total_energy(dyn *b);
vec  total_angular_momentum(dyn *b, vec x = 0, vec v = 0);
real core_radius(dyn *b);
real core_mass(dyn *b);
int  core_number(dyn *b);
real virial_radius(dyn *b);
real tidal_radius(dyn *b);
real core_density(dyn *b);

// Quantities defined (mainly) for individual nodes:

inline real	mass(dyn *b)
		{
		    if (b->get_oldest_daughter()) {
			if (b->get_mass() > 0)
			    return b->get_mass();	// assume OK
			else
			    return get_mass(b);	// recompute
		    } else
			return b->get_mass();
		}

inline vec	pos(dyn *b, vec x = 0)	{return b->get_pos()-x;}
inline vec	vel(dyn *b, vec v = 0)	{return b->get_vel()-v;}

inline real	distance(dyn *b, vec x = 0)
					{return abs(b->get_pos()-x);}
inline real	distance_sq(dyn *b, vec x = 0)
					{return square(b->get_pos()-x);}

inline real	speed(dyn *b, vec v = 0)
					{return abs(b->get_vel()-v);}
inline real	speed_sq(dyn *b, vec v = 0)
					{return square(b->get_vel()-v);}

inline real	v_rad_sq(dyn *b, vec x = 0, vec v = 0)
					{
					    vec dx = b->get_pos()-x;
					    vec dv = b->get_vel()-v;
					    real dr2 = square(dx);
					    if (dr2 > 0)
						return square(dx*dv)/dr2;
					    else
						return 0;
					}
inline real	v_rad(dyn *b, vec x = 0, vec v = 0)
					{return sqrt(v_rad_sq(b, x, v));}

inline real	v_trans_sq(dyn *b, vec x = 0, vec v = 0)
					{
					    vec dx = b->get_pos()-x;
					    vec dv = b->get_vel()-v;
					    real dr2 = square(dx);
					    if (dr2 > 0)
						return dv*dv
							 - square(dx*dv)/dr2;
					    else
						return 0;
					}

inline real	v_trans(dyn *b, vec x = 0, vec v = 0)
					{return sqrt(v_trans_sq(b, x, v));}

inline real	energy(dyn *b, vec v = 0)
					{
					    if (find_qmatch(b->get_dyn_story(),
							    "pot")) {
						if (!b->get_parent())	// root
						    return total_energy(b);
						else {
						    real pot =
						      getrq(b->get_dyn_story(),
							    "pot");
						    vec dv = b->get_vel()-v;
						    return b->get_mass()
							     *(pot+0.5*dv*dv);
						}
					    } else
						return 0;
					}

inline vec	angular_momentum(dyn *b, vec x = 0, vec v = 0)
					{
					    if (b->get_parent()) {
						vec dx = b->get_pos()-x;
						vec dv = b->get_vel()-v;
						return b->get_mass()*dx^dv;
					    } else
						return
						    total_angular_momentum(b);
					}

//----------------------------------------------------------------------

#include  "dyn_kepler.h"

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/dyn.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
