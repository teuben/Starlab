
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file chydro.h  Derived class for hydro systems, with core and envelope.
//
//  version 1:  Jan 1993   Piet Hut
//  version 2:
//
//  This file includes:
//  1) definition of class chydro

#ifndef  STARLAB_CHYDRO_H
#  define  STARLAB_CHYDRO_H

#include  "starlab_vector.h"
#include  "story.h"
#include  "hydrobase.h"
#include  "hydro.h"

/// \a chydro:  The standard class for hydrodynamics, with core and envelope.

class  chydro : public hydro
    {
    protected:

	real  core_mass;
	real  core_radius;

    public:

	chydro(real r = 0, real core_m = 0, real core_r = 0,
	       real mf = 1, real rf = 1, real tf = 1) : hydro(r, mf, rf, tf)
	    {
	    core_mass = core_m;
	    core_radius = core_r;
	    }

    	real  get_core_mass()               {return core_mass;}
	real  get_core_radius()             {return core_radius;}

	void  set_core_mass(const real m)        {core_mass = m;}
	void  set_core_radius(const real r)      {core_radius=r;}

	virtual ostream & print_hydro_story(ostream&);
	virtual istream & scan_hydro_story(istream&);
    };

inline  hydrobase * new_chydro()    {return (hydrobase *) new chydro;}

// Shorthand for conversion from node pointer to chydro pointer:

#define N_CH_PTR ((chydro *)n->get_hydrobase())

inline  real get_core_radius(node * n)
    {return N_CH_PTR->get_r_conv_hydro_to_dyn() * N_CH_PTR->get_core_radius();}
    
inline  void set_core_radius(node * n, const real r)
    {N_CH_PTR->set_core_radius( r / N_CH_PTR->get_r_conv_hydro_to_dyn());}

inline  real get_core_mass(node * n)
    {return N_CH_PTR->get_m_conv_hydro_to_dyn() * N_CH_PTR->get_core_mass();}
    
inline  void set_core_mass(node * n, const real m)
    {N_CH_PTR->set_core_mass( m / N_CH_PTR->get_m_conv_hydro_to_dyn());}
    
void  addchydro(node *, real, real, real);

#endif
 
