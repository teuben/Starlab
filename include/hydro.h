
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  hydro.h: base class for hydro systems, on the same level as dyn
 *.............................................................................
 *    version 1:  Jan 1993   Piet Hut
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class hydro
 *.............................................................................
 */

#ifndef  STARLAB_HYDRO_H
#  define  STARLAB_HYDRO_H

#include  "starlab_vector.h"
#include  "story.h"
#include  "hydrobase.h"
#include  "node.h"       // needed in functions such as get_effective_radius().
                         // It cannot be included in  hydrobase.h  since 
                         // node.h contains  hydrobase.h  and this would
                         // lead to a vicious circle

/*-----------------------------------------------------------------------------
 *  hydro  --  the simplest class for hydrodynamics
 *-----------------------------------------------------------------------------
 */
class  hydro : public hydrobase
    {
    protected:

	real  effective_radius;

	real  m_conv_hydro_to_dyn;   // mass conversion factor
	real  r_conv_hydro_to_dyn;   // length conversion factor
	real  t_conv_hydro_to_dyn;   // time conversion factor

    public:

	hydro(real r = 0, real mf = 1, real rf = 1, real tf = 1)
	    {
	    effective_radius = r;
	    m_conv_hydro_to_dyn = mf;
	    r_conv_hydro_to_dyn = rf;
	    t_conv_hydro_to_dyn = tf;
	    }

//      to convert an internal hydro mass variable m_hydro into the units used
//	by the dyn part, you can use:
//
//          a_dyn_particle->set_mass( m_hydro * get_m_conv_hydro_to_dyn() );

    	real  get_m_conv_hydro_to_dyn()      {return m_conv_hydro_to_dyn;}
    	real  get_r_conv_hydro_to_dyn()      {return r_conv_hydro_to_dyn;}
    	real  get_t_conv_hydro_to_dyn()      {return t_conv_hydro_to_dyn;}

	real  get_effective_radius()         {return effective_radius;}

    	void  set_m_conv_hydro_to_dyn(const real mf)
	    {m_conv_hydro_to_dyn = mf;}
    	void  set_r_conv_hydro_to_dyn(const real rf)
	    {r_conv_hydro_to_dyn = rf;}
    	void  set_t_conv_hydro_to_dyn(const real tf)
	    {t_conv_hydro_to_dyn = tf;}

	void  set_effective_radius(const real r)  {effective_radius=r;}

	virtual ostream & print_hydro_story(ostream&);
	virtual istream & scan_hydro_story(istream&);
    };

inline  hydrobase * new_hydro()    {return (hydrobase *) new hydro;}

// Shorthand for conversion from node pointer to hydro pointer:

#define N_H_PTR ((hydro *)n->get_hydrobase())

// note: automatic conversion from hydro to dyn scaling

inline  real get_effective_radius(node * n)
    {return N_H_PTR->get_r_conv_hydro_to_dyn() *
	    N_H_PTR->get_effective_radius();}
    
inline  real get_m_conv_hydro_to_dyn(node * n)
    {return N_H_PTR->get_m_conv_hydro_to_dyn();}
inline  real get_r_conv_hydro_to_dyn(node * n)
    {return N_H_PTR->get_r_conv_hydro_to_dyn();}
inline  real get_t_conv_hydro_to_dyn(node * n)
    {return N_H_PTR->get_t_conv_hydro_to_dyn();}

// note: automatic conversion from dyn to hydro scaling

inline  void set_effective_radius(node * n, const real r)
    {N_H_PTR->set_effective_radius( r / N_H_PTR->get_r_conv_hydro_to_dyn());}
    
inline  void set_m_conv_hydro_to_dyn(node * n, const real r)
    {N_H_PTR->set_m_conv_hydro_to_dyn(r);}
inline  void set_r_conv_hydro_to_dyn(node * n, const real r)
    {N_H_PTR->set_r_conv_hydro_to_dyn(r);}
inline  void set_t_conv_hydro_to_dyn(node * n, const real r)
    {N_H_PTR->set_t_conv_hydro_to_dyn(r);}
    
void  addhydro(node *, real);

#endif
 
