
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  pdyn.h: Definitions of derived classes _pdyn_ and pdyn, including basic
//	     physical stellar data, for use in applications involving
//           N-body systems using 4-D trees.
//
//.............................................................................
//    version 1:  May 2001   Steve McMillan
//    version 2:
//.............................................................................
//     This file includes:
//  1) definition of class _pdyn_ (physical dyn, pre-tdyn, partiview-dyn?)
//  2) definition of class pdyn (pdyn with extra kepler data)
//.............................................................................

#ifndef  STARLAB_PDYN_H
#  define  STARLAB_PDYN_H

#include "dyn.h"

//-----------------------------------------------------------------------------
//  _pdyn_  --  a derived class of dynamical particles, with additional
//              information on stellar properties.  Serves as a place-
//	        holder for the stellar data used in tdyn applications.
//-----------------------------------------------------------------------------

class  _pdyn_ : public dyn
{
    protected:

	int worldline_index;

	// Stellar stuff:

	int stellar_type;	// (see inc/star/star_support.h for details)
	real temperature;
	real luminosity;

    public:

	inline void _pdyn_init() {
	    worldline_index = -1;
	    stellar_type = -42;
	    temperature = luminosity = 0;
	}

        _pdyn_(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true)
	    : dyn(the_hbpfp, the_sbpfp, use_stories)	{_pdyn_init();}

	virtual ~_pdyn_() {}

	inline void set_worldline_index(int w)	{worldline_index = w;}
	inline int  get_worldline_index()	{return worldline_index;}

	inline void set_temperature(real t)	{temperature = t;}
	inline real get_temperature()		{return temperature;}

	inline void set_luminosity(real t)	{luminosity = t;}
	inline real get_luminosity()		{return luminosity;}

	inline void set_stellar_type(int s)	{stellar_type = s;}
	inline int  get_stellar_type()		{return stellar_type;}

	inline _pdyn_ * get_parent()
	    {return (_pdyn_*) node::get_parent();}
	inline _pdyn_ * get_oldest_daughter()
	    {return (_pdyn_*)node::get_oldest_daughter();}
	inline _pdyn_ * get_younger_sister()
	    {return (_pdyn_*) node::get_younger_sister();}
	inline _pdyn_ * get_elder_sister()
	    {return (_pdyn_*) node::get_elder_sister();}

        inline _pdyn_ * get_root()
            {return (_pdyn_*) node::get_root();}
        inline _pdyn_ * get_top_level_node()
            {return (_pdyn_*) node::get_top_level_node();}
        inline _pdyn_ * get_binary_sister()
            {return (_pdyn_*) node::get_binary_sister();}

	virtual real get_radius();
        virtual  ostream& print_dyn_story(ostream &s,
					  bool print_xreal = true,
					  int short_output = 0);
};

typedef _pdyn_ *_pdyn_ptr;	   // to enable dynamic array declarations

inline  node * new_pdyn_(hbpfp the_hbpfp, sbpfp the_sbpfp,
			 bool use_stories = true)
    {return (node *) new _pdyn_(the_hbpfp, the_sbpfp, use_stories);}

inline  _pdyn_ * get_pdyn_(istream & s = cin,
			   hbpfp the_hbpfp = new_hydrobase,
			   sbpfp the_sbpfp = new_starbase,
			   bool use_stories = true)
    {return  (_pdyn_ *) get_node(s, new_pdyn_, the_hbpfp, the_sbpfp,
				 use_stories);}

#define  put_pdyn_  put_node

// Definition of the pdyn class = _pdyn_ with an extra kepler pointer.

class tdyn;

class  pdyn : public _pdyn_
{
    protected:

	tdyn *kepevent;
	kepler *kep2;
	tdyn *kepevent2;

    public:

	pdyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true)
	    : _pdyn_(the_hbpfp, the_sbpfp, use_stories) {
		kepevent = kepevent2 = NULL;
		kep2 = NULL;
	    }

	inline void rmkepler2() {
	    if (kep2) {
		pdyn* s = (pdyn*)get_binary_sister();

		if (s && s->kep2 == kep2)
		    s->kep2 = NULL;

		delete kep2;
		kep2 = NULL;
	    }
	}

	virtual ~pdyn()
	{
	    if (kep2)
		rmkepler2();	// see note on rmkepler() in dyn.h
	}

	inline kepler * get_kepler2()		{return kep2;}
	void  set_kepler2(kepler * new_kep)	{kep2 = new_kep;}

	inline tdyn *get_kepevent()		{return kepevent;}
	inline void set_kepevent(tdyn *t)	{kepevent = t;}

	inline tdyn *get_kepevent2()		{return kepevent2;}
	inline void set_kepevent2(tdyn *t)	{kepevent2 = t;}

	// The usual...

	inline pdyn * get_parent()
	    {return (pdyn*) node::get_parent();}
	inline pdyn * get_oldest_daughter()
	    {return (pdyn*)node::get_oldest_daughter();}
	inline pdyn * get_younger_sister()
	    {return (pdyn*) node::get_younger_sister();}
	inline pdyn * get_elder_sister()
	    {return (pdyn*) node::get_elder_sister();}

        inline pdyn * get_root()
            {return (pdyn*) node::get_root();}
        inline pdyn * get_top_level_node()
            {return (pdyn*) node::get_top_level_node();}
        inline pdyn * get_binary_sister()
            {return (pdyn*) node::get_binary_sister();}

	// Necessary virtual functions:

        virtual  istream& scan_star_story(istream&, int level = 0);
        virtual  istream& scan_dyn_story(istream&);

#if 0
	virtual ostream& pdyn::print_star_story(ostream& s,
						int short_output = 0)
#endif

};

typedef pdyn *pdynptr;	   // to enable dynamic array declarations

inline  node * newpdyn(hbpfp the_hbpfp, sbpfp the_sbpfp,
			 bool use_stories = true)
    {return (node *) new pdyn(the_hbpfp, the_sbpfp, use_stories);}

inline  pdyn * getpdyn(istream & s = cin,
		       hbpfp the_hbpfp = new_hydrobase,
		       sbpfp the_sbpfp = new_starbase,
		       bool use_stories = true)
    {return  (pdyn *) get_node(s, newpdyn, the_hbpfp, the_sbpfp,
				 use_stories);}

#define  put_pdyn  put_node

#endif
