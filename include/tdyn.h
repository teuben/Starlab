
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file tdyn.h  Derived class for nbody systems using 4-D trees.
//
//  version 1:  Sep 2000   Steve McMillan
//  version 2:
//
//  This file includes:
//  1) definition of class tdyn

#ifndef  STARLAB_TDYN_H
#  define  STARLAB_TDYN_H

#include "pdyn.h"

/// \a tdyn:  A derived class of dynamical particles.  Contains enough
///           information to reconstruct the N-body system.

class  tdyn : public _pdyn_	// base class is _pdyn_, not pdyn, note...
{
    protected:

	xreal   time;
	vec  jerk;           // (d/dt) acc

	// Locators in the 4D hierarchy:

	tdyn *prev;		///< Previous instance of this node.
	tdyn *next;		///< Next instance of this node.

	bool defunct;		///< Node has been terminated.

    public:

        tdyn(hbpfp the_hbpfp = new_hydrobase, sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true)
	   : _pdyn_(the_hbpfp, the_sbpfp, use_stories) {

	    time = 0;
	    jerk = 0;

	    prev = next = NULL;
	    defunct = false;
	}

	virtual ~tdyn() {}

	inline void set_time(const xreal t)	{time = t;}
	inline xreal get_time()			{return time;}

	inline void clear_jerk()		{jerk = 0;}
        inline void set_jerk(const vec& j)	{jerk = j;}
	inline vec get_jerk()                {return jerk;}

	inline void set_prev(tdyn *p)		{prev = p;}
	inline tdyn * get_prev()		{return prev;}

	inline void set_next(tdyn *p)		{next = p;}
	inline tdyn * get_next()		{return next;}

	inline bool is_defunct()		{return defunct;}

	// Convenient:

	inline tdyn * get_parent()
	    {return (tdyn*) node::get_parent();}
	inline tdyn * get_oldest_daughter()
	    {return (tdyn*)node::get_oldest_daughter();}
	inline tdyn * get_younger_sister()
	    {return (tdyn*) node::get_younger_sister();}
	inline tdyn * get_elder_sister()
	    {return (tdyn*) node::get_elder_sister();}

	/// Set or find the root node pointer.

        inline tdyn * get_root()
            {return (tdyn*) node::get_root();}

	/// Return the top-level node of this node.

        inline tdyn * get_top_level_node()
            {return (tdyn*) node::get_top_level_node();}

	/// Find the binary sister of this node.

        inline tdyn * get_binary_sister()
            {return (tdyn*) node::get_binary_sister();}

	// Necessary virtual functions:

	// virtual istream& scan_star_story(istream&, int level = 0);
        virtual istream& scan_dyn_story(istream&);
	virtual bool check_and_correct_node(bool verbose = true);

        virtual ostream& print_dyn_story(ostream &s,
					 bool print_xreal = true,
					 int short_output = 0);
};

typedef tdyn *tdynptr;	   // to enable dynamic array declarations

inline  node * new_tdyn(hbpfp the_hbpfp, sbpfp the_sbpfp,
			bool use_stories = true)
    {return (node *) new tdyn(the_hbpfp, the_sbpfp, use_stories);}

inline  tdyn * get_tdyn(istream & s = cin,
			hbpfp the_hbpfp = new_hydrobase,
			sbpfp the_sbpfp = new_starbase, bool use_stories = true)
    {return  (tdyn *) get_node(s, new_tdyn, the_hbpfp, the_sbpfp, use_stories);}

#define  put_tdyn  put_node

#endif
