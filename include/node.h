
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/*
 *  node.h: base class for nbody systems
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class node
 *.............................................................................
 */

#ifndef  STARLAB_NODE_H
#  define  STARLAB_NODE_H

#include  "starlab_vector.h"
#include  "story.h"
#include  "hydrobase.h"
#include  "starbase.h"

#define __VALID_NODE__		123456789
#define __INVALID_NODE__	-1

/*-----------------------------------------------------------------------------
 *  node  --  a base class for the nodes in a tree of dynamical particles
 *-----------------------------------------------------------------------------
 */
class node
    {
    protected:             // not private, since other classes, such as dyn
	                   // and hdyn should have access to these data.

        // Global variable!

    	static node* root;	  // Address of root node.


	long int node_flag;	  // Experimental attempt to check for
				  // and avoid deleted nodes (Steve, 8/98).
				  // Set only in node constructor.
				  // Reset only in node destructor.

	int  index;               // nodes can be numbered,
	char * name;              // or they can receive individual names.

	real  mass;               // equally accessible by dyn, hydro, star.

	node * parent;            // oya       | As you can see, Japanese is
	node * oldest_daughter;   // choujo    | a more natural language for
	node * elder_sister;      // oneesan   | building a tree structure;
	node * younger_sister;    // imouto    | the names are so much simpler!

        hydrobase * hbase;        // hydrobase is the class underlying all
	                          // hydro classes that handle hydrodynamics.
	                          // This is the only hydro class provided 
	                          // here, in order to allow separate
	                          // recompilation of the inherited classes
				  // based on hydro without necessitating
	                          // compilation of the dynamics part.

        starbase * sbase;         // starbase is the class underlying all star
	                          // classes that handle stellar evolution.
	                          // This is the only star class provided 
	                          // here, in order to allow separate
	                          // recompilation of the inherited classes
				  // based on star without necessitating
	                          // compilation of the dynamics part.

	story * log_story;        // the log story is a generalized scratchpad

	story * dyn_story;        // the dyn story is a placeholder for
                                  // dynamical information not recognized by
                                  // a program -- this allows the information
                                  // to be preserved and passed down a pipe

    public:

	inline void clear_node() {
	    if (name) delete [] name;
	    parent = oldest_daughter = elder_sister = younger_sister = NULL;
	}

	inline void node_init() {
	    node_flag = __VALID_NODE__;
	    index = -1;     // < 0 to show that no number has yet been assigned
	    mass = 1;       // to allow star and hydro to rescale this value
	    name = NULL;
	    clear_node();
	}

	inline void node_set_stories(bool use_stories) {

	    // Potentially very dangerous to set any of the following NULL!

	    if (use_stories) {
		log_story = mk_story_chapter(LOG_ID);
		dyn_story = mk_story_chapter(DYNAMICS_ID);
	    } else {
		log_story = NULL;
		dyn_story = NULL;
	    }
	}

	inline void node_set_hbase(hbpfp the_hbpfp) {
	    if (the_hbpfp)
		hbase = (*the_hbpfp)();
	    else
		hbase = NULL;
	}

	inline void node_set_sbase(sbpfp the_sbpfp) {
	    if (the_sbpfp) {
		sbase = (*the_sbpfp)();
		sbase->set_node(this);
	    } else
		sbase = NULL;
	}

	node(hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true) {
	    node_init();
	    node_set_stories(use_stories);
	    node_set_hbase(the_hbpfp);
	    node_set_sbase(the_sbpfp);
	}

	inline void rmstory() {
	    if (log_story) {
		delete log_story;
		log_story = NULL;
	    }
	    if (dyn_story) {
		delete dyn_story;
		dyn_story = NULL;
	    }
	}

	inline void rmhydrobase() {
	    if (hbase) {
		delete hbase;
		hbase = NULL;
	    }
	}

	inline void rmstarbase() {
	    if (sbase) {
		delete sbase;
		sbase = NULL;
	    }
	}

	virtual ~node() {
	    node_flag = __INVALID_NODE__;
	    if (name) delete [] name;
	    rmstory();
	    rmhydrobase();
	    rmstarbase();
	    if (this == root) root = NULL;
	}

	inline bool is_valid()
	    {return (node_flag == __VALID_NODE__);}
	inline void set_invalid()
	    {node_flag = __INVALID_NODE__;}

	void  set_label(int number)             // to set the index to a number
	    {index = number;}
	void  set_label(char * a_string)        // to set the name to a string
	    {
	    if(name != NULL)
	        delete [] name;
	    name = new char[strlen(a_string)+1];
	    strcpy(name, a_string);
	    }

        // Alternate names for the same functions:

	void  set_index(int number)             // to set the index to a number
	    {index = number;}

	void  set_name(char * a_string)        // to set the name to a string
	{
	    if (name)
	        delete [] name;
	    if (!a_string)
		name = NULL;
	    else {
		name = new char[strlen(a_string)+1];
		strcpy(name, a_string);
	    }
	}

	void clear_name()		    {set_name(NULL);}
	void clear_label()		    {clear_name();}

        void  set_mass(const real new_mass) {mass = new_mass;}

	void  set_parent(node * b)          {parent = b;}
	void  set_oldest_daughter(node * b) {oldest_daughter = b;}
	void  set_elder_sister(node * b)    {elder_sister = b;}
	void  set_younger_sister(node * b)  {younger_sister = b;}
#if 0
	void  set_log_story(story * s)      {log_story = s;}	
#else
	// note that some functions (e.g. node::scan_log_story)
	// simply set the story, and don't bother deleting
	// the old one, which is normally made with new_dyn()
	void  set_log_story(story * s) {
	    if (log_story != NULL) delete log_story;
	    log_story = s;
	}	
#endif
	void  set_dyn_story(story * s)      {dyn_story = s;}

        void  log_comment(char *);
        void  log_history(int, char **);

	void  inc_mass(const real d_mass)          {mass += d_mass;}
	void  scale_mass(const real scale_factor)  {mass *= scale_factor;}

	int  get_index()                    {return index;}
	char * get_name()                   {return name;}

	inline real  get_mass()             {return mass;}

	inline node * get_parent()          {return parent;}
	inline node * get_oldest_daughter() {return oldest_daughter;}
	inline node * get_younger_sister()  {return younger_sister;}
	inline node * get_elder_sister()    {return elder_sister;}

	void set_root(node * b)		    {root = b;}
	node * get_root();

	// Potential time sink: define here and make inline!

	// node * get_top_level_node();

	inline node* get_top_level_node() {

	    if (parent == NULL) return NULL;	// Root node

	    node* n = this;
	    node* g = parent->get_parent();

	    while (g) {
		n = n->get_parent();
		g = g->get_parent();
	    }
	    return n;
	}

	node * get_binary_sister();

	void  set_hydrobase(hydrobase * hb) {hbase = hb;}
	void  set_starbase(starbase * sb)  {sbase = sb;}

	hydrobase * get_hydrobase()         {return hbase;}
	starbase * get_starbase()           {return sbase;}

	story * get_log_story()             {return log_story;}
	story * get_dyn_story() const       {return dyn_story;}
	story * get_hydro_story()           {return hbase->get_hydro_story();}
	story * get_star_story()            {return sbase->get_star_story();}

	virtual void null_pointers();
	virtual void print_static(ostream &s = cerr);

	istream& scan_log_story(istream&, char *);
	istream& scan_hydro_story(istream&);
	virtual istream& scan_star_story(istream&, int level = 0);
	virtual istream& scan_dyn_story(istream&);

	virtual bool check_and_correct_node(bool verbose = false);

	ostream& print_log_story(ostream &s = cout);
	ostream& print_hydro_story(ostream &s = cout);
	ostream& print_star_story(ostream &s = cout,
				  int short_output = 0);
	virtual ostream& print_dyn_story(ostream &s = cout,
					 bool print_xreal = true,
					 int short_output = 0);

	inline bool is_isolated()
	    {return (parent == NULL && oldest_daughter==NULL);}
        inline bool is_root()
	    {return (parent == NULL);}
	inline bool is_leaf()
	    {return (parent != NULL && oldest_daughter==NULL);}

	inline bool is_top_level_node()	    {if (!root) get_root();
					     if (parent != root) return false;
					     return true;
					    }
	inline bool is_top_level_leaf()	    {if (!root) get_root();
					     if (parent != root) return false;
					     return (oldest_daughter == NULL);
					    }
	inline bool is_low_level_node()	    {if (!root) get_root();
					     if (parent == root) return false;
					     return true;
					    }
	inline bool is_low_level_leaf()	    {if (!root) get_root();
					     if (parent == root) return false;
					     return (oldest_daughter == NULL);
					    }

	inline bool is_parent()
	    {return oldest_daughter != NULL;}

        bool is_grandparent();

        node* next_node(node*);
        node* orig_next_node(node*);

        bool name_is(char*);
	char* format_label();

	void print_label(ostream&);
	void pretty_print_node(ostream& s = cerr);
	void pretty_print_tree(ostream& s = cerr, int level = 0);

	void interpret_line(char*, char*);
	int  n_leaves();
	int  n_daughters();

	// SeBa counter access functions.
        //seba_counters* get_seba_counters() {
	//     return starbase->get_seba_counters();}
        //void set_seba_counters(seba_counters *sb) {
	//     starbase->set_seba_counters(sb);}

};

typedef  node *(*npfp)(hbpfp, sbpfp, bool);

typedef node * nodeptr;  // to enable dynamic array declarations such as
                         //    node** node_ptr_array = new nodeptr[n];
                         // (note that the following expression is illegal:
                         //    node** node_ptr_array = new (node *)[n];)

inline  node * new_node(hbpfp the_hbpfp,
			sbpfp the_sbpfp ,
			bool use_stories)
    {return  new node(the_hbpfp, the_sbpfp, use_stories);}

node * mk_flat_tree(int, npfp, hbpfp, sbpfp, bool use_stories = true);
inline node * mknode(int n, hbpfp the_hbpfp = new_hydrobase,
	                    sbpfp the_sbpfp = new_starbase,
		     	    bool use_stories = true)
    {return  mk_flat_tree(n, new_node, the_hbpfp, the_sbpfp, use_stories);}

// Modified first argument (default) 5/03 (Steve):

node * get_node(istream &s = cin,
		npfp the_npfp = new_node,
		hbpfp the_hbpfp = new_hydrobase,
		sbpfp the_sbpfp = new_starbase,
		bool use_stories = true);

// Modified arguments 5/03 (Steve):
//
//  void put_node(ostream &s, node &b,
//  	      bool print_xreal = true,
//  	      int short_output = 0);
//  void put_single_node(ostream &s, node &b,
//  		     bool print_xreal = true,
//  		     int short_output = 0);

void put_node(node *b,
	      ostream &s = cout,
	      bool print_xreal = true,
	      int short_output = 0);
void put_single_node(node *b,
		     ostream &s = cout,
		     bool print_xreal = true,
		     int short_output = 0);

bool forget_node(istream &s = cin);

bool node_contains(node * b, int i);
bool node_contains(node * b, char* s);
bool clump_contains(node * b, int i);
bool clump_contains(node * b, char *s);

void pp(node *, ostream & s = cerr);
void pp2(node *, ostream & s = cerr, int level = 0);

#define  for_all_daughters(dyntype, mother, daughter_name)                    \
         for (dyntype* daughter_name = mother->get_oldest_daughter();         \
	      daughter_name != NULL;                                          \
	      daughter_name = daughter_name->get_younger_sister())

// Note: for_all_nodes and for_all_leaves INCLUDE the base node.

#define  for_all_nodes(dyntype, base, node_name)                              \
         for (dyntype* node_name = base;                                      \
	      node_name != NULL;                                              \
	      node_name = (dyntype*) node_name->next_node(base))

#define  for_all_leaves(dyntype, base, node_name)                             \
         for (dyntype* node_name = base;                                      \
	      node_name != NULL;                                              \
	      node_name = (dyntype*) node_name->next_node(base))              \
	      if (node_name->get_oldest_daughter() == NULL)

node * mknode_mass(int n, real m = 1.0);

// Declaration of functions defined in node_tt.C

real total_mass(node *);
void rmtree(node *b, bool delete_b = true);

void detach_node_from_general_tree(node *n);
void remove_node_with_one_daughter(node *n);
void detach_node_from_binary_tree(node *n);
void extend_tree(node *, node *);
void add_node(node *n, node *parent);
void add_node_before(node *n, node *m);
void insert_node_into_binary_tree(node *, node *, node *);

int is_descendent_of(node *, node *, int);
node * common_ancestor(node *, node *);
node * node_with_index(int i, node * top = NULL);
node * node_with_name(char* s, node * top = NULL);
int depth_of_node(node *);
char * construct_binary_label(node * ni, node * nj);
void label_binary_node(node*);
void label_merger_node(node*);

void print_normal_form(node*, ostream&);
char* get_normal_form(node*);

// From node/util:

void renumber(node* b, int istart, bool mass_order, bool name_nodes = false,
	      bool single_number = false);
void construct_node_name(node* b);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/node.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

