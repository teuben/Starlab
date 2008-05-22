
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                       
//=======================================================//              /|\ ~


// kira_MPI_id_manager:  Manage IDs for MPI parallel kira.
//
// Externally visible functions:
//
//	int add_to_MPI_id_array(hdyn *b, bool check = false)
//	bool remove_from_MPI_id_array(hdyn *b, bool check = false)
//	void recompute_MPI_id_array(hdyn *b)
//	hdynptr get_MPI_hdynptr(int id)

//#ifdef USEMPI

#include "hdyn.h"
#include <vector>
#include <algorithm>

#define DEBUG 0

// Maintain a set of nearly contiguous integer IDs for use in the MPI
// parallelization scheme, stored as hdyn member data, and an array
// indexed by these IDs allowing easy conversion from ID to hdyn pointer.

// List is initialized at system startup, and must be updated or
// recomputed whenever the tree structure changes.

// Use static arrays local to this function for now.
// Later may add to the hdyn class...

static vector<hdynptr> *id_array = NULL;

local inline void initialize_id_array(hdyn *b)		// count nodes
{
    static char *func = "initialize_id_array 1";
    if (DEBUG) cerr << "in " << func << endl;
    if (id_array) return;

    int n = 0;
    for_all_nodes(hdyn, b->get_root(), bb) n++;
    n = (int) (1.2*n);					// arbitrary

    id_array = new vector<hdynptr>;
    id_array->reserve(n);
}

local inline void initialize_id_array(int n)		// specify array length
{
    static char *func = "initialize_id_array 2";
    if (DEBUG) cerr << "in " << func << endl;
    if (id_array) return;

    id_array = new vector<hdynptr>;
    id_array->reserve(n);
}

local inline void reinitialize_id_array(hdyn *b)	// count nodes
{
    if (id_array) delete id_array;
    id_array = NULL;
    initialize_id_array(b);
}

local inline void reinitialize_id_array(int n)		// specify array length
{
    if (id_array) delete id_array;
    id_array = NULL;
    initialize_id_array(n);
}

local inline int find_array_id(hdyn *b)	// should return same as b->get_MPI_id()
{
#if 0

    // The dumb way:

    for (int i = 0; i < id_array->size(); i++)
	if (id_array[i] == b) return i;
    return -1;

#else

    // Better:

    vector<hdynptr>::iterator i
	= find(id_array->begin(), id_array->end(), b);

    if (i == id_array->end())
	return -1;
    else
	return i - id_array->begin();

#endif
}


//----------------------------------------------------------------------
//
// Global functions:

bool check_MPI_id_array(hdyn *b)
{
    static char *func = "check_id_array";
    bool ret = true;

    // Check list consistency, and repair if broken
    //
    //	 - check that all nodes have valid IDs and that
    //	   the IDs refer back to the node
    //	 - non-null entries are legal
    //	 - no duplicates
    //	 - node id matches location in list
    //
    // A false return means we should recompute the list.

    cerr << "in " << func << " at time "
	 << b->get_system_time() << "...";

    for_all_nodes(hdyn, b, bi) {
	int id = bi->get_MPI_id();
	if (id < 0 || id >= id_array->size()) {
	    cerr << "ID " << id << " out of range for "
		 << bi->format_label() << endl;
	    ret = false;
	} else if ((*id_array)[id] != bi) {
	    cerr << "ID pointer incorrect for " << bi->format_label() << endl;
	    hdyn *bj = (*id_array)[id];
	    if (bj == NULL)
		cerr << "    (NULL pointer)" << endl;
	    else if (!bj->is_valid())
		cerr << "    (invalid node)" << endl;
	    else
		cerr << "    (points to " << bj->format_label() << ")" << endl;
	    ret = false;
	}
    }

    if (!ret) return ret;

    // Nodes and id_list are consistent.  Check list integrity.

    for (vector<hdynptr>::iterator i = id_array->begin();
	 i < id_array->end(); i++) {

	if (*i != NULL) {

	    hdynptr bi = *i;
	    int id = i - id_array->begin();

	    if (!bi->is_valid()) {

		cerr << "warning: invalid node at "; PRL(id);
		*i = NULL;

	    } else {

		// id_array points to a valid node, but node may not
		// point back if there are duplicates.

		if (bi->get_MPI_id() != id) {

		    // Node doesn't point back to id.

		    cerr << "warning: ID mismatch at "; PRL(id);
		    PRI(4); PRC(bi->format_label()); PRL(bi->get_MPI_id());

		    // We know that bi->get_MPI_id() is OK, so we must be
		    // a duplicate.  Hence just clear this node.

		    *i = NULL;
		}
	    }
	}
    }

    if (ret) cerr << "OK" << endl;
    return ret;
}

int add_to_MPI_id_array(hdyn *b,
			bool check)		// default = false
{
    static char *func = "add_to_id_array";
    if (DEBUG) {
	cerr << "in " << func << " for " << b->format_label() << ", ";
	PRL(b->get_MPI_id());
    }

    // Add b to the array if necessary, and return its index.
    // Note that indices start at 0, so b is id_array[id].

    // Procedure: First look for b on the list.  If it is there, do
    // nothing.  If it is not on the list, then we have to find a place
    // for it.  Look for an empty slot on the current list.  If one is
    // found, use it.  Otherwise, add b to the end of the list.

    // We could possibly accelerate this by keeping track of the number
    // of empty slots (pointer = NULL) on the list, and maybe even the
    // location of the first.

    if (!id_array) initialize_id_array(b);

    if (DEBUG) PRL(id_array->size());

    int id = -1;
    if (check)
	id = find_array_id(b);			// search the array
    else
	id = b->get_MPI_id();			// use the index

    if (DEBUG) PRL(id);

    if (id >= 0) {
	if (DEBUG)
	    cerr << "ID 1 for " << b->format_label() << " = " << id << endl;
	return id;
    }

    // Node b isn't on the list.  Insert it in the first empty slot,
    // or add it to the end.

    id = find_array_id(NULL);

    if (DEBUG) PRL(id);

    if (id < 0) {
	id_array->push_back(b);
	id = id_array->size() - 1;
    } else
	(*id_array)[id] = b;

    b->set_MPI_id(id);

    if (DEBUG) cerr << "ID 2 for " << b->format_label() << " = " << id << endl;
    return id;
}

bool remove_from_MPI_id_array(hdyn *b,
			      bool check)	// default = false
{
    if (!id_array) return false;

    int id = -1;
    if (check)
	id = find_array_id(b);			// search the array
    else
	id = b->get_MPI_id();			// use the index

    bool ret = false;
    if (id >= 0) {

	// Remove id from the list (set hdynptr NULL).

	(*id_array)[id] = NULL;
	ret = true;
    }

    b->set_MPI_id(-1);
    return ret;
}

void recompute_MPI_id_array(hdyn *b)		// rebuild from scratch
{
    b = b->get_root();
    int n = 0;
    reinitialize_id_array(b);
    for_all_nodes(hdyn, b, bb) {		// root node will have ID 0
	bb->set_MPI_id(-1);
	add_to_MPI_id_array(bb);
    }
}

hdynptr get_MPI_hdynptr(int id)
{
    if (!id_array) return NULL;
    if (id < 0 || id >= id_array->size()) return NULL;
    return (*id_array)[id];
}

//#endif
