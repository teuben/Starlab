
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//  hdyn_grape4.C: functions to use GRAPE-4 (HARP-3)
//
//  **********************************************************
//  ***  Excess Debugging lines removed by Steve 6/16/00.  ***
//  ***  See 000616hhyn_harp3.C                            ***
//  **********************************************************
//.............................................................................
//    version 1:  Aug 1996   Jun Makino
//    version 2:  Aug 1998   Steve McMillan
//.............................................................................
//
//	Externally visible functions (note that GRAPE number is unspecified):
//
//	void check_release_grape
//	void grape_calculate_energies
//	void grape_calculate_acc_and_jerk
//	void grape_calculate_densities
//	void clean_up_hdyn_grape
//
//.............................................................................
//
//  This file includes:
//
//    -     grape_calculate_acc_and_jerk(next_nodes, n_next)
//          A single function which takes care of everything...
//
//    -     check_release_grape()
//          Releases the GRAPE hardware so that someone else can use it
//          after a prescribed amount of CPU time has passed)
//
//    -	    Routines to determine neighbor lists (called from within
//	    grape_calculate_acc_and_jerk)
//
//  When called, it performs the following:
//
// -- In the first call, initialize the GRAPE hardware
//
// -- If the tree structure is changed (or in the first call),
//    send all particles anew and store the number of particles
//
// -- If there are particles in "previous blockstep", which are not in
//    present block step, update these particles as well
//
// -- update the particles in the present block
//
// -- For block of particles,
//    --- Perform force calculation
//    --- store the force
//    --- retrieve the neighbor list
//    --- If no neighbor is found for any of particles, change
//        its neighbor sphere radius and try over again
//    --- In the case of CM particle, store the neighbors to the
//        perturber list
//    --- In the case of single particle, find the nearst neighbor
//    --- and "closest surface" particle. In order to be able to
//        find the closest surface particle, the search radius of
//        must be larger than 2*radius (at least the star with
//        larger radius can find colliding star)
//
// -- store the "previous block" particles

// Class variables related to GRAPE
//
// grape_rnb_sq : neighbour radius (squared) sent to GRAPE hardware
//
// It must be initialized to some value when a particle/node
// is created (or at the first call).  A reasonable value is
// r_min_sq (could be too small, but will be adjusted anyway;
// the details depend on the information required by kira).
//
// - For SINGLE particles, the internal routine should take care
//   of the updating of this variable.
//
// - For CENTER_OF_MASS, there are two possibilities
//   a) if it has a non-null valid perturber list, the perturbation
//	radius can be used.
//   b) if it does not have a non-null valid perturber list, ....
//   c) if it has an overflowed perturber list, I guess you
// 	should still set the perturbation radius as the neighbour
// 	sphere radius...

//.............................................................................
//
// Known Problems
//
// 1) Determination of h2 for unperturbed binary is not quite okay.
//    it should maintain the nearest neighbor? --- Well, seems to
//    do that already for most cases
//
// Question: To what extent should GRAPE be integrated into kira?
//           To have a few extra variables should be no problem.
//
//           One could define a derived class which might make
//           it easier to move to newer (if any...) hardware or
//           different software interface
//
//.............................................................................

// Note from Steve to Steve (7/9/97):
//
// force_by_grape4()	     is called only from grape_calculate_acc_and_jerk()
// force_by_grape4_on_leaves()		"	 force_by_grape4_on_all_leaves()
// force_by_grape4_on_all_leaves()	"	 grape_calculate_energies()

// force_by_grape4 computes forces between TOP-LEVEL NODES in the
// POINT-MASS approximation (like flat_calculate_acc_and_jerk).

// force_by_grape4_on_all_leaves loads all LEAVES into the hardware and
// hence effectively resolves all multiples into components.

//.............................................................................

// Note the extensive use of static arrays within this file, to allow
// management of the GRAPE interface.
//
// Made these global within the file, even though most are only really
// local to specific functions to allow the possibility of clean-up
// (helpful when using ccmalloc to find real memory leaks).
//
//						Steve, 7/99
//.............................................................................

#include "hdyn.h"
#include "grape4.h"
#include "hdyn_inline.C"

// Static declarations:
// -------------------

static bool grape_was_used_to_calculate_potential = false;
static bool grape_is_open = false;
static bool grape_first_attach = true;

static int nboards;

#define MINIMUM_GRAPE_RNB_PERT 1e-20

// Old static memory allocation:

//#   define HARPNMAX 100000
//#   define HARPNBMAX HARPNMAX
//
// static hdyn* nodes[HARPNMAX];
// static hdyn* next_top[HARPNMAX];
// static hdyn* previous_nodes[HARPNMAX];
// static int nb_check_counter[HARPNMAX];
// static int h3nb[HARPNBMAX];



//  **************************************************************************
//  *                                                                        *
//  * Local functions used by more than one other function (global or local) *
//  *                                                                        *
//  **************************************************************************
//
//	jpdma_nodes
//	force_by_grape4
//	put_grape_index_to_top_level_nodes
//	initialize_array

// Static data ("j-arrays"):
//
// Initialized by function initialize_node_lists(), called from
//
//	put_grape_index_to_top_level_nodes()
// and
//	put_grape_index_to_leaf_nodes()

static int grape_n_max = 0;
static int n_previous_nodes = 0;

static hdyn** nodes = NULL;
static hdyn** next_top = NULL;
static hdyn** previous_nodes = NULL;
static int* nb_check_counter = NULL;

static vector * pxj = NULL;
static vector * pvj = NULL;
static vector * paj = NULL;
static vector * pjj = NULL;
static real   * ptj = NULL;
static real   * pmj = NULL;
static int    * ppj = NULL;
static int    * h3nb = NULL;

local void initialize_node_lists()
{
    if (!nodes) {
	nodes = new hdynptr[grape_n_max];
	next_top = new hdynptr[grape_n_max];
	previous_nodes = new hdynptr[grape_n_max];

	// Initialize nb_check_counter[] to zero when first created.

	nb_check_counter = new int [grape_n_max];
	for (int i = 0; i < grape_n_max; i++)
	    nb_check_counter[i] = 0;

	pxj = new vector[grape_n_max];
	pvj = new vector[grape_n_max];
	paj = new vector[grape_n_max];
	pjj = new vector[grape_n_max];
	ptj = new real[grape_n_max];
	pmj = new real[grape_n_max];
	ppj = new int[grape_n_max];
	h3nb = new int[grape_n_max];
    }
}

local void jpdma_nodes(int nnodes, hdyn * nodelist[], bool predicted,
		       real predicted_time)

// Load the j-nodes on the list into the GRAPE.
//
// Called by:	initialize_array()				local
//		grape_calculate_acc_and_jerk() [twice]		global
//		grape_calculate_densities()			global

{
    int jpmax;
    int j, jp;
    hdyn * b;
    jpmax = h3jpmax_();

    jp = 0;
    for (j = 0; j < nnodes; j++) {
	int jj, jm ;
	b = nodelist[j];

	jj =  b->get_grape_index();
	jm = jj-1;
	if (! predicted) {

	    pxj[jm] = b->get_pos();
	    pvj[jm] = b->get_vel();
	    paj[jm] = b->get_acc()*0.5;
	    pjj[jm] = b->get_jerk()*0.1666666666666666666666666666666;
	    ptj[jm] = b->get_time();
	    pmj[jm] = b->get_mass();

	} else {

	    pxj[jm] = b->get_pred_pos();
	    pvj[jm] = b->get_pred_vel();
	    ptj[jm] = predicted_time;
	    pmj[jm] = b->get_mass();

	}
	ppj[jp] = jj;
	
	jp ++;
	int bufid = 0;
	if ((jp == jpmax) || (j == nnodes - 1)) {
	    int zero = 0;
	    int one = 1;
	    h3mjpdma_indirect_(&jp, ppj, pxj, pvj, paj, pjj, pmj, ptj, &one,
			       &bufid);
	    h3mjpdma_start_(&bufid);
	    bufid ++;
	    if (bufid > 5) bufid = 0;  // this 5 should be < JPDMA_BUFFMAX...
	    h3wait_();
	    h3wait_();

	    jp = 0;
	}
    }
}

// More static data:

static vector * px = NULL;
static vector * pv = NULL;
static vector * pa = NULL;
static vector * pj = NULL;
static real   * ppot = NULL;
static real   * peps2 = NULL;
static real   * ph2 = NULL;

local void force_by_grape4(real time, int ni, hdyn * nodes[], int nj)

// Calculate the force on top-level nodes.
//
// Called by:	grape_calculate_acc_and_jerk()			global
//		grape_calculate_densities()			global

{
    int npipe;
    int i, ip,k;
    hdyn * b;
    npipe = h3npipe_();

    if (px == NULL) {
	px = new vector[npipe];
	pv = new vector[npipe];
	pa = new vector[npipe];
	pj = new vector[npipe];
	peps2 = new real[npipe];
	ph2 = new real[npipe];
	ppot = new real[npipe];
	for (i=0; i<npipe; i++) {
	    peps2[i] = ph2[i] = 0.0;	// hmmm... looks like eps = 0 always
	}
    }

#ifdef USE_HALF_CHIP
    if (ni*2 > npipe)
	err_exit("force by grape4: ni too large");
#else
    if (ni > npipe)
	err_exit("force by grape4: ni too large\n");
#endif

    static real ti;
    ti = time;
    int dbg_mode = 0;
    if (dbg_mode) {cerr << "force by grape4 "; PRL(time);}

    h3setti_(&ti);

    int j;
    for (i = 0; i < ni; i++) {
	b = nodes[i];

#ifndef USE_HALF_CHIP
	j = i;
#else	
	j= 2*i;
#endif	

	px[j] = b->get_pred_pos();
	pv[j] = b->get_pred_vel();
	ph2[j] = b->get_grape_rnb_sq();

	if (dbg_mode) {
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << "getting acc on particle " << j << " " << b->format_label()
		 << endl;
	    PRI(8); PRL(px[j]);
	    PRI(8); PRL(pv[j]);
//	    PRI(8); PRL(ph2[j]);
	    cerr.precision(p);
	}

#ifdef USE_HALF_CHIP
	j++;
	px[j] = b->get_pred_pos();
	pv[j] = b->get_pred_vel();
	ph2[j] = b->get_grape_rnb_sq();
#endif	

    }
    for (int jj = j + 1; jj < npipe; jj++) {
	px[jj] = b->get_pred_pos();
	pv[jj] = b->get_pred_vel();
	ph2[jj] = b->get_grape_rnb_sq();
    }

//  h3calc_(&nj, &npipe, px, pv, peps2, ph2, pa, pj, ppot);

    h3calc_firsthalf_(&nj, &npipe, px, pv, peps2, ph2);
    h3wait_();
    h3calc_lasthalf_(&npipe, pa, pj, ppot);

#if 0
    {
	int i;
	for (i = 0; i < npipe; i++) {

	    fprintf(stderr,
	"##### %d %21.14e %21.14e %21.14e %21.14e %21.14e %21.14e %21.14e\n",
		    i,
		    pa[i*3+0], pa[i*3+1], pa[i*3+2],
		    pj[i*3+0], pj[i*3+1], pj[i*3+2],
		    ppot[i]);
	}
    }
#endif

    for (k = 0; k < ni; k++) {

#ifndef USE_HALF_CHIP
	int j = k;
#else
	int j= 2*k;
#endif

	nodes[k]->set_acc_and_jerk_and_pot(pa[j], pj[j], ppot[j]);
    }

    nboards = h3getnboards_();
    h3nbread_(&nboards);
    if (dbg_mode) {
	cerr << "exit force by grape4 "; PRL(time);
    }
}

local int put_grape_index_to_top_level_nodes(hdyn* b)

// Called by:	initialize_array()				local

{
    int index = 0;
    int node_count = 0;
    for_all_daughters(hdyn, b, bb) node_count++;

    if (grape_n_max == 0) {

//	grape_n_max = (int) (node_count * 1.5) + 10;
	grape_n_max = (int) (node_count * 3.0) + 10;

    } else if (grape_n_max < node_count) {

//	grape_n_max = (int) (node_count * 1.5) + 10;
	grape_n_max = (int) (node_count * 3.0) + 10;

	delete nodes;
	delete next_top;
	delete previous_nodes;
	delete nb_check_counter;
	delete pxj;
	delete pvj;
	delete paj;
	delete pjj;
	delete ptj;
	delete pmj;
	delete ppj;
	delete h3nb;
	nodes = NULL;
    }

    if (nodes == NULL)
	initialize_node_lists();

    for_all_daughters(hdyn, b, bbb) {
	nodes[index] = bbb;
	index++;
	bbb->set_grape_index(index);
    }
    return index;
}

local int initialize_array(hdyn * root)

// Called by:	grape_calculate_acc_and_jerk()			global
//		grape_calculate_densities()			global

{
    // Make list of and index top-level nodes.

    int nj = put_grape_index_to_top_level_nodes(root);

    for (int i = 0; i < (nj+10); i++)
	nb_check_counter[i] = 0;

    // Copy the nodes to the GRAPE.

    jpdma_nodes(nj, nodes, false, 0.0);

    // cerr << "init_array "; PRL(nj);

    return nj;
}



//  **********************************************************************
//  *									 *
//  * Globally visible GRAPE-4 functions (and dedicated local helpers).  *
//  * (Probably should be in separate files, but convenient to combine.) *
//  *									 *
//  **********************************************************************

//========================================================================
// check_release_grape:  Accessor for GRAPE release/attach.
//========================================================================

void check_release_grape(kira_options *ko, xreal time)
{
#ifdef SHARE_GRAPE

    // cerr << "GRAPE CPU check: "; PRL(cpu_time());

    if (cpu_time() - ko->grape_last_cpu > ko->grape_max_cpu) {

	cerr << "\nReleasing GRAPE-4 at time " << time << " after "
	     << cpu_time() - ko->grape_last_cpu <<" CPU sec\n";

	h3close_();
	grape_is_open = false;
	grape_first_attach = false;
    }

#endif
}



//======================================================================
// grape_calculate_energies:  Calculate total energy of the system
//			      (requires GRAPE reset after use).
//======================================================================

local int put_grape_index_to_leaf_nodes(hdyn* b,
					bool cm = false) // use CM approximation

// Called by:	jpdma_all_leaves()

{
    if (grape_n_max == 0) {
	int node_count = 0;
	for_all_daughters(hdyn, b, bb) node_count++;
	grape_n_max = (int) (node_count * 3.0) + 10;
    }

    if (nodes == NULL)
	initialize_node_lists();

    // Operation in case cm = true is similar to that of
    // put_grape_index_to_top_level_nodes, but details and
    // data structures may differ, so don't reuse.

    int index = 0;
    if (!cm) {
	for_all_leaves(hdyn, b, bb) {
	    nodes[index] = bb;
	    index++;
	    bb->set_grape_index(index);
	}
    } else {
	for_all_daughters(hdyn, b, bb) {
	    nodes[index] = bb;			// replicated code...
	    index++;
	    bb->set_grape_index(index);
	}
    }
    return index;
}

local inline void jpdma_node(hdyn *b,
			     bool predicted,
			     real predicted_time,
			     int jpmax, int nj, int& jp)

// Called by: jpdma_all_leaves()

{
    int jj, jm ;
    jj =  b->get_grape_index();
    jm = jj-1;

    // Allow the possibility of using already predicted pos and vel
    // rather than having the GRAPE do the prediction.

    if (!predicted) {

	pxj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_pos);
	pvj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_vel);
	paj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_acc)*0.5;
	pjj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_jerk)
       			   * 0.1666666666666666666666666666666;
	ptj[jm] = b->get_time();
	pmj[jm] = b->get_mass();

    } else {

	pxj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_pred_pos);
	pvj[jm] = hdyn_something_relative_to_root(b, &hdyn::get_pred_vel);
	ptj[jm] = predicted_time;
	pmj[jm] = b->get_mass();

    }
    ppj[jp] = jj;
	
    jp++;
    int bufid = 0;
    if ((jp == jpmax) || (jj == nj)) {
	int zero = 0;
	int one = 1;
	h3mjpdma_indirect_(&jp, ppj, pxj, pvj, paj, pjj, pmj, ptj,&one,
			   &bufid);
	h3mjpdma_start_(&bufid);
	bufid++;
	if (bufid > 5)	// this 5 should be < JPDMA_BUFFMAX...
	    bufid = 0;
	h3wait_();
	h3wait_();
	
	jp = 0;
    }
}

local int jpdma_all_leaves(hdyn *root,
			   bool predicted,
			   real predicted_time,
			   bool cm = false)		// use CM approximation

// Called by:	grape_calculate_energies()

{
    int jpmax = h3jpmax_();
    int nj = put_grape_index_to_leaf_nodes(root, cm);
    int jp = 0;

    if (!cm)
	for_all_leaves(hdyn,root,b)
	    jpdma_node(b, predicted, predicted_time, jpmax, nj, jp);
    else
	for_all_daughters(hdyn,root,b)
	    jpdma_node(b, predicted, predicted_time, jpmax, nj, jp);

    return nj;
}

// Yet more static data:

static vector * pxl = NULL;
static vector * pvl = NULL;
static vector * pal = NULL;
static vector * pjl = NULL;
static real   * ppl = NULL;
static real   * peps2l = NULL;
static real   * ph2l = NULL;

local void force_by_grape4_on_leaves(real time, int ni, hdyn * nodes[], int nj)

// Note: All pipelines are used; leaves taken from nodes[].
//
// Called by:	force_by_grape4_on_all_leaves()

{
    int npipe;
    int i, ip,k;
    hdyn * b;
    npipe = h3npipe_();

    if (pxl == NULL) {
	pxl = new vector[npipe];
	pvl = new vector[npipe];
	pal = new vector[npipe];
	pjl = new vector[npipe];
	peps2l = new real[npipe];
	ph2l = new real[npipe];
	ppl = new real[npipe];
	for (i=0; i<npipe; i++) {
	    peps2l[i] = ph2l[i] = 0.0;
	}
    }

    if (ni > npipe)
	err_exit("force by grape4: ni too large\n");

    static real ti;

    ti = time;
    int dbg_mode = 0;
    if (dbg_mode) {cerr << "force by grape4 "; PRL(time);}
    h3setti_(&ti);

    int j;
    for (i = 0; i < ni; i++) {
	b = nodes[i];
	
	j = i;
	pxl[j] = hdyn_something_relative_to_root(b, &hdyn::get_pred_pos);
	pvl[j] = hdyn_something_relative_to_root(b, &hdyn::get_pred_vel);
	ph2l[j] = 0;
	if (dbg_mode) {
	    cerr << "force on  particle "; b->print_label(cerr);
	    PRL(j);
	    PRL(pxl[j]);
	    PRL(pvl[j]);
	    PRL(ph2l[j]);
	}
    }

    for (int jj = j + 1; jj < npipe; jj++) {
	pxl[jj] = b->get_pred_pos();
	pvl[jj] = b->get_pred_vel();
	ph2l[jj] = b->get_grape_rnb_sq();
    }

    h3calc_firsthalf_(&nj, &npipe, pxl, pvl, peps2l, ph2l);
    h3wait_();
    h3calc_lasthalf_(&npipe, pal, pjl, ppl);

    for (k = 0; k < ni; k++) {
	int j = k;

	// We *don't* want to change/set acc and jerk here!

	// nodes[k]->set_acc_and_jerk_and_pot(pal[j], pjl[j], ppl[j]);

	nodes[k]->set_pot(ppl[j]);

    }
    if (dbg_mode) {cerr << "exit force by grape4 "; PRL(time);}
}

// More static data:

static hdyn ** allnodes = NULL;

local void force_by_grape4_on_all_leaves(real time, hdyn * b, int nj,
					 bool cm = false)	// CM approx

// Called by:	grape_calculate_energies()

{
    int npipe = h3npipe_();

    if (!allnodes)
	allnodes = new hdynptr[npipe];

    int i = 0;
    int ip = 0;

    if (!cm) {
	for_all_leaves(hdyn,b,bb) {
	    allnodes[ip] = bb;
	    i++;
	    ip++;
	    if (ip == npipe) {
		force_by_grape4_on_leaves(time, ip, allnodes, nj);
		ip = 0;
	    }
	}
    } else {
	for_all_daughters(hdyn,b,bb) {
	    allnodes[ip] = bb;			// more replicated code...
	    i++;
	    ip++;
	    if (ip == npipe) {
		force_by_grape4_on_leaves(time, ip, allnodes, nj);
		ip = 0;
	    }
	}
    }

    if (ip)
	force_by_grape4_on_leaves(time, ip, allnodes, nj);
}

//  *****************************
//  *****************************
//  ***                       ***
//  ***  The global function  ***
//  ***                       ***
//  *****************************
//  *****************************

void grape_calculate_energies(hdyn * b,
			      real &epot,
			      real &ekin,
			      real &etot,
			      bool cm)		// default = false
{
    if (!grape_is_open) {

	cerr << endl << "grape_calculate_energies: ";
	if (!grape_first_attach) cerr << "re";
	cerr << "attaching GRAPE" << endl;
	h3open_();

	set_time_check_mode(0);
	if (b->get_kira_options())
	    b->get_kira_options()->grape_last_cpu = cpu_time();
	int dbgl = 0;
	h3setdebuglevel_(&dbgl);
	grape_is_open = true;
    }

    real time = b->get_system_time();
    int nj =  jpdma_all_leaves(b, true, time, cm);

    //                            ^^^^  use predicted pos and vel (unpert. OK)

    force_by_grape4_on_all_leaves(time, b,  nj, cm);	// does *not* change
							// acc and jerk

    grape_was_used_to_calculate_potential = true;	// trigger a reset
							// next time around

    epot =  ekin = etot = 0;

    if (!cm) {
	for_all_leaves(hdyn,b,bb) {
	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	    ekin += 0.5*mi*vel*vel;
	}
    } else {
	for_all_daughters(hdyn,b,bb) {
	    real mi = bb->get_mass();		// replicated code again...
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = bb->get_vel();
	    ekin += 0.5*mi*vel*vel;
	}
    }
    etot = ekin + epot;

#if 0
    cerr << "grape: "; PRC(etot); PRC(ekin); PRL(epot);
    calculate_energies(b, epot, ekin, etot);
    cerr << "host:  "; PRC(etot); PRC(ekin); PRL(epot);
#endif

}


//======================================================================
// grape_calculate_acc_and_jerk: Use the GRAPE hardware to compute the
//				  accs and jerks on a list of nodes.
//======================================================================

local bool get_neighbors_and_adjust_h2(int chip, hdyn * b)

// This function actually sets nn, coll, d_nn_sq and d_coll_sq.

// As in the case of flat_calculate_..., these values overwrite
// previous ones.  Previously set values are ignored.

// Called by:	grape_calculate_acc_and_jerk()

{
    int nnb;
    bool no_nb = false;

    nnb = h3nblist_(&nboards, &chip, h3nb);	// get neighbor list for
						// current neighbor radius

    // We will determine the nn of all stars (and coll, for leaves),
    // and the perturber lists for (the top-level) nodes.

    if (b->is_parent()) {

	if (b->get_oldest_daughter()->get_slow())
	    clear_perturbers_slow_perturbed(b);

	b->new_perturber_list();
    }

    if (nnb < 2) {

	no_nb = true;

    } else {

	// Found at least one neighbor -- find the nearest and
	// determine perturber list for a center-of-mass node.

	real dmin_sq = VERY_LARGE_NUMBER;
	hdyn * bmin = NULL;
	real dcmin_sq = VERY_LARGE_NUMBER;
	hdyn * cmin = NULL;
	int npl = 0;

	hdyn ** pl;
	real rpfac;

	if (b->is_parent()) {
	    pl = b->get_perturber_list();
	    rpfac = b->get_perturbation_radius_factor();
	}

	for (int j = 0; j < nnb; j++) {

	    hdyn *bb = nodes[h3nb[j]];

	    // bb is the j-th neighbor of b.

	    if (bb != b) {

		// Note that it is possible that pred_pos of bb
		// is garbage (never updated)...  Use get_pos
		// instead of get_pred_pos here.

		// Now, get_pred_pos performs Just-In-Time prediction
		// so the above case cannot occur (as of Dec 12, 1996)

		vector diff = b->get_pred_pos() - bb->get_pred_pos();
		real d2 = diff * diff;

		real sum_of_radii = get_sum_of_radii(b, bb);
		update_nn_coll(b, 100,		// (100 = ID)	    // inlined
			       d2, bb, dmin_sq, bmin,
			       sum_of_radii,
			       dcmin_sq, cmin);

		// Recompute the perturber list for parent nodes.
		// See equivalent code for use without GRAPE in
		// hdyn_ev.C/flat_calculate_acc_and_jerk.

		if (b->is_parent()) {

		    if (is_perturber(b, bb->get_mass(),
				     d2, rpfac)) {		    // inlined

			if (npl < MAX_PERTURBERS) {
			    pl[npl] = bb;
			}

			npl++;
		    }
		}
	    }
	}

	if (b->is_parent())
	    b->set_n_perturbers(npl);

	if (bmin == NULL)
	    no_nb = true;
	else {
	    b->set_nn(bmin);
	    b->set_d_nn_sq(dmin_sq);
	}

	if (cmin != NULL) {
	    b->set_coll(cmin);
	    b->set_d_coll_sq(dcmin_sq);
	}	

    }

    if (no_nb) {

    	// No nearest neighbor found --  enlarge the neighbour-sphere
	// radius and try again...

	b->set_grape_rnb_sq(b->get_grape_rnb_sq()*2);

	return false;
    }

    return true;
}

// set_grape4_neighbour_radius: adjust the grape4 neighbour radius
//			       to some reasonable value.

local void set_grape4_neighbour_radius(hdyn * b)

// Called by:	grape_calculate_acc_and_jerk()

{
    int hindex = b->get_grape_index()-1;

    if (b->is_leaf()) {

	// For a single particle, adjust the nnb radius so that it will
	// contain just 1-2 neighbours (set r = sqrt(2)*d_nn, if known).

	if (b->get_nn() != NULL
	    && b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER) {

	    // Node seems to have a valid nearest neighbor pointer.
	    // We can perhaps use the distance as well.
	
	    // Set zero neighbor sphere radius if nb_check_counter is non-zero.

	    if (nb_check_counter[hindex] == 0) {

		// Radius information included here to allow coll
		// criterion to be applied elsewhere.

		b->set_grape_rnb_sq(max(b->get_d_nn_sq(),
					4*b->get_radius()*b->get_radius()));

		// In this case, if the actual value set is zero,
		// one might have to do something.

		if (b->get_grape_rnb_sq() < MINIMUM_GRAPE_RNB_PERT) {

		    // (Jun says this should never happen...)

		    cerr << "h2 set to zero for \n";
		    pp3(b,cerr);
		    PRL(b->get_d_nn_sq());
		    b->set_grape_rnb_sq(b->get_d_min_sq());
		}

	    } else {

		b->set_grape_rnb_sq(0.0);
	    }

	} else {

	    // Node does not know its nearest neighbor.

	    nb_check_counter[hindex] = 0;
	    b->set_grape_rnb_sq(pow(b->get_d_min_sq(), 1.0/3.0));

	    // Note: d_min_sq^(1/3) ~ square of the average interparticle
	    // spacing for standard N-body units -- OK for leaves and nn.
	
	}

    } else {

	// For a node, we will want to compute the perturbers, so
	// we need a larger value of grape_rnb_sq.

	// If the node does not have a valid perturber list,
	// either it has overflowed or nothing has been set before.

	// If it does have non-zero valid perturber list,
	// the perturbation radius should be used.

	real r_bin = 2*binary_scale(b);

	// Use of binary_scale() here is potentially quite inefficient...

	real r_pert2 = max(b->get_perturbation_radius_factor(),
			   r_bin*r_bin);

	// Probably want to clean up the extra distance limits applied
	// to perturbation_radius_factor()...  The extra condition here
	// is to ensure that the harp search radius includes the entire
	// binary.  The choice may only be relevant to GRAPE code, as
	// the non-GRAPE version will check all top-level nodes, and the
	// perturbation_radius condition *ought* to pick up perturbers.

	if (b->get_nn() != NULL
	    && b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER)

	    b->set_grape_rnb_sq(max(r_pert2, b->get_d_nn_sq()));

	else

	    b->set_grape_rnb_sq(r_pert2);

	// Note that this can cause overflow, resulting in
	// not_valid_perturber anyway...

	nb_check_counter[hindex] = 0;
    }

    if (nb_check_counter[hindex] == 0 &&
	b->get_grape_rnb_sq() < MINIMUM_GRAPE_RNB_PERT) {
	cerr << "set_grape4_neighbor ... error \n",
	pp3(b,cerr);
	PRL(nb_check_counter[hindex]);
	PRL(sqrt(b->get_grape_rnb_sq()));
	PRL(b->get_nn());
	PRL(b->get_d_min_sq());
    }
}

//  *****************************
//  *****************************
//  ***                       ***
//  ***  The global function  ***
//  ***                       ***
//  *****************************
//  *****************************

//  If restart is true, we must reinitialize the GRAPE interface
//  after a change in the tree or other kira configuration.
//
//  This function is called from kira_calculate_top_level_acc_and_jerk,
//  which is called only from calculate_acc_and_jerk_for_list.

void grape_calculate_acc_and_jerk(hdyn ** next_nodes,
				  int n_next,
				  xreal time,
				  bool restart)
{
    static int nj_on_grape4;

    hdyn * root = next_nodes[0]->get_root();
    kira_options *ko = root->get_kira_options();

    // if (restart) {
    // 	   cerr << "grape_calculate...: restart\n";
    // }

    // Test the state of the GRAPE and open it if necessary.
    //
    // (The GRAPE release check is now performed externally.  The main
    //  advantage to doing the check here was that we only had to do it
    //  once.  However a major disadvantage was that the hardware could
    //  get tied up unnecessarily by a process that was stuck elsewhere
    //  in the code (e.g. in a multiple encounter.)
    //
    // It is necessary to know when a restart has been triggered ONLY by
    // the release/reattachment of the GRAPE hardware.  Indicator is:
    //
    //		grape_reattached = true

    bool grape_reattached = false;

    if (!grape_is_open) {

	cerr << endl << "grape_calculate_acc_and_jerk: ";
	if (!grape_first_attach) cerr << "re";
	cerr << "attaching GRAPE" << endl;
	h3open_();

	if (!restart) grape_reattached = true;

	set_time_check_mode(0);
	ko->grape_last_cpu = cpu_time();
	int dbgl = 0;
	h3setdebuglevel_(&dbgl);
	grape_is_open = true;

	// If newly opened, restart irrespective of the actual argument.

	restart = true;
    }

    if (grape_was_used_to_calculate_potential) {
	restart = true;
	grape_was_used_to_calculate_potential = false;
    }

    if (restart) {

        // cerr << "grape4 restart\n";

	if (grape_reattached)		// save the nb_check_counter[] array
					// -- unclear why we can't simply modify
					//    the action of initialize_array...
	    for_all_daughters(hdyn, root, bi) {

		putiq(bi->get_dyn_story(), "nb_check_counter",
		      nb_check_counter[bi->get_grape_index()-1]);

	    }

	nj_on_grape4 = initialize_array(root);

	if (grape_reattached)		// Restore the nb_check_counter[] array
	    for_all_daughters(hdyn, root, bi) {

		if (find_qmatch(bi->get_dyn_story(), "nb_check_counter")) {
		    nb_check_counter[bi->get_grape_index()-1] =
			getiq(bi->get_dyn_story(), "nb_check_counter");
		    rmq(bi->get_dyn_story(), "nb_check_counter");
		}
	    }

	n_previous_nodes = 0;
    }

    int i, j;

    // Create the list of TOP_LEVEL nodes in the present block step.

    for (i = j = 0; i < n_next; i++) {
	if (next_nodes[i]->is_top_level_node()) {

	    next_top[j] = next_nodes[i];

	    // Whether prediction should be performed here
	    // or not is pretty unclear....

	    next_top[j]->predict_loworder(time);
	    j++;
	}
    }

    int n_top = j;	// number of top-level nodes in current list

    // Store the predicted positions of top-level
    // current-block nodes to GRAPE memory.

    jpdma_nodes(n_top, next_top, true, time);

    // Store the particles in the previous block which are not in
    // the present block (update GRAPE for previous step).

    for (i = j = 0; i < n_previous_nodes; i++) {
	if (previous_nodes[i]->get_next_time() > time) {
	    previous_nodes[j] = previous_nodes[i];
	    j++ ;
	}
    }

    n_previous_nodes = j;
    jpdma_nodes(n_previous_nodes, previous_nodes, false, 0.0);

    for (i = 0; i < n_top; i++) {

	// store current list

	hdyn * b = previous_nodes[i] = next_top[i];

	// Set some reasonable h2 value.

	set_grape4_neighbour_radius(b);

	if (b->is_parent())
	    b->set_valid_perturbers(true);

    }

    n_previous_nodes = n_top;

    // Now we actually calculate the force...

#ifndef USE_HALF_CHIP
    int nimax = h3npipe_();
#else
    int nimax = h3npipe_()/2;
#endif

    real h2_critical = next_top[0]->get_d_min_sq()*8192;

    // We will stop expanding the GRAPE neighbor sphere once its size
    // exceeds this critical value.  However, it is legal to set
    // grape_rnb_sq greater than h2_critical -- the neighbor sphere
    // then simply won't be expanded if no neighbors are found.

    // *** Should contain a factor of ~(m_max/<m>)^2, but not crucial...

    // Note:  for equal-mass systems and standard units, this critical
    //	      radius is less than the interparticle spacing for N > ~1000.

    real d2_max = h2_critical * 2;

    i = 0;

    while (i < n_top) {

	int inext = i;
	int ip = min(nimax, n_top - i);

	force_by_grape4(time, ip, next_top+i, nj_on_grape4);

	for (int ichip = 0; ichip < ip ; ichip ++) {
	    int real_ichip;

#ifndef USE_HALF_CHIP
	    real_ichip = ichip;
#else
	    real_ichip = ichip*2;
#endif

	    hdyn * b = next_top[i+ichip];
	    int hindex = b->get_grape_index()-1;

	    // PRC(ichip); PRL(real_ichip);

	    if (nb_check_counter[hindex] == 0) {

		if (b->get_grape_rnb_sq() <= 0) {

		    // Should not happen!  About to enter an infinite loop...

		    err_exit("grape_calculate_acc_and_jerk: rnb_sq = 0");

		}

		// Compute neighbors/perturbers iff nb_check_counter = 0.

		if (next_top[i+ichip]->get_grape_rnb_sq() < 1e-20) {
		    cerr << "h2 = 0 for \n";
		    pp3(b,cerr);
		    PRC(hindex);
		    PRL(nb_check_counter[hindex]);
		}

		int inew = i+ichip;
		int h2_too_large = 0;

		while (!(h2_too_large ||
			 get_neighbors_and_adjust_h2(real_ichip, b))) {

		    if (b->get_grape_rnb_sq() > h2_critical) {

			h2_too_large = 1;	// ==> exit while loop

			b->set_nn(b);
			b->set_d_nn_sq(d2_max);

			// cerr << "no nb found for particle ";
			// b->print_label(cerr);
			// cerr << " at position " << b->get_pos();
			// cerr << " with radius " << b->get_grape_rnb_sq()
			//      << endl;

		    } else {

			// Expand the neighbor sphere -- must recompute
			// the force.

			// cerr << "RETRY: no nb found for particle ";
			// b->print_label(cerr);
			// cerr << " at (" << b->get_pos() << ") with radius "
			//      << b->get_grape_rnb_sq()<< endl;

			force_by_grape4(time, 1, next_top+inew, nj_on_grape4);

			real_ichip = 0;
			ichip = ip;

			// Setting ichip = ip here forces exit from the ichip
			// loop when we eventually leave the current while
			// loop.  We will then restart the outermost (i) loop
			// with the *next* particle after this one.
			//
			// Note that, if we somehow get to the "if" part of
			// this structure without first going through this
			// "else" clause, then we will likely make an error.

		    }			// end of if (..) else (...)
		}			// end of while (h2...)
	    }				// end of if (check_...)
	    inext ++;
	}				// end of for (ichip...)
	i = inext;
    }					// end of while (i...)

    for (i = 0; i < n_top ; i++) {

	hdyn * b = next_top[i];
	int hindex = b->get_grape_index()-1;

	// Reduce frequency of nn checks (every fourth force calculation).

	if (b->is_leaf())
	    nb_check_counter[hindex] = (nb_check_counter[hindex] + 1)%4;
	else
	    nb_check_counter[hindex] = 0;

    }
}


//======================================================================
//  grape_calculate_densities:  Determine particle densities, assigning
//				zero density to particles with no
//				neighbor within sqrt(h2_crit).
//======================================================================

#define DEBUG 0

local void set_grape4_density_radius(hdyn * b, real h2_max)
{
    int hindex = b->get_grape_index() - 1;

    // For a single particle, adjust the radius so that it will
    // contain just ~10 neighbours (set r = 4*d_nn, if known).

    if (b->get_nn() != NULL
	&& b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER
	&& b->get_d_nn_sq() > 0) {

	// Node seems to have a valid nearest neighbor pointer.

	b->set_grape_rnb_sq(9 * b->get_d_nn_sq());

    } else

	// Node does not know its nearest neighbor.

	b->set_grape_rnb_sq(9 * pow(b->get_d_min_sq(), 1.0/3.0));

}

local bool count_neighbors_and_adjust_h2(int chip, hdyn * b)
{
    int nnb = h3nblist_(&nboards, &chip, h3nb);		// get neighbor list

    if (nnb < 15) {

	if (DEBUG)
	    cerr << "increasing grape_rnb_sq for " << b->format_label()
		 << " (nnb = " << nnb << ", grape_rnb = "
		 << sqrt(b->get_grape_rnb_sq()) << ")" << endl;

	real fac = 2;
	if (nnb < 4) fac = 4;
	b->set_grape_rnb_sq(fac * b->get_grape_rnb_sq());

	return false;
    }

    // Make a list of nodes to send to compute_com().

    dyn** dynlist = new dynptr[nnb];

    real d_max = 0;
    for (int j = 0; j < nnb; j++) {

	hdyn * bb = nodes[h3nb[j]];
	dynlist[j] = (dyn*)bb;

	// bb is the j-th neighbor of b.

	if (DEBUG) {
	    if (bb != b) {
		vector diff = b->get_pred_pos() - bb->get_pred_pos();
		real d2 = diff * diff;
		d_max = max(d_max, d2);
	    }
	}
    }

    if (DEBUG) {
	real grape_rnb = sqrt(b->get_grape_rnb_sq());
	d_max = sqrt(d_max);
	cerr << b->format_label() << ": ";
	PRC(nnb), PRC(grape_rnb), PRL(d_max);
    }

    compute_density(b, 12, dynlist, nnb);	// writes to dyn story

    // Strangely, compute_density sometimes fails to write the proper
    // info to the dyn story.  Reason and circumstances still unknown.

    if (find_qmatch(b->get_dyn_story(), "density_time")
	&& find_qmatch(b->get_dyn_story(), "density_k_level")
	&& find_qmatch(b->get_dyn_story(), "density")) {

	real density_time = getrq(b->get_dyn_story(), "density_time");
	int  density_k    = getiq(b->get_dyn_story(), "density_k_level");
	real density      = getrq(b->get_dyn_story(), "density");

	if (DEBUG)
	    PRL(density);
    }

    delete [] dynlist;
    return true;
}

//  *****************************
//  *****************************
//  ***                       ***
//  ***  The global function  ***
//  ***                       ***
//  *****************************
//  *****************************

void grape_calculate_densities(hdyn* b,
			       real h2_crit)	// default = 4
{
    if (!grape_is_open) {

	cerr << endl << "grape_calculate_densities: ";
	if (!grape_first_attach) cerr << "re";
	cerr << "attaching GRAPE" << endl;
	h3open_();

	set_time_check_mode(0);
	b->get_kira_options()->grape_last_cpu = cpu_time();
	int dbgl = 0;
	h3setdebuglevel_(&dbgl);
	grape_is_open = true;
    }

    // Copy all nodes to the GRAPE hardware.
    // No need to save counters here, as we will force a restart later...

    int nj_on_grape4 = initialize_array(b);

    // Make a list of all top-level nodes.

    hdyn** top_nodes = new hdynptr[b->n_daughters()];

    int n_top = 0;
    for_all_daughters(hdyn, b, bb)
	top_nodes[n_top++] = bb;

    // Store the predicted positions of top-level
    // current-block nodes to GRAPE memory.

    real time = b->get_system_time();
    jpdma_nodes(n_top, top_nodes, true, time);

    // Set h2 values.

    for (int i = 0; i < n_top; i++)
	set_grape4_density_radius(top_nodes[i], h2_crit);

#ifndef USE_HALF_CHIP
    int nimax = h3npipe_();
#else
    int nimax = h3npipe_()/2;
#endif

    int n_retry = 0;

    int i = 0;
    while (i < n_top) {

	int inext = i;
	int ip = min(nimax, n_top - i);

	force_by_grape4(time, ip, top_nodes+i, nj_on_grape4);

	for (int ichip = 0; ichip < ip ; ichip ++) {

#ifndef USE_HALF_CHIP
	    int real_ichip = ichip;
#else
	    int real_ichip = ichip*2;
#endif

	    hdyn * bb = top_nodes[i+ichip];
	    int hindex = bb->get_grape_index()-1;

	    // Just count neighbors, for now.

	    int inew = i+ichip;
	    int h2_too_large = 0;

	    while (!(h2_too_large ||
		     count_neighbors_and_adjust_h2(real_ichip, bb))) {

		if (bb->get_grape_rnb_sq() > h2_crit) {

		    // Write zero density to dyn story.

		    putrq(bb->get_dyn_story(), "density_time",
			  (real)bb->get_system_time());
	            putrq(bb->get_dyn_story(), "density", 0.0);

		    if (DEBUG) {
			PR(bb->get_grape_rnb_sq());
			cerr << " too large for "
			     << bb->format_label() << endl;
		    }

		    h2_too_large = 1;	// ==> exit while loop

		} else {

		    // Expand the neighbor sphere -- must recompute
		    // the force.

		    force_by_grape4(time, 1, top_nodes+inew, nj_on_grape4);
		    n_retry++;

		    real_ichip = 0;
		    ichip = ip;
		}
	    }

	    if (DEBUG) {
		real dens = getrq(bb->get_dyn_story(), "density");
		cerr << i << "  (" << bb->format_label() << ")  "; PRL(dens);
	    }

	    inext++;
	}
	i = inext;
    }

    // Timestamp the root node.

    putrq(b->get_dyn_story(), "density_time", (real)b->get_system_time());

    if (n_retry > 10) {
	cerr << "grape_calculate_densities:  ";
	PRL(n_retry);
    }

    // Force cleanup later.

    grape_was_used_to_calculate_potential = true;
    delete [] top_nodes;

}


//========================================================================
// External cleanup -- delete local static arrays.
//========================================================================

void clean_up_hdyn_grape()
{
    // Set in put_grape_index_to_top_level_nodes():

    if (nodes) delete [] nodes;
    if (next_top) delete [] next_top;
    if (previous_nodes) delete [] previous_nodes;
    if (nb_check_counter) delete [] nb_check_counter;

    if (pxj) delete [] pxj;
    if (pvj) delete [] pvj;
    if (paj) delete [] paj;
    if (pjj) delete [] pjj;
    if (ptj) delete [] ptj;
    if (pmj) delete [] pmj;
    if (ppj) delete [] ppj;
    if (h3nb) delete [] h3nb;

    // Set in force_by_grape4():

    if (px) delete [] px;
    if (pv) delete [] pv;
    if (pa) delete [] pa;
    if (pj) delete [] pj;
    if (ppot) delete [] ppot;
    if (peps2) delete [] peps2;
    if (ph2) delete [] ph2;

    // Set in force_by_grape4_on_leaves():

    if (pxl) delete [] pxl;
    if (pvl) delete [] pvl;
    if (pal) delete [] pal;
    if (pjl) delete [] pjl;
    if (ppl) delete [] ppl;
    if (peps2l) delete [] peps2l;
    if (ph2l) delete [] ph2l;

    // Set in force_by_grape4_on_all_leaves():

    if (allnodes) delete [] allnodes;
}
