
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//  hdyn_grape6.C:  functions to use GRAPE-6.
//
//    version 1:  Jul 2000   Steve McMillan, Jun Makino
//    version 2:  Dec 2001   Steve McMillan
//.............................................................................
//
//	Externally visible functions (note that the GRAPE version
//	is not specified):
//
//	void check_release_grape
//	void grape_calculate_energies
//	void grape_calculate_acc_and_jerk
//	void grape_calculate_densities
//	void clean_up_hdyn_grape
//
//.............................................................................

#include "hdyn.h"
#include "grape6.h"

#include "kira_debug.h"	// (a handy way to turn on blocks of debugging)
#ifndef T_DEBUG_hdyn_grape6
#   undef T_DEBUG
#endif

#include "hdyn_inline.C"

#define EXPAND_TO_FIND_COLL
#undef EXPAND_TO_FIND_COLL

#define SORT_NODE_LIST
//#undef SORT_NODE_LIST

// Set DEBUG =	0 for no debugging info
//		1 for a substantial amount of high-level info
//		2 for tons of output!

#define DEBUG	0

// Allow possibility of no neighbor list (former SPZ convention
// reversed by Steve, 6/01):

//#define NO_G6_NEIGHBOUR_LIST

// Convenient to allow inline functions to be separated out for
// debugging and profiling purposes.

#define INLINE	inline
//#define INLINE

#define ONE2	0.5
#define ONE6	0.16666666666666666667

#define NRETRY	10

// Certain GRAPE functions will call themselves up to NRETRY times in
// the event of a hardware error, to try to fix the problem.  A nonzero
// return thus means that NRETRY hardware resets have already been
// attempted, and have failed.  The calling function should probably
// then stop and call Jun...

// Static declarations:
// -------------------

static int cluster_id = 0;		// GRAPE-6 cluster to use
static real grape_time_offset = 0;	// offset time to maintain precision

// Static j-arrays, created in initialize_grape_arrays:

static hdyn **node_list = NULL;		// "reverse index pointers":
					// node_list[i]->get_grape_index() = i

static hdyn **current_nodes = NULL;	// current top-level nodes
static hdyn **previous_nodes = NULL;	// previous top-level nodes

// Arrays current_nodes and previous_nodes are used only in
// grape_calculate_acc_and_jerk().

static int grape_nj_max = 0;		// maximum size of static j-arrays
static int n_previous_nodes = 0;	// number of "previous" nodes

// GRAPE state indicators:

static bool grape_is_open = false;
static bool grape_first_attach = true;
static bool grape_was_used_to_calculate_potential = false;



//  **************************************************************************
//  *                                                                        *
//  * Local functions used by more than one other function (global or local) *
//  *                                                                        *
//  **************************************************************************
//
//		reattach_grape
//		send_j_node_to_grape
//		initialize_grape_arrays
//		force_by_grape
//		hw_err_exit



local void reattach_grape(real time, char *id, kira_options *ko)

// Reattach the GRAPE-6.

// Called by:	grape_calculate_energies()			// global
//		grape_calculate_acc_and_jerk()			// global
//		grape_calculate_densities()			// global

{
	cerr << id << ":  ";
	if (!grape_first_attach) cerr << "re";
	cerr << "attaching GRAPE at time "
	     << time << endl << flush;

	g6_open_(&cluster_id);
	if (ko) ko->grape_last_cpu = cpu_time();
	grape_is_open = true;

	if (DEBUG) {
	    cerr << "GRAPE successfully attached"
		 << endl << flush;
	}
	PRL(cpu_time());
}



local INLINE void send_j_node_to_grape(hdyn *b,
				       bool computing_energy = false)

// Send node b to the GRAPE-6 j-particle list.

// Called by:	initialize_grape_arrays				// local
//		send_all_leaves					// local
//		grape_calculate_acc_and_jerk			// global

// The flag computing_energy will be used only in the energy calculation,
// to ensure that no hardware prediction is done on the GRAPE during the
// force computation, and that all positions are expressed relative to the
// root node.

{
    if (DEBUG > 1) {
	cerr << "  send_j_node_to_grape:  ";
	PRL(b->format_label());
    }

    int grape_index = b->get_grape_index();

    // If we are computing the energy, we may be doing it at an odd time,
    // which may cause the GRAPE to complain.  Since we never do prediction
    // in this case, just send a time of zero to avoid prediction and warning
    // messages. Note that this convention *must* be consistent with the actions
    // taken by force_by_grape.

    real t = 0;
    if (!computing_energy)
	t = b->get_time() - grape_time_offset;			// 0 <= t <= 1

    real dt = b->get_timestep();
    if (dt <= 0) dt = 1;
    real m = b->get_mass();

    vector pos;
    if (computing_energy)
	pos = hdyn_something_relative_to_root(b, &hdyn::get_pos);
    else
	pos = b->get_pos();

    vector vel = b->get_vel();

    vector a2 = ONE2 * b->get_acc();		// would be more efficient
    vector j6 = ONE6 * b->get_jerk();		// to store these in hdyn...

    vector k18 = b->get_k_over_18();		// is stored in hdyn

    if (DEBUG > 1) {
	PRI(4); PRC(cluster_id); PRL(grape_index);
	PRI(4); PRC(t); PRC(dt); PRL(m);
	PRI(4); PRL(j6);
	PRI(4); PRL(a2);
	PRI(4); PRL(vel);
	PRI(4); PRL(pos);
    }
    
    g6_set_j_particle_(&cluster_id,
		       &grape_index,				// address
		       &grape_index,				// index
		       &t, &dt, &m,
		       &k18[0],					// for now
		       &j6[0],
		       &a2[0],
		       &vel[0],
		       &pos[0]);

    if (DEBUG > 1) {
	cerr << "  ...sent to GRAPE-6" << endl << flush;
    }
}



local int initialize_grape_arrays(hdyn *b,		// root node
				  bool reset_counters = true)

// Initialize arrays associated with GRAPE-6 j-particles.
// Return the total number of j-particles sent to the GRAPE.

// Called by:	grape_calculate_acc_and_jerk()			// global
//		grape_calculate_densities()			// global

{
    if (DEBUG) {
	cerr << "initialize_grape_arrays at time "
	     << b->get_real_system_time()<< ":  ";
	PRL(reset_counters);
    }

    int nj = b->n_daughters();

    if (DEBUG)
	PRL(nj);

    if (grape_nj_max <= 0 || grape_nj_max < nj) {

	// Reset the node_list array.

	if (node_list) delete [] node_list;
	if (current_nodes) delete [] current_nodes;
	if (previous_nodes) delete [] previous_nodes;

	grape_nj_max = 3*nj + 10;				// cautious!
	node_list = NULL;					// flag to set
    }

    if (!node_list) {

	if (DEBUG)
	    PRL(grape_nj_max);

	node_list = new hdynptr[grape_nj_max];
	current_nodes = new hdynptr[grape_nj_max];
	previous_nodes = new hdynptr[grape_nj_max];
    }

    // See if we need to increase the time offset.

    while (b->get_real_system_time() - grape_time_offset > 1)
	grape_time_offset += 1;

    // Set up the j-particle data.

    nj = 0;
    for_all_daughters(hdyn, b, bb) {

	// For now, let GRAPE index = address.

	bb->set_grape_index(nj);
	node_list[nj++] = bb;
	
	// Copy the node to the GRAPE with index = address.

	send_j_node_to_grape(bb);

	if (reset_counters)
	    bb->set_grape_nb_count(0);
    }

    if (DEBUG) {
	cerr << endl << "Initialized GRAPE-6 arrays, ";	PRL(nj);
	cerr << "...leaving initialize_grape_arrays" << endl;
    }

    return nj;
}



local hdyn *find_and_print_nn(hdyn *b)
{
    real d2min = VERY_LARGE_NUMBER;
    hdyn *bmin = NULL;
    for_all_daughters(hdyn, b->get_root(), bb) {
	if (bb != b) {
	    real d2 = square(bb->get_pred_pos() - b->get_pred_pos());
	    if (d2 < d2min) {
		d2min = d2;
		bmin = bb;
	    }
	}
    }
    if (bmin) {
	cerr << "    found true nn of " << b->format_label();
	cerr << " = " << bmin->format_label()
	     << ";  grape_index = " << bmin->get_grape_index()
	     << "  at distance " << sqrt(d2min) << endl;
    } else
	cerr << "    can't find true nn for " << b->format_label() << endl;

    return bmin;
}



local void reset_grape(hdyn *b)					// root node

// Perform a "hard" reset of the grape interface to try to recover
// from a hardware error.

// Called by:	force_by_grape					// global

{
    // First two lines will result in a slow reopening of the GRAPE
    // (as when the hardware is attached for the first time).

    g6_reset_(&cluster_id);
    g6_reset_fofpga_(&cluster_id);

    g6_close_(&cluster_id);
    g6_open_(&cluster_id);

    cerr << endl << "*** Reset GRAPE-6 at time "
	 << b->get_real_system_time() << endl;

    // Reload the j-arrays.

    initialize_grape_arrays(b, false);
}



// These arrays are actually local to force_by_grape, but it is convenient
// to make them globally accessible in order to permit cleanup.

static int    *iindex = NULL;
static vector *ipos   = NULL;
static vector *ivel   = NULL;
static vector *iacc   = NULL;
static vector *ijerk  = NULL;
static real   *ipot   = NULL;
static real   *ih2    = NULL;
static int    *inn    = NULL;
static real   eps2    = 0;

//-------------------------------------------------------------------------

local INLINE int force_by_grape(xreal xtime,
				hdyn *nodes[], int ni,
				int nj, int n_pipes,
				bool pot_only = false,
				int level = 0)

// Calculate the force on ni nodes due to nj particles already
// stored in GRAPE memory, using up to n_pipes hardware pipes.
// Return 0 iff no problems occurred.

// Called by:	force_by_grape_on_all_leaves			// local
//		get_force_and_neighbors				// local
//		grape_calculate_densities()			// global

// Flag pot_only = true means that we are calculating the total energy.
// As of 11/02, we ignore the value of xtime in this case, although it may
// still be used for debugging purposes (Steve).

{
    static char *func = "force_by_grape";

    if (ni <= 0) return 0;
    if (ni > n_pipes) err_exit("force_by_grape: ni too large\n");

    if (DEBUG > 1) {
	cerr << "  " << func << ":  ";
	PRC(xtime); PRC(ni); PRC(nj); PRL(pot_only);
    }

    if (!iindex) {

	// Create the i-particle arrays.

	iindex = new int[n_pipes];
	ipos   = new vector[n_pipes];
	ivel   = new vector[n_pipes];
	iacc   = new vector[n_pipes];
	ijerk  = new vector[n_pipes];
	ipot   = new real[n_pipes];
	ih2    = new real[n_pipes];
	inn    = new int[n_pipes];
	eps2   = nodes[0]->get_eps2();
    }

    // Send the current time to the GRAPE.

    real time = 0;
    if (!pot_only) time = xtime - grape_time_offset;
    g6_set_ti_(&cluster_id, &time);

    // Pack the i-particle data and start the GRAPE calculation.

    if (DEBUG > 1) {
	cerr << "  loading nodes..."
	     << endl << flush;
    }

//    PRL(ni);

#ifdef T_DEBUG
    real sys_t = xtime;
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 1 << endl << flush;
    }
#endif

    for (int i = 0; i < ni; i++) {

//	PRC(i); PRL(nodes[i]);

	iindex[i] = nodes[i]->get_grape_index();

	// New version suggested by Jun (May 8 2001) to prevent the
	// "call Jun" message.  Implemented by SPZ.  Corrected by Steve.
	// We must avoid sending zeroes to g6calc_firsthalf() in acc,
	// jerk, or pot, as these are used for scaling purposes.

	if (pot_only) {

            // Nodes in this case are leaves, not necessarily top-level.
	    // Don't extrapolate in this case (see 11/02 note below).

            ipos[i] = hdyn_something_relative_to_root(nodes[i],
//						      &hdyn::get_pred_pos);
						      &hdyn::get_pos);

        } else
            ipos[i]   = nodes[i]->get_pred_pos();

        ivel[i]   = nodes[i]->get_pred_vel();
	iacc[i]   = nodes[i]->get_old_acc();
	ijerk[i]  = ONE6*nodes[i]->get_old_jerk();
        ipot[i]   = nodes[i]->get_pot();
        ih2[i]    = nodes[i]->get_grape_rnb_sq();

	// Must take care of pot, acc, and jerk.  Quite expensive to
	// check the abs of a vector, so try more indirect tests first.
	// Normally,  these values will be zero only at initialization,
	// but we may also encounter this immediately after a new CM
	// node is created.  Also, if we are calculating the energy,
	// there may be problems with low-level nodes.

	if (nodes[i]->is_top_level_node()) {

	    if (abs(ipot[i]) <= VERY_SMALL_NUMBER
		|| abs(abs(iacc[i][0])) <= VERY_SMALL_NUMBER) {

		// Assume that pot, acc, jerk are not set properly.
		// For now, assume that the average numbers appropriate
		// to the mean field of a standard cluster are OK.

		ipot[i] = 1;

		if (abs(iacc[i]) <= VERY_SMALL_NUMBER ||
		    abs(ijerk[i]) <= VERY_SMALL_NUMBER) {

		    // cerr << "WARNING: Initializing acc and jerk from zero."
		    //      << endl;

		    // Small values here will cause overflow and annoying
		    // messages from the GRAPE.  Large values will cause
		    // underflow and no messages.

		    iacc[i]   = vector(1);
		    ijerk[i]  = vector(1);

		    // Hmmm.  Increasing these numbers seems to increase
		    // the number of error messages...
		}
	    }

	} else {

	    // Only interested in a low-level node during the potential
	    // calculation.  In this case, it is quite unlikely that
	    // any of pot, acc, or jerk is usable for scaling purposes.
	    // Choose quantities appropriate to the nearest neighbor.

	    hdyn *sis = nodes[i]->get_younger_sister();
	    if (!sis) sis = nodes[i]->get_elder_sister();
	    real m = sis->get_mass();
	    real r, r2;

	    kepler *k = nodes[i]->get_kepler();
	    if (k) {
		r = k->get_separation();
		r2 = r*r;
	    } else {
		r2 = square(nodes[i]->get_pos() - sis->get_pos());
		r = sqrt(r);
	    }

	    real j = m/(r*r2);
	    ipot[i] = j*r2;
	    iacc[i] = j*r;
	    ijerk[i] = j;
	}

	if (DEBUG > 1) {
	    if (i > 0) cerr << endl;
	    PRI(4); PRC(i); PRC(iindex[i]); PRL(nodes[i]->format_label());
	    PRI(4); PRL(ipos[i]);
	    PRI(4); PRL(ivel[i]);
	    PRI(4); PRL(iacc[i]);
	    PRI(4); PRL(ijerk[i]);
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 2 << endl << flush;
    }
#endif

    // Clear neighbor radii for unused pipes (added by SPZ, May 8 2001).

    for (int i = ni; i < n_pipes; i++)
        ih2[i] = 0;

    g6calc_firsthalf_(&cluster_id, &nj, &ni, iindex,
		      ipos, ivel, iacc, ijerk, ipot,
		      &eps2, ih2);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 3 << endl << flush;
    }
#endif

    int error = g6calc_lasthalf2_(&cluster_id, &nj, &ni, iindex,
				  ipos, ivel, &eps2, ih2,
				  iacc, ijerk, ipot, inn);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 4 << endl << flush;
	PRL(error);
    }
#endif

    if (DEBUG > 1) {
	cerr << endl << "  ...results:"
	     << endl << flush;
    }

    if (!error) {

	for (int i = 0; i < ni; i++) {

	    // Hmmm...  Looks like it is possible for data to return
	    // wrong even though the error flag is not set.  The best
	    // we can do for now is to check that inn is within range
	    // and pot < 0.  Not perfect, but... (Steve, 7/00)

             if (inn[i] < 0 ||  ipot[i] > 0) {
                 error = 42;
                 break;
             }
             if (inn[i] >= nj) {
                 inn[i] = 0;
                 cerr << "warning: NN forced to zero for particle "<<endl;  
                 PRI(4); PRC(i); PRC(iindex[i]);
                 PRL(nodes[i]->format_label());
                 PRL(nodes[i]->get_pos());
	     }

	    if (pot_only)

		nodes[i]->set_pot(ipot[i]);

	    else {

		hdyn *nn = node_list[inn[i]];
		real d_nn_sq = 0;

		if (DEBUG > 1) {
		    if (i > 0) cerr << endl;
		    PRI(4); PRC(i); PRC(iindex[i]);
		    	    PRL(nodes[i]->format_label());
		    PRI(4); PRL(iacc[i]);
		    PRI(4); PRL(ijerk[i]);
		    PRI(4); PRL(ipot[i]);
		    PRI(4); PRC(inn[i]); PR(nn);
		}

		if (nn) {

		    d_nn_sq = square(ipos[i] - nn->get_pred_pos());

		    if (DEBUG > 1) {
			PRI(2); PRL(nn->format_label());
		    }

		} else {

		    // Set nn = nodes[i] for handling elsewhere.

		    nn = nodes[i];

		    if (DEBUG > 1)
			cerr << endl;
		}

		if (DEBUG > 1) {
		    hdyn *true_nn = find_and_print_nn(nodes[i]);
		    if (true_nn != nn)
			cerr << "    *** error: nn != true_nn" << endl;
		}
		
		nodes[i]->set_acc_and_jerk_and_pot(iacc[i], ijerk[i], ipot[i]);
		nodes[i]->set_nn(nn);
		nodes[i]->set_d_nn_sq(d_nn_sq);


#if 0000
		if (nodes[i]->name_is("(5394,21337)")
		     || nodes[i]->name_is("(21337,5394)")) {
		    cerr << func << ": ";
		    PRC(nodes[i]); PRC(inn[i]); PRC(nn); PRL(d_nn_sq);
		}
#endif


		if (DEBUG > 1) {

		    // Cross-check:

		    bool found_index = false;
		    for_all_daughters(hdyn, nodes[0]->get_root(), bb)
			if (bb->get_grape_index() == inn[i]) {
			    cerr << "    check: found grape_index = " << inn[i]
				 << " for node " << bb->format_label()
				 << endl;
			    found_index = true;
			    if (bb != nn)
				cerr << "    *** error: nn mismatch" << endl;
			}
		    if (!found_index)
			cerr << "    *** error: grape_index " << inn[i]
			     << " not found" << endl;
		}
	    }
	}
    }

    // Retest error because it could have been set in the i-loop above.

    if (error) {

	// Make NRETRY attempts (recursively) to reset and correct the
	// error before returning a "true" error condition.

	cerr << "*** " << func << "(" << level << "):  hardware error "
	     << error << " at time " << time << ",  ni = " << ni << endl;

	if (level < NRETRY) {

	    cerr << "Resetting GRAPE and retrying..." << endl;
	    reset_grape(nodes[0]->get_root());

	    error = force_by_grape(xtime, nodes, ni, nj, n_pipes, pot_only,
				   level+1);
	}
    }

    if (DEBUG > 1) {
	cerr << "  ...leaving " << func << "(" << level <<"):  ";
	PRL(error);
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 5 << endl << flush;
	PRL(error);
    }
#endif

    return error;
}



local void hw_err_exit(char *func, int id, hdyn *b)

// Exit following a serious hardware error...

// Called by:	grape_calculate_energies()			// global
//		grape_calculate_acc_and_jerk()			// global
//		grape_calculate_densities()			// global

{
    char buf[256];
    sprintf(buf, "%s[%d]:  time to call Jun!", func, id);

    // Dump out a copy of the system.

    char *dumpfile = "hw_error_dump";

    ofstream dump(dumpfile);
    if (dump) {
	put_node(dump, *b->get_root(), b->get_kira_options()->print_xreal);
	dump.close();
	cerr << "Data written to file " << dumpfile << endl;
    }
    
    err_exit(buf);
}



//  *************************************************************************
//  *************************************************************************
//  **									   **
//  **  Globally visible GRAPE-6 functions (and dedicated local helpers).  **
//  **									   **
//  *************************************************************************
//  *************************************************************************


//  **********************************************************************
//  *                                                                    *
//  *  check_release_grape:  Accessor for GRAPE release/attach.          *
//  *                                                                    *
//  **********************************************************************


void check_release_grape(kira_options *ko, xreal time, bool verbose)
{
#ifdef SHARE_GRAPE

    if (DEBUG) {
	cerr << "GRAPE CPU check:  ";
	PRL(cpu_time());
    }

    if (cpu_time() - ko->grape_last_cpu > ko->grape_max_cpu) {

	int p = cerr.precision(STD_PRECISION);
	cerr << endl << "Releasing GRAPE-6 at time " << time << " after ";
	cerr.precision(2);
	real current_cpu = cpu_time();
	cerr << current_cpu - ko->grape_last_cpu <<" CPU sec" << endl;
	cerr.precision(p);
	PRL(current_cpu);

	g6_close_(&cluster_id);
	grape_is_open = false;
	grape_first_attach = false;

	cerr << endl;
    }

#endif
}



//  *********************************************************************
//  *                                                                   *
//  *  grape_calculate_energies:  Calculate total energy of the system  *
//  *			          (requires GRAPE reset after use).     *
//  *                                                                   *
//  *********************************************************************

local inline bool use_cm_approx(hdyn *bb, real d_crit)
{
    // Check if a binary can be treated as unperturbed for purposes of
    // computing its energy if the estimated tidal effect of its neighbors
    // is negligible.  Looks as though accepting any unperturbed binary
    // may lead to unacceptably large tidal errors.  Use d_crit to define
    // a "close" binary to be treated as unperturbed -- to be reviewed.
    //
    // 							(Steve, 5/02)

    bool use_cm = false;

    if (bb->is_low_level_node()) {

	kepler *k = bb->get_kepler();

	if (k) {

	    if (k->get_semi_major_axis() < d_crit) {

		use_cm = true;

#if 0
		cerr << "use_cm = true (kepler) for " << bb->format_label()
		     << endl;
#endif

	    }

	} else {

	    hdyn *sis = bb->get_younger_sister();
	    if (sis) {

		vector dx = bb->get_pos() - sis->get_pos();

		if (abs(dx[0]) < d_crit
		    && abs(dx[1]) < d_crit
		    && abs(dx[2]) < d_crit) {

		    use_cm = true;

#if 0
		    cerr << "use_cm = true (close) for " << bb->format_label()
			 << endl;
		    PRL(dx);
#endif

		}

#if 0
		if (!use_cm) {
		    cerr << "no CM for " << bb->format_label() << endl;
		    PRC(d_crit); PRL(abs(dx));
		} else {
		    cerr << "using CM for " << bb->format_label() << endl;
		    PRC(d_crit); PRL(abs(dx));
		}
#endif

	    }
	}

#if 0
	if (use_cm) {
	    cerr << "CM approx: ";
	    hdyn *par = bb->get_parent();
	    PRC(par->format_label()); PRL(bb->format_label());
	}
#endif

    }

    return use_cm;
}




//*************************************************************************
//*************************************************************************
//
// **** Note from Steve, 11/02: ****
//
// From now on, we *never* extrapolate positions or velocities when
// computing energies.  This is consistent with the use of the non-GRAPE
// version of this function (calculate_energies) and, in normal use, the
// system should always be synchronized when the energy is computed, so
// extrapolation is unnecessary.  Only unperturbed binaries may be out of
// sync.  However, (1) they cannot be predicted in the usual way and (2)
// their internal energies should be correct anyway, so not extrapolating
// is OK here too.  In addition, other programs (e.g. sys_stats) using this
// function may not know the acc and jerk, so again prediciton is dangerous.
//
//*************************************************************************
//*************************************************************************
//
// Experimental code for grape_calculate_energies().
//
// Basic compile-time options:
//
//	OLD_MAY02:	the original version, prior to May 26, 2002
//	NEW_MAY02:	the new version, with various "enhancements"
//
// In the latter case, choose between an intermediate version, cosmetically
// similar to the new code but functionally more like the original, and the
// true new version, by use of E_NODES (set ==> new).
//
// All versions share the function use_cm_approx() and produce the same
// results.  The NEW/E_NODES version is cleanest, but the OLD version still
// seems marginally fastest...  
//
//							(Steve, 5/02)

#define OLD_MAY02
#define NEW_MAY02

// NOTE: NEW takes precedence over OLD if both are set...
//	 E_NODES is only relevant in the NEW case.

#define E_NODES
//#undef E_NODES

//=========================================================================
//
// New code.

#ifdef NEW_MAY02

#ifdef E_NODES
local int send_all_e_nodes_to_grape(hdyn *b,		// root node
				    bool cm,		// use CM approx
				    real d_crit,	// CM cutoff
				    hdyn **e_nodes,	// list of nodes
				    real& e_unpert,	// unperturbed energy
				    int& n_leaves)	// leaves on list
#else
local int send_all_leaves_to_grape(hdyn *b,		// root node
				   bool cm,		// use CM approx
				   real d_crit, 	// CM cutoff
				   real& e_unpert)	// unperturbed energy
#endif

// Send predicted positions (velocities, etc. are not needed to
// compute the potentials) of all leaves to GRAPE j-particle memory.
// Return the total number of leaves sent to the GRAPE.

// Called by:	grape_calculate_energies()			// global

{
    if (DEBUG) {
	cerr << "  send_all_e_nodes_to_grape"
	     << endl << flush;
    }

    int nj = 0;
    e_unpert = 0;
#ifdef E_NODES
    n_leaves = 0;
#endif

    if (!cm) {

	for_all_leaves(hdyn, b, bb) {

	    // In center of mass approximation, replace bb by its parent.
	    // The order of traversal of the tree means that this should
	    // occur at the elder component of a binary, and will skip the
	    // other component.  The while loop should also take care of
	    // unperturbed multiples!  For efficiency, keep and return a
	    // running sum of all internal energies excluded.

	    bool reset_bb = false;

	    while (use_cm_approx(bb, d_crit)) {
		hdyn *sis = bb->get_younger_sister();
		hdyn *par = bb->get_parent();
		real reduced_mass = bb->get_mass() * sis->get_mass()
						   / par->get_mass();
		if (bb->get_kepler())
		    e_unpert += reduced_mass * bb->get_kepler()->get_energy();
		else {
		    vector dx = bb->get_pos() - sis->get_pos();
		    vector dv = bb->get_vel() - sis->get_vel();
		    e_unpert += reduced_mass * (0.5*square(dv)
						- par->get_mass()/abs(dx));
		}
		bb = par;
		reset_bb = true;
	    }

	    // Must be careful to avoid an infinite loop (because of the
	    // logic used by for_all_leaves():

	    if (reset_bb) {
		bb = bb->get_oldest_daughter()->get_younger_sister()
					      ->next_node(b);
		if (!bb) break;
	    }

	    // For now, let GRAPE index = address on list.

#ifdef E_NODES
	    e_nodes[nj] = bb;
	    n_leaves += bb->is_leaf();
#endif
	    bb->set_grape_index(nj++);
	
	    // Copy the leaf to the GRAPE with index = address.
	    // Send system_time to avoid prediction.

	    send_j_node_to_grape(bb, true);
	}

    } else {

	for_all_daughters(hdyn, b, bb) {

	    // For now, let GRAPE index = address on list.

#ifdef E_NODES
	    e_nodes[nj] = bb;
	    n_leaves += bb->is_leaf();
#endif
	    bb->set_grape_index(nj++);
	
	    // Copy the node to the GRAPE with index = address.
	    // Send system_time to avoid prediction.

	    send_j_node_to_grape(bb, true);
	}
    }

    if (DEBUG) {
	cerr << "  ...leaving send_all_e_nodes_to_grape:  ";
	PRL(nj);
    }

    return nj;
}


#ifdef E_NODES

local bool force_by_grape_on_all_e_nodes(hdyn **e_nodes,    // node list
					 int nj)	    // number on list

// Compute the forces on all e_nodes due to all other e_nodes.
// Only interested in determining the potential energies.

// Called by:	grape_calculate_energies()			// global

{
    if (DEBUG) {
	cerr << "  force_by_grape_on_all_e_nodes:  ";
	PRL(nj);
    }

    bool status = false;

    if (nj > 0) {

	xreal sys_t = e_nodes[0]->get_system_time();
	int n_pipes = g6_npipes_();

	// Compute particle forces, n_pipes at a time.

	for (int ip = 0; ip < nj; ip += n_pipes)
	    status |= force_by_grape(sys_t,
				     e_nodes+ip,Starlab:: min(n_pipes, nj-ip),
				     nj, n_pipes,
				     true);	// "true" ==> compute pot only

    }

    if (DEBUG) {
	cerr << "  ...leaving force_by_grape_on_all_e_nodes:  ";
	PRL(status);
    }

    return status;
}

#else


local bool force_by_grape_on_all_leaves(hdyn *b,		// root node
					int nj,
					bool cm)		// CM approx

// Compute the forces on all particles due to all other particles.
// Only interested in determining the potential energies.

// Called by:	grape_calculate_energies()			// global

{
    if (DEBUG) {
	cerr << "  force_by_grape_on_all_leaves:  ";
	PRL(nj);
    }

    int n_pipes = g6_npipes_();
    hdyn **ilist = new hdynptr[n_pipes];

    // Compute particle forces, n_pipes at a time.

    bool status = false;

    int ip = 0;
    if (!cm) {
	for_all_leaves(hdyn, b, bb) {
	    ilist[ip++] = bb;
	    if (ip == n_pipes) {
		status |= force_by_grape(b->get_system_time(),
					 ilist, ip, nj, n_pipes,
					 true);	// "true" ==> compute pot only
		ip = 0;
	    }
	}
    } else {
	for_all_daughters(hdyn, b, bb) {
	    ilist[ip++] = bb;			// replicated code...
	    if (ip == n_pipes) {
		status |= force_by_grape(b->get_system_time(),
					 ilist, ip, nj, n_pipes,
					 true);	// "true" ==> compute pot only
		ip = 0;
	    }
	}
    }

    if (ip)
	status |= force_by_grape(b->get_system_time(),
				 ilist, ip, nj, n_pipes,
				 true);

    if (DEBUG) {
	cerr << "  ...leaving force_by_grape_on_all_leaves:  ";
	PRL(status);
    }

    delete [] ilist;
    return status;
}

#endif


//						*****************************
//						*****************************
//						***                       ***
//						***  The global function  ***
//						***                       ***
//						*****************************
//						*****************************


local inline real critical_sep(hdyn *b)
{
    // Standard choice of radius to treat a binary as unperturbed,
    // for purposes of computing the energy.  Parameters are fairly
    // arbitrary, but GRAPE precision is such that d_crit can be chosen
    // conservatively small to minimize the tidal effects.

    real d_crit = 1.e-9;
    if (b->get_d_min_sq() > 0)
	d_crit = Starlab::min(0.01*sqrt(b->get_d_min_sq()), d_crit);
    return d_crit;
}


void grape_calculate_energies(hdyn *b,			// root node
			      real &epot,
			      real &ekin,
			      real &etot,
			      bool cm)			// default = false
{
    // Note: cm = true means use the CM approximation for all top-level
    // nodes.  If cm = false, then still use the CM approximation for
    // sufficiently close binaries to avoid roundoff errors on the GRAPE.

//    if (DEBUG) {
	cerr << endl << "entering grape_calculate_energies... ";
	PRL(cm);
//    }

    if (!grape_is_open)
	reattach_grape(b->get_real_system_time(),
		       "grape_calculate_energies", b->get_kira_options());

    real cpu0 = cpu_time();
    real d_crit = critical_sep(b);
    real e_unpert;

#ifdef E_NODES

    // New code (Steve, 5/02): make a list of nodes to use, taking
    // into account the possible use of the CM approximation.

    int n_top_level = b->n_daughters();
    hdyn **e_nodes = new hdynptr[3*n_top_level+10];	// (? same as above)

    int n_leaves;
    int nj =  send_all_e_nodes_to_grape(b, cm, d_crit,
					e_nodes, e_unpert, n_leaves);


    if (force_by_grape_on_all_e_nodes(e_nodes, nj)) {

#else

    int nj =  send_all_leaves_to_grape(b, cm, d_crit, e_unpert);

    if (force_by_grape_on_all_leaves(b, nj, cm)) {

#endif

	cerr << "grape_calculate_energies: "
	     << "error on return from force_by_grape_on_all_xxx()"
	     << endl;

	hw_err_exit("grape_calculate_energies", 1, b);
    }

    epot = ekin = 0;

#ifdef E_NODES

	for (int i = 0; i < nj; i++) {
	    hdyn *bb = e_nodes[i];
	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	    ekin += 0.5*mi*vel*vel;

#if 0
	    if (bb->is_low_level_node() && bb->get_kepler()) {
		PRC(bb->format_label());
		PRC(bb->get_pot());
		real pe = -bb->get_binary_syster()->get_mass()
		    / bb->get_kepler()->get_separation();
		PRL(pe);
	    }
#endif

	}
#else

    if (!cm) {

	for_all_leaves(hdyn, b, bb) {

	    // Logic here follows that in send_all_e_nodes_to_grape().

	    bool reset_bb = false;

	    while (use_cm_approx(bb, d_crit)) {
		bb = bb->get_parent();
		reset_bb = true;
	    }
	    if (reset_bb) {
		bb = bb->get_oldest_daughter()->get_younger_sister()
					      ->next_node(b);
		if (!bb) break;
	    }

	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	    ekin += 0.5*mi*vel*vel;
	}
    } else {
	for_all_daughters(hdyn, b, bb) {
	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = bb->get_vel();
	    ekin += 0.5*mi*vel*vel;
	}
    }

#endif

    etot = ekin + epot + e_unpert;

    grape_was_used_to_calculate_potential = true;	// trigger a reset
							// next time around

//    cerr << "CPU time for grape_calculate_energies() = "
//	 << cpu_time() - cpu0 << endl;

//    if (DEBUG) {
	cerr << "...leaving grape_calculate_energies...  ";
	PRL(etot);
	cerr << endl;
//    }

#ifdef E_NODES
    delete [] e_nodes;
#endif
}


//=========================================================================

#else

// Old code.

local int send_all_leaves_to_grape(hdyn *b,		// root node
				   real& e_unpert,	// unperturbed energy
				   bool cm = false)	// use CM approx

// Send predicted positions (velocities, etc. are not needed to
// compute the potentials) of all leaves to GRAPE j-particle memory.
// Return the total number of leaves sent to the GRAPE.

// Called by:	grape_calculate_energies()			// global

{
    if (DEBUG) {
	cerr << "  send_all_leaves_to_grape"
	     << endl << flush;
    }

    int nj = 0;
    e_unpert = 0;

    if (!cm) {

	real d_crit = critical_sep(b);

	for_all_leaves(hdyn, b, bb) {

	    // In center of mass approximation, replace bb by its parent.
	    // The order of traversal of the tree means that this should
	    // occur at the elder component of a binary, and will skip the
	    // other component.  The while loop should also take care of
	    // unperturbed multiples!  For efficiency, keep and return a
	    // running sum of all internal energies excluded.

	    bool reset_bb = false;

	    while (use_cm_approx(bb, d_crit)) {
		hdyn *sis = bb->get_younger_sister();
		hdyn *par = bb->get_parent();
		real reduced_mass = bb->get_mass() * sis->get_mass()
						   / par->get_mass();
		if (bb->get_kepler())
		    e_unpert += reduced_mass * bb->get_kepler()->get_energy();
		else {
		    vector dx = bb->get_pos() - sis->get_pos();
		    vector dv = bb->get_vel() - sis->get_vel();
		    e_unpert += reduced_mass * (0.5*square(dv)
						- par->get_mass()/abs(dx));
		}
		bb = par;
		reset_bb = true;
	    }

	    // Must be careful to avoid an infinite loop (because of the
	    // logic used by for_all_leaves():

	    if (reset_bb) {
		bb = bb->get_oldest_daughter()->get_younger_sister()
					      ->next_node(b);
		if (!bb) break;
	    }

	    // For now, let GRAPE index = address.

	    bb->set_grape_index(nj++);
	
	    // Copy the leaf to the GRAPE with index = address.
	    // Send system_time to avoid prediction.

	    send_j_node_to_grape(bb, true);
	}

    } else {

	for_all_daughters(hdyn, b, bb) {

	    // For now, let GRAPE index = address.

	    bb->set_grape_index(nj++);
	
	    // Copy the node to the GRAPE with index = address.
	    // Send system_time to avoid prediction.

	    send_j_node_to_grape(bb, true);
	}
    }

    if (DEBUG) {
	cerr << "  ...leaving send_all_leaves_to_grape:  ";
	PRL(nj);
    }

    return nj;
}



local bool force_by_grape_on_all_leaves(hdyn *b,		// root node
					int nj,
					bool cm = false)	// CM approx

// Compute the forces on all particles due to all other particles.
// Only interested in determining the potential energies.

// Called by:	grape_calculate_energies()			// global

{
    if (DEBUG) {
	cerr << "  force_by_grape_on_all_leaves:  ";
	PRL(nj);
    }

    int n_pipes = g6_npipes_();
    hdyn **ilist = new hdynptr[n_pipes];

    // Compute particle forces, n_pipes at a time.

    bool status = false;

    int ip = 0;
    if (!cm) {
	for_all_leaves(hdyn, b, bb) {
	    ilist[ip++] = bb;
	    if (ip == n_pipes) {
		status |= force_by_grape(b->get_system_time(),
					 ilist, ip, nj, n_pipes,
					 true);	// "true" ==> compute pot only
		ip = 0;
	    }
	}
    } else {
	for_all_daughters(hdyn, b, bb) {
	    ilist[ip++] = bb;			// replicated code...
	    if (ip == n_pipes) {
		status |= force_by_grape(b->get_system_time(),
					 ilist, ip, nj, n_pipes,
					 true);	// "true" ==> compute pot only
		ip = 0;
	    }
	}
    }

    if (ip)
	status |= force_by_grape(b->get_system_time(),
				 ilist, ip, nj, n_pipes,
				 true);

    if (DEBUG) {
	cerr << "  ...leaving force_by_grape_on_all_leaves:  ";
	PRL(status);
    }

    delete [] ilist;
    return status;
}



//						*****************************
//						*****************************
//						***                       ***
//						***  The global function  ***
//						***                       ***
//						*****************************
//						*****************************


void grape_calculate_energies(hdyn *b,			// root node
			      real &epot,
			      real &ekin,
			      real &etot,
			      bool cm)			// default = false
{
    // Note: cm = true means use the CM approximation for all top-level
    // nodes.  If cm = false, then still use the CM approximation for
    // sufficiently close binaries to avoid roundoff errors on the GRAPE.

    if (DEBUG) {
	cerr << endl << "grape_calculate_energies..."
	     << endl << flush;
	PRL(cm);
    }

    if (!grape_is_open)
	reattach_grape(b->get_real_system_time(),
		       "grape_calculate_energies", b->get_kira_options());

    real cpu0 = cpu_time();

    real e_unpert;
    int nj =  send_all_leaves_to_grape(b, e_unpert, cm);



    PRC(nj); PRL(e_unpert); 



    if (force_by_grape_on_all_leaves(b, nj, cm)) {

	cerr << "grape_calculate_energies: "
	     << "error on return from force_by_grape_on_all_leaves()"
	     << endl;

	hw_err_exit("grape_calculate_energies", 1, b);
    }

    epot = ekin = 0;

    if (!cm) {

	real d_crit = critical_sep(b);

	for_all_leaves(hdyn, b, bb) {

	    // Logic here follows that in send_all_leaves_to_grape().

	    bool reset_bb = false;

	    while (use_cm_approx(bb, d_crit)) {
		bb = bb->get_parent();
		reset_bb = true;
	    }
	    if (reset_bb) {
		bb = bb->get_oldest_daughter()->get_younger_sister()
					      ->next_node(b);
		if (!bb) break;
	    }

	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = hdyn_something_relative_to_root(bb, &hdyn::get_vel);
	    ekin += 0.5*mi*vel*vel;
	}

	PRC(epot); PRL(ekin);

    } else {
	for_all_daughters(hdyn, b, bb) {
	    real mi = bb->get_mass();
	    epot += 0.5*mi*bb->get_pot();
	    vector vel = bb->get_vel();
	    ekin += 0.5*mi*vel*vel;
	}
    }

    etot = ekin + epot + e_unpert;

    grape_was_used_to_calculate_potential = true;	// trigger a reset
							// next time around

//    cerr << "CPU time for grape_calculate_energies() = "
//	 << cpu_time() - cpu0 << endl;

    if (DEBUG) {
	cerr << "...leaving grape_calculate_energies...  ";
	PRL(etot);
	cerr << endl;
    }
}

#endif

//*************************************************************************
//*************************************************************************



//  ***********************************************************************
//  *                                                                     *
//  *  grape_calculate_acc_and_jerk: Use the GRAPE hardware to compute    *
//  *				     accs and jerks on a list of nodes.   *
//  *                                                                     *
//  ***********************************************************************


local INLINE bool set_grape_neighbor_radius(hdyn * b, int nj_on_grape)

// Adjust the GRAPE neighbor radius to some reasonable value.

// Called by:	grape_calculate_acc_and_jerk()			// global

{
    static char *func = "set_grape_neighbor_radius";

    b->set_grape_rnb_sq(0.0);

    if (b->is_leaf()) {

	// For a single particle, since GRAPE-6 always computes the nn,
	// we only need to deal with neighbor lists if our radius is
	// nonzero (which we assume implies that all stellar radii are
	// nonzero and hence we must compute colls) and grape_nb_count
	// is zero.

	// Note that the nn from the GRAPE depends on distance only, so
	// is not necessarily as good as the criterion used elsewhere
	// in kira.  May be necessary also to use the neighbor lists in
	// cases where masses can differ greatly from the mean.
	// *** To be studied (Steve, 7/00). ***

	if (b->get_grape_nb_count() == 0 && b->get_radius() > 0) {

	    real rnb_sq = 4*b->get_radius()*b->get_radius();

	    // Probably not necessary to expand the neighbor lists in search
	    // of a coll.  Instead, confine the search to twice the stellar
	    // radius and use the nn if no coll is found.  Assume that one
	    // star or the other will pick up its collision partner in time
	    // for merger to occur.

#ifdef EXPAND_TO_FIND_COLL

	    // Old code:

	    // If the node has a valid nearest neighbor pointer, try
	    // to use the nn distance as well.

	    if (b->get_nn() && b->get_nn() != b
		&& b->get_d_nn_sq() < 0.1*VERY_LARGE_NUMBER) {

		rnb_sq = Starlab::max(rnb_sq, b->get_d_nn_sq());

	    } else {

		// Node does not know its nearest neighbor.
		// Note connections between d_min, r90, and rnn.

		real r90_sq = b->get_d_min_sq() / square(b->get_d_min_fac());
		real rnn_sq = r90_sq * pow((real)nj_on_grape, 4.0/3);

		// NB: rnn_sq ~ square of the average interparticle spacing.

		rnb_sq = Starlab::max(rnb_sq, rnn_sq);
	    }
#endif

	    // New code simply uses twice the radius (Steve, 5/02).

	    b->set_grape_rnb_sq(rnb_sq);
	}

    } else {

	// For a node, we will want to compute the perturbers, so we
	// need a larger value of grape_rnb_sq.  As with coll, only
	// bother to compute the list if grape_nb_count = 0, unless
	// we are forced to rebuild by an invalid perturber list.

	if (b->get_grape_nb_count() == 0 || !b->get_valid_perturbers()) {

	    // Old code:
	    //
	    // real rnb_sq = b->get_perturbation_radius_factor();

	    // If the node has a valid nearest neighbor pointer, use
	    // the nn distance as well.

	    // if (b->get_nn() && b->get_nn() != b
	    //	&& b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER)
	    //	rnb_sq = Starlab::max(rnb_sq, b->get_d_nn_sq());

	    // b->set_grape_rnb_sq(rnb_sq);

	    // New code (GRAPE-4 and GRAPE-6 versions; Steve, 12/01):

	    real r_pert2 = perturbation_scale_sq(b, b->get_gamma23());

//*****
	    // Limit r_pert2 to avoid GRAPE problems...

	    if (r_pert2 > 100) r_pert2 = 100;	// ***** arbitrary -- should
//*****						// ***** tie to system scale

	    // Note that r_pert2 may be very small for a tightly
	    // bound binary.

	    hdyn *nn = b->get_nn();
	    if (nn && nn != b && b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER)
		r_pert2 = Starlab::max(r_pert2, b->get_d_nn_sq());

	    // Note that this may cause overflow...

	    b->set_grape_rnb_sq(r_pert2);

#ifdef T_DEBUG
	    real sys_t = b->get_system_time();
	    if (IN_DEBUG_RANGE(sys_t)) {
		cerr << "DEBUG: " << func << ": computing perturbers for node "
		     << b->format_label() << endl << flush;
		PRC(binary_scale(b)); PRL(r_pert2);
		if (T_DEBUG_LEVEL > 0) {
		    PRC(b->get_gamma23()); PRL(b->get_d_nn_sq());
		    PRL(b->get_n_perturbers());
		    hdyn *od = b->get_oldest_daughter();
		    PRL(od->get_perturbation_squared());
		    PRL(od->get_pos());
		    PRL(od->get_pred_pos());
		}
	    }
#endif

	}
    }

    return (b->get_grape_rnb_sq() > 0);		// NB using rnb > 0 as a flag
}



local INLINE bool get_force_and_neighbors(xreal xtime,
					  hdyn *nodes[], int ni,
					  int nj_on_grape, int n_pipes,
					  bool need_neighbors,
					  int level = 0)

// Calculate the forces on the specified i-list of nodes, then read the
// GRAPE neighbor list if needed.   Deal with hardware neighbor list
// problems before returning.  Return true iff neighbor list overflow
// occurs.

// Called by:	get_coll_and_perturbers				// local
//		grape_calculate_acc_and_jerk			// global

{
    static char *func = "get_force_and_neighbors";

    if (ni <= 0) return false;

#ifdef T_DEBUG
    real sys_t = xtime;
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 1 << "  ";
	PRC(ni); PRC(nj_on_grape); PRL(n_pipes);
	PRC(need_neighbors); PRL(nodes[0]);
#if 0
	for (int i = 0; i < ni; i++) {
	    cerr << i << " " << nodes[i] << " "; PRL(nodes[i]->is_valid());
	    if (nodes[i]->is_valid())
		cerr << "    " << nodes[i]->format_label() << " "
		     << nodes[i]->get_pos() << endl
		     << "    " << nodes[i]->get_valid_perturbers() << " "
		     << nodes[i]->get_grape_rnb_sq() << endl << flush;
	}
#endif
    }
#endif

    if (force_by_grape(xtime, nodes, ni, nj_on_grape, n_pipes)) {

	// Hardware error has persisted despite NRETRY GRAPE reset(s).
	// Give up...

	hw_err_exit(func, 1, nodes[0]);
    }

    bool error = false;

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 2 << endl << flush;
    }
#endif

    if (need_neighbors) {

	// At least one i-particle needs coll or perturber information.
	// Bring all neighbor lists from the GRAPE to the front end.

      int status = 0;
#ifndef NO_G6_NEIGHBOUR_LIST
      status = g6_read_neighbour_list_(&cluster_id);
#endif

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 3 << endl << flush;
    }
#endif

	if (status) {

	    // An error has occurred.  Flag it and take appropriate action...

	    cerr << func << ":  error getting GRAPE neighbor data: ";
	    PRL(status);

	    if (status > 0) {

		// Hardware perturber lists on the GRAPE have overflowed.
		// Calling function must reduce neighbor radii and repeat
		// this group of particles.

		error = true;

	    } else {

		// An internal error has occurred -- do a hard
		// reset, repeat the force calculation and reread
		// the neighbor lists.

		// Make NRETRY attempts to reset to correct the problem
		// before giving up.

		cerr << "*** " << func << "(" << level
		     << "):  hardware error " << error << " at time "
		     << xtime << ",  ni = " << ni << endl;

		if (level < NRETRY) {

		    cerr << "Resetting GRAPE and retrying..." << endl;
		    reset_grape(nodes[0]->get_root());

		    error = get_force_and_neighbors(xtime, nodes, ni,
						    nj_on_grape, n_pipes,
						    need_neighbors,
						    level+1);
		} else

		    hw_err_exit(func, 2, nodes[0]);
	    }
	}
    }

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 4 << endl << flush;
    }
#endif

    return error;
}



local INLINE void swap(hdynptr ilist[], int i, int j)
{
    hdyn *tmp = ilist[i];
    ilist[i] = ilist[j];
    ilist[j] = tmp;
}

local INLINE int sort_nodes_and_reduce_rnb(hdynptr ilist[], int ni)

// Reorder the list of i-nodes to place those with rnb = 0 at the start.
// Reduce rnb for the remaining nodes and return the location on the
// new list of the first node with rnb > 0.

// Called by:	grape_calculate_acc_and_jerk()			// global

{
    static char *func = "sort_nodes_and_reduce_rnb";

    cerr << "In " << func << "() after neighbor-list overflow at time "
	 << ilist[0]->get_system_time() << endl;

    int inext = ni;

    // First move all the rnb>0 nodes to the end of the list.
    // The rnb>0 nodes will start at inext.

    int imax = ni;
    real rnb_max = 0;
    for (int i = ni-1; i >= 0; i--) {
	if (ilist[i]->get_grape_rnb_sq() > 0) {
	    if (i < --inext) swap(ilist, i, inext);
	    if (ilist[inext]->get_grape_rnb_sq() > rnb_max) {
		imax = inext;
		rnb_max = ilist[inext]->get_grape_rnb_sq();
	    }
	}	    
    }

    // ...then place the CM node with the biggest neighbor radius (a
    // guess at the node that caused the overflow) at inext...

    if (imax != inext) swap(ilist, inext, imax);

    // ...then adjust the neighbor radii of the others (starting at inext),
    // reducing radii for all leaves but only the *first* CM node found
    // (in location inext, by construction).

    for (int i = inext; i < ni; i++) {
	if (ilist[i]->is_leaf() || i == inext)
	    ilist[i]->set_grape_rnb_sq(0.5*ilist[i]->get_grape_rnb_sq());
    }

//    PRC(ni); PRL(inext);
    return inext;
}



// Neighbor list space is shared by get_neighbors_and_adjust_h2 (acc_and_jerk)
// and count_neighbors_and_adjust_h2 (densities) -- will probably be merged.

#define MAX_NEIGHBORS		(16*MAX_PERTURBERS)
#define RNB_INCREASE_FAC	2.0
#define RNB_DECREASE_FAC	0.8

static int max_neighbors = MAX_NEIGHBORS;
static int neighbor_list[MAX_NEIGHBORS];

// The reason for the large choice of MAX_NEIGHBORS here is to ensure
// that, even if the perturber list overflows, we still have room to
// store the full neighbor list, so the nn and coll pointers will still
// be correctly set.  The size of the neighbor sphere must be reduced
// if overflow occurs, or else the nn pointers will be unreliable.

//-------------------------------------------------------------------------

local INLINE int get_neighbors_and_adjust_h2(hdyn * b, int pipe)

// Set nn, coll, d_nn_sq and d_coll_sq, for particle b (pipe specified),
// and compute perturber lists for parent nodes.  Possible returns are:
//
//	0	All OK, no further action needed.
//	1	Array overflow, reduce grape_rnb_sq and retry.
//	2	No neighbors or no coll.
//		Old code: increase grape_rnb_sq and retry.
//		New code: accept NULL coll as is (and use nn).
//
// The value of grape_rnb_sq is changed here.  Retrying is up to the
// calling function.

// Called by:	get_coll_and_perturbers()			// local

{
    static char *func = "get_neighbors_and_adjust_h2";

    // Get the list of b's neighbors from the GRAPE, if available.

    int n_neighbors = 0;
    int status = 0;

#ifndef NO_G6_NEIGHBOUR_LIST
    status = g6_get_neighbour_list_(&cluster_id,
				    &pipe,
				    &max_neighbors,
				    &n_neighbors,
				    neighbor_list);
#endif

    if (status) {

	// GRAPE found too many neighbors:  n_neighbors >= max_neighbors,
	// so the kira perturber list will overflow.  Flag an error for
	// now (maybe unnecessary), and modify the neighbor radius.  We
	// reduce the size of the neighbor sphere to a point that will
	// still ensure too many perturbers, but which will not overflow
	// the neighbor list array.

	// Note from Steve (12/01): Looks like n_neighbors returns equal
	// to max_neighbors when overflow occurs, so we can't use the
	// value of n_neighbors to reduce grape_rnb_sq.  Instead, just
	// reduce the neighbor radius by a factor of 2.

	if (n_neighbors >= max_neighbors) {	// should be redundant...

	    // real rnb_fac = pow(max_neighbors/(real)n_neighbors, 0.6666667);
	    // b->set_grape_rnb_sq(rnb_fac*b->get_grape_rnb_sq());

	    b->set_grape_rnb_sq(0.25*b->get_grape_rnb_sq());

#if 0
	    cerr << func << ":  neighbor list overflow for "
		 << b->format_label() << " (pipe "
		 << pipe << ")" << endl;
	    cerr << "    at time " << b->get_system_time() << "  ";
	    PRC(n_neighbors); PRL(max_neighbors);
	    cerr << "    new rnb = " << sqrt(b->get_grape_rnb_sq())
		 << "  status = 1" << endl;
#else
#if 00
	    // Default short diagnostic message:
	    cerr << func << ": node " << b->format_label()
		 << " time " << b->get_system_time()
		 << " n_n " << n_neighbors
		 << endl << flush;
#endif
#endif

	    return 1;
	}
    }


#if 0000
    if (b->name_is("(5394,21337)") || b->name_is("(21337,5394)")) {
	cerr << func << ": "; PRL(b);
	PRI(4); PRC(b->get_grape_rnb_sq()); PRL(n_neighbors);
    }
#endif


    status = 2;

    if (n_neighbors > 0) {

	// Found at least one neighbor -- find the nearest and
	// determine the perturber list for a center-of-mass node.
	// Note that, if we get here, we recompute the nn pointer
	// for this particle, (possibly) overriding the result
	// returned by the GRAPE hardware.

	status = 0;

	real dmin_sq = VERY_LARGE_NUMBER;
	hdyn *bmin = NULL;
	real dcmin_sq = VERY_LARGE_NUMBER;
	hdyn *cmin = NULL;

	int npl = 0;
	hdyn **pl = NULL;
	real rpfac = 0;

	if (b->is_parent() && !b->get_valid_perturbers()) {

	    // Need to rebuild the perturber list.

	    if (b->get_oldest_daughter()->get_slow())
		clear_perturbers_slow_perturbed(b);

	    b->new_perturber_list();
	    pl = b->get_perturber_list();
	    rpfac = b->get_perturbation_radius_factor();
	}

	for (int j = 0; j < n_neighbors; j++) {

	    hdyn *bb = node_list[neighbor_list[j]];

	    // bb is the j-th neighbor of b (list not ordered).

	    if (bb != b) {		// bb = b shouldn't occur...

		vector diff = b->get_pred_pos() - bb->get_pred_pos();
		real d2 = diff * diff;

		// (Re)compute nn and coll here.

		real sum_of_radii = get_sum_of_radii(b, bb);
		update_nn_coll(b, 100,		// (100 = ID)	    // inlined
			       d2, bb, dmin_sq, bmin,
			       sum_of_radii,
			       dcmin_sq, cmin);

		// Recompute the perturber list for parent nodes.
		// See equivalent code for use without GRAPE in
		// hdyn_ev.C/flat_calculate_acc_and_jerk.

		if (b->is_parent() && !b->get_valid_perturbers()) {

		    // Update the new perturber list.

		    if (is_perturber(b, bb->get_mass(),
				     d2, rpfac)) {		    // inlined

			if (npl < MAX_PERTURBERS)
			    pl[npl] = bb;

			npl++;
		    }
		}
	    }
	}

	if (b->is_parent() && !b->get_valid_perturbers()) {

	    // Found npl perturbers.

	    b->set_n_perturbers(npl);

	    if (npl > MAX_PERTURBERS) {

		// Too many perturbers.

		b->remove_perturber_list();

//		cerr << func << ": too many perturbers for "
//		     << b->format_label()
//		     << " -- deleting list" << endl << flush;

	    } else {

		b->set_valid_perturbers(true);
//		cerr << "valid_perturbers = true #1 for "
//		     << b->format_label()
//		     << endl;
	    }

	}

	if (bmin) {
	    b->set_nn(bmin);				// overwrite nn
	    b->set_d_nn_sq(dmin_sq);


#if 0000
	    if (b->name_is("(5394,21337)") || b->name_is("(21337,5394)")) {
		cerr << "get_nbrs:  ";
		PRC(b); PRC(bmin); PRL(dmin_sq);
	    }
#endif


	} else
	    status = 2;

	if (cmin) {
	    b->set_coll(cmin);
	    b->set_d_coll_sq(dcmin_sq);
	} else
	    status = 2;

    } else {

	// No neighbors found.  In new code this is OK -- just modify
	// the coll pointer appropriately.

	if (b->is_parent()) {

	    // No neighbors OK in this case.  Set up valid zero-length
	    // perturber list and use hardware nearest neighbor as both
	    // nn and coll.

	    b->new_perturber_list();
	    b->set_valid_perturbers(true);
	    b->set_n_perturbers(0);

//	    cerr << "valid_perturbers = true #2 for "
//	         << b->format_label()
//	         << endl;

	}

	// Use nn as coll (always do this in the new version).

#ifdef EXPAND_TO_FIND_COLL
	if (b->is_parent())
#endif
	{
	    if (b->get_nn()) {
		b->set_d_coll_sq(b->get_d_nn_sq());
		b->set_coll(b->get_nn());
	    } else
		b->set_nn(b);
	}

    }


#if 0000
    if (b->name_is("(5394,21337)") || b->name_is("(21337,5394)")) {
	cerr << func << ": ";
	PRC(b->get_grape_rnb_sq()); PRL(n_neighbors);
    }
#endif

#ifdef EXPAND_TO_FIND_COLL

    // Old code:
    //
    // If no nearest neighbor or coll is found, enlarge the neighbor-sphere
    // radius and try again.

    if (status)
	b->set_grape_rnb_sq(RNB_INCREASE_FAC*b->get_grape_rnb_sq());

#endif

    // New code: no neighbors OK in all cases.

    return status;	// note: status indicates success of nn/coll search,
			//	 not the validity of the neighbor list
}



#define MAX_FORCE_COUNT	20

local INLINE int get_coll_and_perturbers(xreal xtime,
					 hdynptr *ilist, int ni,
					 real h2_crit,
					 int nj_on_grape, int n_pipes)

// Compute the colls, nns, and perturber lists for those nodes on the
// i-list that require them.  Return the location following the last
// node successfully treated.

// Called by:	grape_calculate_acc_and_jerk()			// global

{
    static char *func = "get_coll_and_perturbers";

    int inext = 0;

//    cerr << func << ": "; PRL(ni);

    for (int ip = 0; ip < ni; ip++) {

	int pipe = ip;
	hdyn *bb = ilist[ip];

	if (bb->get_grape_rnb_sq() > 0) {

//	    if (bb->is_parent()) {
//		PRI(4); PRC(ip);
//		cerr << bb->format_label() << " "
//		     << bb->get_grape_rnb_sq() << " "
//		     << bb->get_valid_perturbers() << endl;
//	    }

	    // Note: for a CM node, only reach this point if we want to
	    // recompute the neighbor lists.  Nodes with valid lists
	    // will use the hardware nn and skip this section.

	    int count_force = 1;
	    int status;

	    while ((status = get_neighbors_and_adjust_h2(bb, pipe))) {

//		if (bb->is_parent()) {
//		    PRC(bb->get_time()); PRC(bb->format_label()); PRL(status);
//		    PRC(count_force); PRL(bb->get_grape_rnb_sq());
//		}

		// Neighbor list must be recomputed:
		//
		//	status = 1  ==>  too many neighbors; decrease radius
		//	status = 2  ==>  too few neighbors
		//			 old: increase radius
		//			 new: accept as is
		//
		// The value of grape_rnb_sq has already been adjusted if
		// necessary, and nn has been set from the GRAPE in either
		// case.

//		if (bb->is_parent()) PRL(status);

		if (status == 2
#ifdef EXPAND_TO_FIND_COLL
		    && bb->is_parent()
#endif
		    ) {

		    // Didn't find a neighbor, but we started our search at
		    // the perturbation radius or twice the stellar radius.
		    // Use the hardware neighbor as coll and assume zero
		    // perturbers with no perturbation -- already set in
		    // get_neighbors_and_adjust_h2(), so nothing to do here.

		    break;				// go on to next ip

		} else if (count_force > MAX_FORCE_COUNT
			   || (status == 2
			       && bb->get_grape_rnb_sq() > h2_crit)) {

		    // Give up -- can't find a neighbor.  Flag with nn = bb
		    // (not really necessary, but probably harmless).

		    bb->set_nn(bb);			// kira checks
							// for nn = bb
		    bb->set_d_nn_sq(2*h2_crit);

#if 0000
		    if (bb->name_is("(5394,21337)")
			 || bb->name_is("(21337,5394)")) {
			cerr << func << ":  setting nn = bb for ";
			PRC(bb); PRL(bb->get_d_nn_sq());
			PRC(status); PRL(count_force);
		    }
#endif

		    break;				// go on to next ip

		} else {

		    // We have modified the neighbor sphere and must
		    // recompute the force to get the new neighbor list.
		    // Currently, we simply iterate until the neighbors
		    // for this particle are found (or we give up), then
		    // exit, returning the next ip on the list.

		    // Reducing rnb for a CM node means that we
		    // can't construct a valid perturber list from
		    // the neighbor information.

		    if (status == 1 && bb->is_parent()) {

			// bb->set_valid_perturbers(false);	// already set

//			cerr << "valid_perturbers = false, #1 for "
//			     << bb->format_label()
//			     << endl;
		    }

#if 00000
		    if (bb->name_is("(5394,21337)")
			 || bb->name_is("(21337,5394)")) {
			cerr << func << ": recomputing force for ";
			PRC(ip); PRL(bb);
			PRL(ilist[ip]);
		    }
#endif

		    // Recompute just this particle -- it will go in pipe 0.

		    pipe = 0;
		    ip = ni;			// force exit from "for" loop

		    while (get_force_and_neighbors(xtime,
						   &bb,
						   1,
						   nj_on_grape, n_pipes,
						   true)) {

			// The GRAPE neighbor list has overflowed; we must
			// decrease the neighbor radius.  Jun says that
			// there is no useful information in the returned
			// number of neighbors, so just reduce the radius
			// by a fixed factor.  Note that the reduction
			// factor is much closer to 1 than the expansion
			// factor used in get_neighbors_and_adjust_h2(),
			// so we shouldn't run into an infinite loop.
			// Nevertheless, cautiously count iterations
			// and impose a (generous) upper limit.

			bb->set_grape_rnb_sq(RNB_DECREASE_FAC
					     * bb->get_grape_rnb_sq());
			count_force++;

			if (bb->is_parent()) {

			    // bb->set_valid_perturbers(false);	// already set

//			    cerr << "valid_perturbers = false, #2 for "
//				 << bb->format_label()
//				 << endl;
			}

		    }

		}		// if (count_force...) {} else

#if 0000
		    if (bb->name_is("(5394,21337)")
			 || bb->name_is("(21337,5394)")) {
			cerr << func << ": repeating while loop with ";
			PRL(bb->get_grape_rnb_sq());
		    }
#endif

	    }			// while ((status = get_...))

	}			// if (bb->get_grape_rnb_sq()...)

	inext++;
    }				// for (ip...)

    return inext;
}



//						*****************************
//						*****************************
//						***                       ***
//						***  The global function  ***
//						***                       ***
//						*****************************
//						*****************************


void grape_calculate_acc_and_jerk(hdyn **next_nodes,
				  int n_next,
				  xreal xtime,
				  bool restart)

//  This function is called from kira_calculate_top_level_acc_and_jerk,
//  which is called only from calculate_acc_and_jerk_for_list.

{
    static char *func = "grape_calculate_acc_and_jerk";

    static int n_pipes = 0;		// number of pipelines to use

    static int nj_on_grape;		// current number of j-particles
					// in GRAPE memory

    if (DEBUG) {
	cerr << endl << func << "..." << endl << flush;
    }

    if (n_pipes == 0) n_pipes = g6_npipes_();

    if (n_next <= 0) return;
    hdyn *root = next_nodes[0]->get_root(); 
    kira_options *ko = root->get_kira_options();

    int coll_freq = Starlab::max(1, ko->grape_coll_freq);
    int pert_freq = Starlab::max(1, ko->grape_pert_freq);

    //------------------------------------------------------------------

    // Test the state of the GRAPE and open it if necessary.
    // If restart is true, we must reinitialize the GRAPE interface
    // after a change in the tree or other kira configuration.
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

	reattach_grape((real)xtime, func, ko);

	if (!restart) grape_reattached = true;

	// Restart irrespective of the actual argument.

	restart = true;
    }

    if (grape_was_used_to_calculate_potential) {
	restart = true;
	grape_was_used_to_calculate_potential = false;
    }

    if (restart) {
//	cerr << "restarting GRAPE at time " << xtime << endl << flush;
	nj_on_grape = initialize_grape_arrays(root, !grape_reattached);
	n_previous_nodes = 0;
    }

    //------------------------------------------------------------------

    // Store the particles in the previous block to GRAPE memory
    // (i.e. update GRAPE for the previous step).

    for (int i = 0; i < n_previous_nodes; i++)
	send_j_node_to_grape(previous_nodes[i]);

    //------------------------------------------------------------------

    int i, n_top;

    // Create the list of top-level nodes in the present block step.

    for (n_top = i = 0; i < n_next; i++)
	if (next_nodes[i]->is_top_level_node())
	    current_nodes[n_top++] = next_nodes[i];

    // Now n_top is the number of top-level nodes in the current list.

    // Initialize neighbor and perturber information and save a copy of
    // the current_nodes list as previous_nodes (for use next time around).

    n_previous_nodes = n_top;

//    int n_needpertlist = 0;

#ifdef T_DEBUG
    real sys_t = root->get_system_time();
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 1 << endl << flush;
    }
#endif

    // See which nodes need neighbors from the GRAPE and optionally sort
    // the list to place all such nodes at the start.
    //					(Experimental - Steve, 5/02)

    int last_nbr_node = -1;

    for (i = 0; i < n_top; i++) {

	hdyn *bb = previous_nodes[i] = current_nodes[i];

	// Set a reasonable h2 value for this node.  Valid_perturbers
	// is used in set_grape_neighbor_radius() to force a neighbor
	// computation for a node that doesn't have a valid list.

	bool need_nbrs = set_grape_neighbor_radius(bb, nj_on_grape);

	// From here on, use valid_perturbers as a flag to indicate
	// whether the hardware neighbor list should be used to construct
	// a perturber list.

	if (need_nbrs && bb->is_parent()) {
	    bb->set_valid_perturbers(false);
//	    n_needpertlist++;
	}

#ifdef SORT_NODE_LIST
	if (need_nbrs) {
	    hdyn *tmp = current_nodes[++last_nbr_node];
	    current_nodes[last_nbr_node] = current_nodes[i];
	    current_nodes[i] = tmp;
	}
#endif

    }

    if (DEBUG) {
	cerr << func << ":  ";
	PRC(xtime); PRC(n_next); PRL(n_top);
    }


//    PRC(n_top); PRL(n_needpertlist);


#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 2 << endl << flush;
    }
#endif

    //------------------------------------------------------------------

    // Calculate the force on the current_nodes list.

    real h2_crit = 8192 * current_nodes[0]->get_d_min_sq();

    // We will stop expanding the GRAPE neighbor sphere (in those cases
    // where expansion is indicated -- finding nearest neighbors and colls)
    // once its size exceeds the critical value h2_crit.  However, it is
    // OK to set grape_rnb_sq greater than h2_crit -- the neighbor sphere
    // then simply won't be expanded if no neighbors are found.

    // *** Should contain a factor of ~(m_max/<m>)^2, but not crucial...

    // Note:  for equal-mass systems and standard units, this critical
    //	      radius is less than the interparticle spacing for N > ~1000.

    // Compute the i-forces in chunks of maximum size n_pipes.
    // Compute need_neighbors separately for each chunk.

    bool need_neighbors = false;

    i = 0;
    while (i < n_top) {

	int ni = Starlab::min(n_pipes, n_top - i);

	need_neighbors = false;
	for (int ip = 0; ip < ni; ip++)
	    if (current_nodes[i+ip]->get_grape_rnb_sq() > 0) {
		need_neighbors = true;
		break;
	    }

	if (DEBUG > 1) {
	    PRI(2); PRC(i); PRC(ni); PRL(need_neighbors);
	}

	// Get the forces on particles i through i + ni - 1 and
	// read the GRAPE neighbor list, if necessary.  Function
	// get_force_and_neighbors (via force_by_grape) also sets
	// nn pointers.

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 3 << endl << flush;
    }
#endif

	if (get_force_and_neighbors(xtime, current_nodes + i, ni,
				    nj_on_grape, n_pipes, need_neighbors)) {

	    // Hardware neighbor list overflow.  Restructure the list,
	    // reduce neighbor radii, and retry starting with those nodes
	    // for which colls are actually needed.  Forces and nns are
	    // OK at this point, but neighbor lists are not.

	    cerr << "  HW nbr overflow, "; PRL(i);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 4 << endl << flush;
    }
#endif

	    i += sort_nodes_and_reduce_rnb(current_nodes+i, ni);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 5 << endl << flush;
    }
#endif

	} else if (need_neighbors) {

	    // Neighbor lists are OK.  Get colls and perturber lists.

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 6 << endl << flush;
    }
#endif

	    i += get_coll_and_perturbers(xtime, current_nodes+i, ni,
					 h2_crit, nj_on_grape, n_pipes);

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 7 << endl << flush;
    }
#endif

	} else

	    i += ni;

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 8 << endl << flush;
    }
#endif

    }

    //------------------------------------------------------------------

    // Update the grape_nb_count flags.

//    int n_gotpertlist = 0;

#ifdef T_DEBUG
    if (IN_DEBUG_RANGE(sys_t)) {
	cerr << "DEBUG: " << func << " " << 9 << endl << flush;
    }
#endif

    for (i = 0; i < n_top ; i++) {

	hdyn *bb = current_nodes[i];

//	n_gotpertlist += bb->get_valid_perturbers();


#if 0000
	if (bb->name_is("(5394,21337)")
	    || bb->name_is("(21337,5394)"))
	    cerr << "leaving " << func << endl << endl << flush;
#endif


	// Frequency of coll checks is every grape_coll_freq force calculation.
	// Frequency of pert checks is every grape_pert_freq force calculation.
	// Defaults are set in kira_defaults.h.

	if (bb->is_leaf()) {

	    bb->set_grape_nb_count((bb->get_grape_nb_count() + 1) % coll_freq);

	    // Shouldn't be needed, but just in case:

	    if (bb->get_coll() == NULL) bb->set_coll(bb->get_nn());

	} else {

//	    *** If we reduce the frequency of perturber checks, then we
//	    *** must be sure to restore the CMs on the perturber list,
//	    *** as the correction to the CM force depends on it...
//	    ***
//	    *** Some care is needed if we do reduce the frequency, as
//	    *** nodes may vanish or merge.  (However, the list contains
//	    *** only single stars after correction, so CM changes
//	    *** shouldn't be a problem).
//	    ***						(Steve, 6/01)

	    bb->set_grape_nb_count((bb->get_grape_nb_count() + 1) % pert_freq);

//	    (Precise value of pert_freq doesn't seem to have a large effect
//	    on overall CPU time -- Steve, 5/02)

	}
    }

//    PRL(n_gotpertlist);

    if (DEBUG) {
	cerr << "...leaving " << func << endl << endl << flush;
    }
}



//  **********************************************************************
//  *                                                                    *
//  * grape_calculate_densities:  Determine particle densities, giving	 *
//  *				  zero density to particles with no	 *
//  *				  neighbor within sqrt(h2_crit).	 *
//  *                                                                    *
//  **********************************************************************

//  **********************************************************************
//  *****  NOTE that the "new" density algorithm from hdyn_grape4.C  *****
//  *****  should be incorporated here too.                          *****
//  **********************************************************************

local INLINE void set_grape_density_radius(hdyn *b, real rnn_sq)

// For a single particle, try to adjust the initial radius so that
// it will contain just ~10-20 neighbors (set r = 3*d_nn, if known).

// Called by:	grape_calculate_densities()			// global

{
    if (b->get_nn() != NULL && b->get_nn() != b
	&& b->get_d_nn_sq() < 0.1* VERY_LARGE_NUMBER
	&& b->get_d_nn_sq() > 0) {

	// Node seems to have a valid nearest neighbor pointer.

	// Old code:
	//
	// b->set_grape_rnb_sq(9 * b->get_d_nn_sq());

	// New code (from GRAPE-4 version):

	// Modify the initial guess for particles with a close nn.

	if (DEBUG) {
	    cerr << "nn OK for " << b->format_label() << ",  ";
	    PRC(rnn_sq); PRL(sqrt(b->get_d_nn_sq()));
	    PRL(b->get_nn()->format_label());
	}

#if 0
	if (b->get_d_nn_sq() < rnn_sq)
	    rnn_sq = sqrt(rnn_sq * b->get_d_nn_sq());	// empirical
	else
	    rnn_sq = b->get_d_nn_sq();
#endif

	rnn_sq = 5*b->get_d_nn_sq();

    } else {

	// Node does not know its nearest neighbor.  Value of d_nn_sq
	// is where the neighbor search stopped.

	// Old:
	//
	// b->set_grape_rnb_sq(9 * pow(b->get_d_min_sq(), 1.0/3.0));  // (??)

	rnn_sq = 4*Starlab::max(rnn_sq, b->get_d_nn_sq());

	if (DEBUG)
	    cerr << "no nn for " << b->format_label() << ",  ";
    }

    b->set_grape_rnb_sq(rnn_sq);
    if (DEBUG) PRL(sqrt(b->get_grape_rnb_sq()));
}



local void check_neighbors(hdyn *b, real rnb_sq, int indent = 0)

// (The hard way...)

{
    int nn = 0;
    b = b->get_top_level_node();
    for_all_daughters(hdyn, b->get_root(), bb)
	if (bb != b)
	    if (square(bb->get_pos() - b->get_pos()) < rnb_sq) nn++;
    PRI(indent);
    cerr << "check:  " << nn << " neighbors within " << sqrt(rnb_sq)
	 << " of " << b->format_label() << endl;
}



#define N_DENS	12	// density is based on the 12th nearest neighbor

local INLINE bool count_neighbors_and_adjust_h2(hdyn * b, int pipe)

// Determine density from the neighbor list, or increase the radius
// of the neighbor sphere.

// Called by:	grape_calculate_densities()			// global

{
    // Get the list of neighbors from the GRAPE.

    int n_neighbors = 0;
    int status = 0;
#ifndef NO_G6_NEIGHBOUR_LIST
    status = g6_get_neighbour_list_(&cluster_id,
				    &pipe,
				    &max_neighbors,
				    &n_neighbors,
				    neighbor_list);
#endif

//    cerr << "  GRAPE "; PRC(max_neighbors); PRL(n_neighbors);
//    check_neighbors(b, b->get_grape_rnb_sq(), 2);

    if (status) {

	// GRAPE has found too many neighbors:  n_neighbors >= max_neighbors.
	// Attempt to reduce the size of the neighbor sphere to contain
	// (say) 10*N_DENS stars.  Note (Steve, 12/01) that n_neighbors
	// returns equal to max_neighbors in case of overflow.

	if (n_neighbors >= max_neighbors) {	// should be redundant...

	    real rnb_fac = pow(10*N_DENS/(real)n_neighbors, 0.6666667);
	    b->set_grape_rnb_sq(rnb_fac*b->get_grape_rnb_sq());

	    cerr << "count_neighbors_and_adjust_h2: "
		 << "neighbor list overflow for "
		 << b->format_label() << " (pipe "
		 << pipe << ")" << endl;
	    PRC(n_neighbors); PRC(max_neighbors);
	    cerr << "new grape_rnb = " << sqrt(b->get_grape_rnb_sq())
		 << endl;

	    return false;
	}
    }

    if (n_neighbors < N_DENS) {

	// Too few neighbors.  Try again.

	if (DEBUG > 1) {
	    cerr << "  increasing grape_rnb_sq for " << b->format_label()
		 << " (n_neighbors = " << n_neighbors << ", grape_rnb = "
		 << sqrt(b->get_grape_rnb_sq()) << ")" << endl;
	}

	// Old:
	//
	// real fac = 2;
	// if (n_neighbors < 4) fac = 4;

#if 0
	real fac = 1.25;
	if (n_neighbors < 12) fac *= 1.25;
	if (n_neighbors < 8) fac *= 1.6;
	if (n_neighbors < 4) fac *= 1.6;
	if (n_neighbors < 2) fac *= 1.6;
#endif

	real fac = Starlab::min(2.0, pow((N_DENS+3.0)/(1.0+n_neighbors), 1.5));

	b->set_grape_rnb_sq(fac * b->get_grape_rnb_sq());
	return false;
    }

    // Make a list of nodes to send to compute_density().

    dyn **dynlist = new dynptr[n_neighbors];

    real d_max = 0;
    for (int j = 0; j < n_neighbors; j++) {

	hdyn *bb = node_list[neighbor_list[j]];
	dynlist[j] = (dyn*)bb;

	// bb is the j-th neighbor of b.

	if (DEBUG > 1) {
	    if (bb != b) {
		vector diff = b->get_pred_pos() - bb->get_pred_pos();
		real d2 = diff * diff;
		d_max = Starlab::max(d_max, d2);
	    }
	}
    }

    if (DEBUG > 1) {
	real grape_rnb = sqrt(b->get_grape_rnb_sq());
	d_max = sqrt(d_max);
	cerr << "  " << b->format_label() << ": ";
	PRC(n_neighbors), PRC(grape_rnb), PRL(d_max);
    }

    compute_density(b, N_DENS, dynlist, n_neighbors);	// writes to dyn story

    // Strangely, compute_density sometimes fails to write the proper
    // info to the dyn story.  Reason and circumstances still unknown...

    if (DEBUG > 1) {
	if (find_qmatch(b->get_dyn_story(), "density_time")
	    && find_qmatch(b->get_dyn_story(), "density_k_level")
	    && find_qmatch(b->get_dyn_story(), "density")) {

	    real density_time = getrq(b->get_dyn_story(), "density_time");
	    int  density_k    = getiq(b->get_dyn_story(), "density_k_level");
	    real density      = getrq(b->get_dyn_story(), "density");

	    PRI(4); cerr << b->format_label() << ": ";
	    PRC(density_time); PRC(density_k); PRL(density);
	    PRI(4); PRL(b->get_pos());
	    check_neighbors(b, b->get_grape_rnb_sq(), 4);

	} else {

	    PRI(4); cerr << "density data unavailable for "
	                 << b->format_label() << endl;
	}
    }

    delete [] dynlist;
    return true;
}



local INLINE bool get_densities(xreal xtime, hdyn *nodes[],
				int ni, real h2_crit,
				int nj_on_grape, int n_pipes)
{
    // Compute the densities of the ni particles in nodes[].
    // Return true iff an error occurred and densities could
    // not be determined.

    static char *func = "get_densities";
    bool error = false;

    if (DEBUG)
	cerr << func << "..." << nodes[0]->format_label()
	     << "  " << nodes[0]->get_grape_rnb_sq() << endl << flush;

    if (ni < 0) return error;			// should never happen

    // Get the forces and neighbor lists for the current group of particles.

    int status = 0;
    if (status = get_force_and_neighbors(xtime, nodes, ni,
					 nj_on_grape, n_pipes, true)) {

	// Neighbor list overflow.  Return with error = true will cause
	// the calling function to reduce neighbor radii and retry.

        cerr << func << ": " << "error 1 getting GRAPE neighbor data: "; 
	PRL(status);

	error = true;

    } else {

	// Determine densities for the present block of particles.

        for (int ip = 0; ip < ni; ip++) {

	    int n_retry = 0;
	    hdyn *bb = nodes[ip];

	    if (bb->get_grape_rnb_sq() > 0) {	    // neighbors still needed

		while (!count_neighbors_and_adjust_h2(bb, ip)) {

		    // (Function sets density on success.)

		    if (bb->get_grape_rnb_sq() > h2_crit) {

			// Write zero density to the dyn story.

			putrq(bb->get_dyn_story(), "density_time",
			      (real)bb->get_system_time());
			putrq(bb->get_dyn_story(), "density", 0.0);

			if (DEBUG > 1) {
			    PRI(2); PR(bb->get_grape_rnb_sq());
			    cerr << " too large for "
			         << bb->format_label() << endl;
			    PRI(2); PRL(bb->get_pos());
			    check_neighbors(bb, bb->get_grape_rnb_sq(), 2);
			}

			break;

		    } else {

			// Changed the neighbor sphere size; recompute all
			// forces and neighbor lists (probably overkill).

			// If this becomes an issue, can do better by setting
			// neighbor radii for successfully handled particles
			// to zero.  Also, we could restart after those
			// particles, as in grape_calculate_acc_and_jerk().
			//					[Steve, 6/01]

			if (++n_retry > 20) 
			    hw_err_exit(func, 2, nodes[0]);

			if (DEBUG > 1
			    || (DEBUG > 0 && n_retry > 4 && n_retry%5 == 0)
			    || (n_retry > 9 && n_retry%10 == 0) ) {
			    PRI(2); cerr << func << ": recomputing forces for "
				<< nodes[0]->format_label() << " + " << ni-1
				<< endl;
			    cerr << "                 first rnb_sq = "
				 << nodes[0]->get_grape_rnb_sq()
				 << ",  n_retry = " << n_retry << endl;
			}

			int status = 0;
			if (status = get_force_and_neighbors(xtime, nodes, ni,
						nj_on_grape, n_pipes, true)) {

			    cerr << func << ": "
				 << "error 2 getting GRAPE neighbor data: "; 
			    PRL(status);

			    error = true;

			    ip = ni;	// force exit from for loop
			    break;
			}
		    }
		}

		if (!error)
		    bb->set_grape_rnb_sq(0);	// don't recalculate density
	    }
	}
    }

    if (DEBUG) {
	cerr << "leaving " << func << "...";
	PRL(error);
    }

    return error;
}



//						*****************************
//						*****************************
//						***                       ***
//						***  The global function  ***
//						***                       ***
//						*****************************
//						*****************************


void grape_calculate_densities(hdyn* b,			// root node
			       real h2_crit)		// default = 4
{
    static char *func = "grape_calculate_densities";

    if (DEBUG) {
	cerr << endl << func << "..."; PRL(h2_crit);
    }

#ifdef NO_G6_NEIGHBOUR_LIST
    cerr << "No hardware neighbor list available..." << endl;
    return;
#endif

    static int max_neighbors = MAX_PERTURBERS;
    static int neighbor_list[MAX_PERTURBERS];

    if (!grape_is_open)
	reattach_grape(b->get_real_system_time(), func, b->get_kira_options());

    // Copy all top-level nodes to the GRAPE hardware.
    // We will compute the forces, etc. using the same functions
    // as grape_calculate_acc_and_jerk.

    int nj_on_grape = initialize_grape_arrays(b);

    // Make a list of top-level nodes.

    hdyn **top_nodes = new hdynptr[nj_on_grape];

    int n_top = 0;
    for_all_daughters(hdyn, b, bb)
	top_nodes[n_top++] = bb;

    // Set h2 values.

    // New code from GRAPE-4 version:

    // Determine rnn_sq ~ square of the average interparticle spacing,
    // for use as a scale in setting h2.  Note the connections between
    // d_min, r90, and rnn.

    real r90_sq = b->get_d_min_sq() / square(b->get_d_min_fac());
    real rnn_sq = r90_sq * pow((real)nj_on_grape, 4.0/3);

    for (int j = 0; j < n_top; j++)
//	set_grape_density_radius(top_nodes[j], h2_crit);
	set_grape_density_radius(top_nodes[j], rnn_sq);

    int n_pipes = g6_npipes_();

    int n_retry = 0;
    int count = 0;

    // Compute the densities in chunks of maximum size n_pipes.

    int i = 0;

    while (i < n_top) {

	int inext = i;
	int ni = Starlab::min(n_pipes, n_top - i);

#ifdef SPZ_GRAPE6

	// Use one pipeline only to (try to) avoid neighbor list overflow...
	// Implemented by SPZ on May 8 2001.
	// May no longer be necessary (Steve, 6/01).

	ni = 1;
#endif

	// Get the forces on particles i through i + ni - 1, determine
	// their neighbors and hence their densities.

	// Code follows that in grape_calculate_acc_and_jerk.
	// May be possible to combine the two...

	if (get_densities(b->get_system_time(),
			  top_nodes + i, ni, h2_crit,
			  nj_on_grape, n_pipes)) {

	    // The neighbor list overflowed.  Reduce all current neighbor
	    // radii and try again.

	    if (DEBUG) {
	        cerr << "reducing neighbor radii: i = " << i
		     << ", first rnb_sq = " << top_nodes[i]->get_grape_rnb_sq()
		     << endl;
	    }

	    // Reduction factor is 0.9 for now...

	    for (int j = i; j < i + ni; j++) {
	        hdyn *bb = top_nodes[j];
		bb->set_grape_rnb_sq(0.9 * bb->get_grape_rnb_sq());
	    }

	    n_retry++;
	    if (++count > 20) {

	        // Too many retries.  Quit in confusion...

	        cerr << func << ": "
		     << "error getting GRAPE neighbor data: " << endl;

		hw_err_exit(func, 1, b);
	    }

	} else {

	    // Move on to the next block of particles.

	    i += ni;
	    count = 0;
	}
    }

    // Timestamp the root node.

    putrq(b->get_dyn_story(), "density_time", (real)b->get_system_time());

    if (n_retry > 10) {
	cerr << func << ":  ";
	PRL(n_retry);
    }

    if (DEBUG)
        cerr << "...leaving " << func << endl << endl << flush;

    // Force cleanup later.

    grape_was_used_to_calculate_potential = true;
    delete [] top_nodes;
}



//  **********************************************************************
//  *                                                                    *
//  * External cleanup -- delete local static arrays.                    *
//  *                                                                    *
//  **********************************************************************

void clean_up_hdyn_grape()

// Explicitly delete static local data.

{
    // From initialize_grape_arrays()/grape_calculate_acc_and_jerk():

    if (node_list) delete [] node_list;
    if (current_nodes) delete [] current_nodes;
    if (previous_nodes) delete [] previous_nodes;

    // From force_by_grape():

    if (iindex) delete [] iindex;
    if (ipos) delete [] ipos;
    if (ivel) delete [] ivel;
    if (iacc) delete [] iacc;
    if (ijerk) delete [] ijerk;
    if (ipot) delete [] ipot;
    if (ih2) delete [] ih2;
    if (inn) delete [] inn;
}
