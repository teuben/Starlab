
//// make_tree:  hierarchically decompose a flat tree into substructure.
////
//// Options:    -b    use binary decomposition only [false]
////             -d    require that a tuple be bound [false]
////             -D    specify debug level [0]
////             -k    specify maximum k-tuple sought [no maximum]
////             -o    pipe restructured system to cout [true]
////             -s    require stability of tuple at periastron [false]
////
//// Example:    mkscat -t -p | make_tree -D 1 | pretty_print_tree

#include "sdyn.h"

// pp: Recursively pretty-print a node:

local void pp(sdyn* b, ostream & s, int level = 0) {

    s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    b->pretty_print_node(s);
    s << " \t"<< b->get_mass() << " \t"
      << b->get_pos() << "   "
      << b->get_vel() <<endl;

    for_all_daughters(sdyn, b, daughter) pp(daughter, s, level + 1);	
}

#ifndef TOOLBOX

static char id_string[50];

// id: Return a node's name or number as a character string.

char* id(sdyn* bb)
{
    if (bb->get_name() != NULL)
	return bb->get_name();
    else {
	sprintf(id_string, "%d", bb->get_index());
	return id_string;
    }
    return NULL;
}


// count_daughters: Return the number of daughters of the specified node.

local int count_daughters(sdyn* b)
{
    int n = 0;
    for_all_daughters(sdyn, b, bb) n++;

    return n;
}


// find_neighbors: Determine near-neighbors and next-nearest neighbors of
//		   all particles.

local void find_neighbors(sdyn* b)
{
    for_all_daughters(sdyn, b, bb) {
	bb->set_nn_dr2(VERY_LARGE_NUMBER);
	bb->set_nn_ptr(NULL);
	bb->set_nnn_dr2(VERY_LARGE_NUMBER);
	bb->set_nnn_ptr(NULL);
    }

    for_all_daughters(sdyn, b, bi)
	for (sdyn* bj = bi->get_younger_sister(); 
	     bj != NULL; bj = bj->get_younger_sister()) {

	    real rij2 = square(bi->get_pos() - bj->get_pos());

	    if (rij2 < bi->get_nn_dr2()) {
		bi->set_nn_dr2(rij2);
		bi->set_nn_ptr(bj);
	    } else if (rij2 < bi->get_nnn_dr2()) {
		bi->set_nnn_dr2(rij2);
		bi->set_nnn_ptr(bj);
	    }

	    if (rij2 < bj->get_nn_dr2()) {
		bj->set_nn_dr2(rij2);
		bj->set_nn_ptr(bi);
	    } else if (rij2 < bj->get_nnn_dr2()) {
		bj->set_nnn_dr2(rij2);
		bj->set_nnn_ptr(bi);
	    }
	}
}


// detach_from_tree: Cut all links between the specified node and the tree,
//		     correct the tree structure accordingly.

local void detach_from_tree(sdyn* bi)
{
    sdyn* parent = bi->get_parent();
    sdyn* elder = bi->get_elder_sister();
    sdyn* younger = bi->get_younger_sister();

    if (younger != NULL) younger->set_elder_sister(elder);

    if (elder == NULL)
	parent->set_oldest_daughter(younger);
    else
	elder->set_younger_sister(younger);
}


// add_node: Insert a new node as the eldest daughter of the root node.

local void add_node(sdyn* root, sdyn* new_node)
{
    new_node->set_parent(root);
    new_node->set_elder_sister(NULL);

    sdyn* od = root->get_oldest_daughter();
    new_node->set_younger_sister(od);

    if (od) od->set_elder_sister(new_node);
    root->set_oldest_daughter(new_node);
}


// m_sum: Calculate the total mass of the nodes on the list.

local real m_sum(sdyn* list[], int k_tuple)
{
    real sum = 0;
    for (int i = 0; i < k_tuple; i++) sum += list[i]->get_mass();

    return sum;
}


// pos_sum: Calculate the weighted m * pos sum for the nodes on the list.

local vector pos_sum(sdyn* list[], int k_tuple)
{
    vector sum = vector(0,0,0);

    for (int i = 0; i < k_tuple; i++)
	sum += list[i]->get_mass() * list[i]->get_pos();

    return sum;
}


// vel_sum: Calculate the weighted m * vel sum for the nodes on the list.

local vector vel_sum(sdyn* list[], int k_tuple)
{
    vector sum = vector(0,0,0);

    for (int i = 0; i < k_tuple; i++)
	sum += list[i]->get_mass() * list[i]->get_vel();

    return sum;
}


// tuple_size: Determine the instantaneous radius of the specified clump,
//	       relative to its center of mass.

local real tuple_size(sdyn* list[], int k_tuple)
{
    real max_dist = 0;
    vector cmpos = pos_sum(list, k_tuple) / m_sum(list, k_tuple);

    for (int i = 0; i < k_tuple; i++)
	max_dist = max(max_dist, abs(list[i]->get_pos() - cmpos));

    return max_dist;
}


// tuple_virial_radius: Return the virial radius of the specified clump.

local real tuple_virial_radius(sdyn* list[], int k_tuple)
{
    real total_mass = m_sum(list, k_tuple);
    real potential = 0;

    for (int i = 0; i < k_tuple; i++)
	for (int j = i+1; j < k_tuple; j++)
	    potential -= list[i]->get_mass() * list[j]->get_mass()
		          / abs(list[i]->get_pos() - list[j]->get_pos());

    return -0.5 * total_mass * total_mass / potential;
}


// distance_sq: Return the square of the minimum distance between any
//		element of the first list and any element of the second.

local real distance_sq(sdyn* list[], int k_tuple, sdyn* rest[], int n_rest)
{
    real min_dist_sq = VERY_LARGE_NUMBER;

    for (int i = 0; i < k_tuple; i++)
	for (int j = 0; j < n_rest; j++)
	    min_dist_sq = min(min_dist_sq,
			      square(list[i]->get_pos() - rest[j]->get_pos()));

    return min_dist_sq;
}


// is_bound: Return TRUE iff the specified clump is bound,

local bool is_bound(sdyn* list[], int k_tuple)
{
    vector cmvel = vel_sum(list, k_tuple) / m_sum(list, k_tuple);

    real kinetic = 0;
    for (int i = 0; i < k_tuple; i++)
	kinetic += list[i]->get_mass() * square(list[i]->get_vel());
    kinetic /= 2;

    real potential = 0;
    for (int i = 0; i < k_tuple; i++)
	for (int j = i+1; j < k_tuple; j++)
	    potential -= list[i]->get_mass() * list[j]->get_mass()
		          / abs(list[i]->get_pos() - list[j]->get_pos());

    if (kinetic + potential >= 0) return FALSE;
    return TRUE;
}


// make_tuple_cm: Combine the specified nodes, placing their center of mass
//		  at the start of the daughter list.

local sdyn* make_tuple_cm(sdyn* list[], int k_tuple, real radius,
			  bool meta = FALSE)
{
    if (k_tuple < 2) return NULL;

    sdyn* root = list[0]->get_parent();

    // Remove the list elements (and any daughters) from the tree.

    for (int i = 0; i < k_tuple; i++)
	detach_from_tree(list[i]);

    // Place the new CM at the start of the list (as the oldest daughter of
    // root) to prevent it from being picked up in the loop through the system.

    sdyn* cm = new sdyn;
    add_node(root, cm);

    // Reattach list elements;

    cm->set_oldest_daughter(list[0]);
    list[0]->set_elder_sister(NULL);
    for (int i = 0; i < k_tuple; i++) {
	list[i]->set_parent(cm);
	if (i > 0) list[i]->set_elder_sister(list[i-1]);
	if (i < k_tuple-1) list[i]->set_younger_sister(list[i+1]);
    }
    list[k_tuple-1]->set_younger_sister(NULL);

    // Determine CM mass, position, and velocity.

    real m_total = m_sum(list, k_tuple);
    vector cmpos = pos_sum(list, k_tuple) / m_total;
    vector cmvel = vel_sum(list, k_tuple)/ m_total;

    cm->set_mass(m_total);
    cm->set_pos(cmpos);
    cm->set_vel(cmvel);
    cm->set_time(root->get_time());

    cm->set_radius(radius);

    // Set the CM neighbor pointers to NULL for convenience in the
    // "find_and_make" loop.

    cm->set_nn_ptr(NULL);
    cm->set_nnn_ptr(NULL);

    // Offset the new components relative to the center of mass:

    for (int i = 0; i < k_tuple; i++) {
	list[i]->inc_pos(-cmpos);
	list[i]->inc_vel(-cmvel);
    }

    // New CM ID:

    bool index = TRUE;
    for (int i = 0; i < k_tuple; i++)
	if (list[i]->get_index() <= 0) index = FALSE;

    if (index) {
	int sum = 100 * list[0]->get_index();
	for (int i = 1; i < k_tuple; i++) sum += list[i]->get_index();
	cm->set_index(sum);				       // (For now...)
    }

    bool name = TRUE;
    for (int i = 0; i < k_tuple; i++)
	if (list[i]->get_name() == NULL) name = FALSE;

    if (name) {
	char temp[100];
	temp[0] = '\0';
	strcat(temp, (meta ? "{" : "("));
	for (int i = 0; i < k_tuple; i++) {
	    strcat(temp, list[i]->get_name());
	    if (i < k_tuple-1) strcat(temp, ",");
	}
	strcat(temp, (meta ? "}" : ")"));
	cm->set_name(temp);
    }

    return cm;
}


// identify_node: Print node address and ID.

local void identify_node(char* string, sdyn* bb)
{
    cerr << string << " = " << bb;
    if (bb != NULL) cerr << " = " << id(bb);
    cerr << endl;
}


// Specify the following as parameters?

#define MIN_PERI_FACTOR 2	// Unconditionally unstable within this peri
#define MAX_PERI_FACTOR 6	// Unconditionally stable outside this peri
#define ISOL_FACTOR 4		// For multiples only
#define ISOL_FACTOR_SQ		(ISOL_FACTOR * ISOL_FACTOR)

#define QUARANTINE_TIME_LIMIT	100	// Unit = outer orbit period
#define QSMA_TOL		0.01	// Acceptable variation in outer sma
#define QECC_TOL		0.1	// Acceptable variation in outer ecc


// set_tqflag: Recursively set the temporary quarantine flag of a node and
//	       its descendents.

local void set_tqflag(sdyn* b, int flag) {

    b->set_temp_quarantine_flag(flag);
    for_all_daughters(sdyn, b, daughter) set_tqflag(daughter, flag);	
}


// check_tqflag: Recursively clear all quarantine flags that aren't current.

local void check_tqflag(sdyn* b) {

    if (b->get_temp_quarantine_flag() == 0) b->set_quarantine_flag(0);
    for_all_daughters(sdyn, b, daughter) check_tqflag(daughter);
}


// set_quar: Recursively set the quarantine flag, time, and other parameters
//	     of the leaves lying under a node.

local void set_quar(sdyn* b, int flag, real time, real sma, real ecc) {

    if (b->get_oldest_daughter() == NULL) {
	b->set_quarantine_flag(flag);
	b->set_quarantine_time(time);
	b->set_quarantine_sma(sma);
	b->set_quarantine_ecc(ecc);
    }else
	for_all_daughters(sdyn, b, daughter)
	    set_quar(daughter, flag, time, sma, ecc);
}


// check_quar: Recursively check all quarantine flags, returning FALSE
//	       if any flag is not equal to the specified flag, or if any
//	       time, semi-major axis, or eccentricity is not equal to the
//	       specified value.  (This is a multiply redundant check,
//	       since all leaves in a quarantined system should contain the
//	       same information.)

local bool check_quar(sdyn* b, int flag, real time, real sma, real ecc) {

    // Only check leaves.  CM nodes are created at each invocation of
    // make_tree, so they do not contain any useful information.

    if (b->get_oldest_daughter() == NULL) {
	if (b->get_quarantine_flag() < 0
	    || b->get_quarantine_flag() != flag
	    || b->get_quarantine_time() != time
	    || b->get_quarantine_sma() != sma
	    || b->get_quarantine_ecc() != ecc) return FALSE;
    } else
	for_all_daughters(sdyn, b, daughter)
	    if (!check_quar(daughter, flag, time, sma,ecc)) return FALSE;

    return TRUE;
}


// check_quarantine: check quarantine status of the potentially stable
//		     multiple consisting of node b and neighbor n.  If not
//	             already in quarantine, begin new quarantine.  If already
//		     in quarantine, return TRUE if time limit is exceeded.
// 		     Extend and/or return FALSE in all other cases. In all
//		     cases, set the temp_quarantine_flag of all particles to 1.

local bool check_quarantine(sdyn* b, sdyn* n, real sma, real ecc, int debug)
{
    set_tqflag(b, 1);
    set_tqflag(n, 1);

    // Check that all particles under b and n have same flag and time.
    // If they do, check current time, return TRUE iff quarantine limit is up.
    // If they don't, start new quarantine.  A FALSE return means that the
    // system is not (yet) stable.

    // Note that b and/or n may be compound. In order to preserve quarantine
    // information between invocations of make_tree, copies of all data are 
    // stored in all nodes involved in the subsystem.

    sdyn* first = first_leaf(b);

    int  qflag = first->get_quarantine_flag();
    real qtime = first->get_quarantine_time();
    real qsma = first->get_quarantine_sma();
    real qecc = first->get_quarantine_ecc();

    if (check_quar(b, qflag, qtime, qsma, qecc)
	&& check_quar(n, qflag, qtime, qsma, qecc)) {

	// System is already in quarantine.

	if (abs(sma - qsma) < QSMA_TOL * qsma
	     && abs(ecc - qecc) < QECC_TOL * qecc) {

	    // Outer orbit parameters are acceptable.  Check the time since
	    // the system was placed in quarantine.

	    if (b->get_time() - qtime > QUARANTINE_TIME_LIMIT *
		TWO_PI * sqrt( (b->get_mass() + n->get_mass()) / pow(sma, 3)))
		return TRUE;
	    else
		return FALSE;
	}
    }

    // Start a new quarantine if this is a new system, or if the
    // parameters of an old system have drifted too much.

    qflag = (int)b / 4 + (int)n / 4;
    qtime = b->get_time();

    set_quar(b, qflag, qtime, sma, ecc);
    set_quar(n, qflag, qtime, sma, ecc);

    if (debug) cerr << ": new quarantine";
    return FALSE;
}


// find_and_make_binaries: Search for isolated, mutual near-neighbors.  Use the
//                         additional constraint that the pair be bound if the
//                         dynamics flag is set.  Also insist that the two
//                         components be isolated *at periastron* if the
//                         stability flag is set.
//
//			   Keep searching and creating binary nodes until the
//			   system has been entirely checked, or until some
//			   neighbor pointer has been corrupted (in which case
//			   we reset the pointers externally and repeat).

local int find_and_make_binaries(sdyn* b, bool dynamics, bool stability, 
				 int debug)
{
    int n_bin = 0;

    sdyn* bb = b->get_oldest_daughter();
    sdyn* last_bb = NULL;

    while (bb) {
	sdyn* nn = bb->get_nn_ptr();

	if (debug > 1) {
	    identify_node("bb", bb);
	    identify_node("nn", nn);
	}

	if (nn != NULL && nn->get_nn_ptr() == bb) {		// mutual

	    real radius = abs(bb->get_pos() - nn->get_pos());
	    bool is_binary = TRUE;
	    bool meta = FALSE;

	    if (debug) cerr << "mutual pair " << id(bb)
			    << " " << id(nn) << flush;

	    if (dynamics) {

		real m_total = bb->get_mass() + nn->get_mass();
		real eij = 0.5 * square(bb->get_vel() - nn->get_vel())
		           - m_total / radius;

		if (eij < 0) {					// bound

		    real sma = -0.5 * m_total / eij;	        // (> 0, note)
		    if (debug) cerr << ": bound, a = " << sma << flush;

		    if (stability) {

			// Determine periastron and check stability.

			vector rel_pos = bb->get_pos() - nn->get_pos();
			vector rel_vel = bb->get_vel() - nn->get_vel();

			real rdotv = rel_pos * rel_vel;
			real temp = 1 - radius / sma;

			real ecc = sqrt(max(0.0, rdotv * rdotv
					            / (m_total * sma)
					          + temp * temp));
			real periastron = sma * (1 - ecc);
			real size = max(bb->get_radius(), nn->get_radius());
	
			// Stability is based on periastron distance.
			// Implement quarantine here too.

			if (periastron < MIN_PERI_FACTOR * size) {

			    // Unconditionally unstable

			    is_binary = FALSE;
			    nn->set_nn_ptr(NULL);   // Avoid repeated effort

			    if (debug) cerr << ": unstable, e = " << ecc
				            << flush;

			} else if (periastron < MAX_PERI_FACTOR * size) {

			    // Set or check quarantine before deciding
			    // stability in this case.

			    if ( (is_binary = check_quarantine(bb, nn,
							       sma, ecc,
							       debug))
				   == FALSE) {
				nn->set_nn_ptr(NULL);
				if (debug) cerr << ": quarantined, e = " << ecc
						<< flush;
			    } else {
				if (debug) cerr << ": metastable, e = " << ecc
						<< flush;
				meta = TRUE;
			    }
			} else

			    // Unconditionally stable

			    if (debug) cerr << ": stable, e = " << ecc
				            << flush;
		    }

		    radius = sma;	// binary radius = semi-major axis

		} else {

		    is_binary = FALSE;
		    nn->set_nn_ptr(NULL);	    // Avoid repeated effort
		    if (debug) cerr << ": not bound, separation = " << radius
				    << flush;
		}
	    }

	    if (is_binary) {

		radius += bb->get_radius() + nn->get_radius();
		real r_crit_sq = ISOL_FACTOR_SQ * radius * radius;

		if (bb->get_nnn_dr2() > r_crit_sq
		    && nn->get_nnn_dr2() > r_crit_sq) {		// isolated

		    sdyn* list[2];
		    list[0] = bb;
		    list[1] = nn;
		    make_tuple_cm(list, 2, radius, meta);

		    n_bin++;
		    bb = last_bb;

		    // Don't consider neighbors of bb or nn this time around

		    for_all_daughters(sdyn, b, bbb)
			if (bbb->get_nn_ptr() == bb
			    || bbb->get_nn_ptr() == nn
			    || bbb->get_nnn_ptr() == bb
			    || bbb->get_nnn_ptr() == nn)
			    bbb->set_nn_ptr(NULL);
		} else {
		    nn->set_nn_ptr(NULL);
		    if (debug) cerr << ": not isolated" << flush;
		}
	    }
	    if (debug) cerr << endl;
	}

	last_bb = bb;

	// Find the next node to condsider:

	if (bb != NULL)

	    bb = bb->get_younger_sister();

	else {

	    // We have just merged the first non-center-of-mass element
	    // on the list.  Find the next node by looking for a particle
	    // whose nn pointer is non-null (CM pointers are NULL by
	    // construction).  Note that some non-CM particles will also
	    // have NULL nn pointers, if they lay too close to a binary,
	    // so also terminate if the last daughter is reached.

	    bb = b->get_oldest_daughter();
	    while (bb != NULL && bb->get_nn_ptr() == NULL)
		bb = bb->get_younger_sister();
	}
    }

    return n_bin;
}


// find_and_make_tuple: Search for an isolated k-tuple.  Use the additional
//			constraint that a tuple be bound if the dynamics
//			flag is set.  Return when the first tuple is found.
//                      Treat binaries (k_tuple = 2) as a special case.

local int find_and_make_tuple(sdyn* b, int k_tuple, bool dynamics,
			      bool stability, int debug)
{
    if (debug > 1) cerr << "find_and_make_tuple: k = " << k_tuple << endl;

    if (k_tuple < 2)

	return 0;

    else if (k_tuple == 2)

	return find_and_make_binaries(b, dynamics, stability, debug);

    else {

	sdyn** list = new sdyn*[k_tuple+1];

	for_all_daughters(sdyn, b, bi) {
	    list[0] = bi;
	    int nlist = 1;
	    for (int i = 1; i <= k_tuple; i++) {
		sdyn* nn = list[i-1]->get_nn_ptr();
		bool dup = FALSE;
		for (int j = 0; j < i; j++)		// nn already on list?
		    if (list[j] == nn) dup = TRUE;
		if (!dup)
		    list[nlist++] = nn;
		else
		    break;
	    }

	    if (debug > 1) cerr << "bi = " << id(bi)
		                << "  nlist = " << nlist << endl;

	    // All near-neighbors but the were put on the list.

	    if (nlist == k_tuple) { 

		// Now list(0,...,k_tuple-1) is a candidate k-tuple set,
		// since the nearest neighbor of each element is also a
		// member of the set.

		if (debug) {
		    cerr << "candidate " << k_tuple << "-tuple";
		    for (int i = 0; i < k_tuple; i++)
			cerr << " " << id(list[i]);
		    cerr << flush;
		}

		// Make a list of the other particles in the system.

		int n = count_daughters(b);
		sdyn** rest = new sdyn*[n-k_tuple+1];
		int n_rest = 0;

		for_all_daughters(sdyn, b, bj) {
		    bool dup = FALSE;
		    for (int i = 0; i < k_tuple; i++)
			if (bj == list[i]) {
			    dup = TRUE;
			    break;
			}
		    if (!dup) rest[n_rest++] = bj;
		}

		// Optionally check bound-ness:

		bool is_binary = TRUE;
		if (dynamics && !is_bound(list, k_tuple)) is_binary = FALSE;

		if (is_binary) {

		    // Check isolation (NOTE definition of tuple radius):

		    real radius = max(tuple_size(list, k_tuple),
				      tuple_virial_radius(list, k_tuple));

		    for (int i = 0; i < k_tuple; i++)
			radius += list[i]->get_radius(); // [or max(radius)?]

		    real r_crit_sq = ISOL_FACTOR_SQ * radius * radius;

		    real dist_sq = distance_sq(list, k_tuple, rest, n_rest);
		    delete rest;

		    if (r_crit_sq < dist_sq) {
			sdyn* cm = make_tuple_cm(list, k_tuple, radius);
			if (debug) {
			    cerr << ", radius = " << radius << endl;
			    pp(cm, cerr);
			}
			return 1;
		    }
		    if (debug) cerr << ": not isolated\n";

		} else {
		    delete rest;
		    if (debug) cerr << ": not bound\n";
		}
	    }
	}
	delete list;
    }

    return 0;
}


// make_tree: Recursively seek k-tuples, for k >= 2.  Continue search at
//	      the k = 2 level if any structure is found.

void make_tree(sdyn* b, bool dynamics, bool stability, int k_max, int debug)
{
    b->flatten_node();  // (Just in case...)

    // Look for tuples and create a tree.

    if (k_max <= 0) k_max = 1000000;
    set_tqflag(b, 0);

    int nt;
    do {
	nt = 0;
	for (int k_tuple = 2;
	     k_tuple < min(k_max + 1, count_daughters(b));
	     k_tuple++) {

	    // If k_tuple = 2, initialize the neighbor pointers.
	    // If k_tuple = 3, reconstruct the pointers, in case some
	    //		       were corrupted during the binary search

	    if (k_tuple < 4) find_neighbors(b);

	    if ( (nt = find_and_make_tuple(b, k_tuple, dynamics,
					   stability, debug)) > 0)
		break;
	}
    } while (nt > 0);

    check_tqflag(b);
}

#else

local void delete_node(sdyn* b)
{

    // Recursively delete node b and its descendents.

    sdyn* bi = b->get_oldest_daughter();
    while (bi) {
	sdyn* tmp = bi->get_younger_sister();
	delete_node(bi);
	bi = tmp;
    }
    delete b;
}

#define USAGE "    Usage:  make_tree [-b] [-d] [-D [#]] [-k #] [-o] [-s]"

main(int argc, char **argv)
{
    sdyn* root;             // pointer to the N-body system

    int  debug = 0;
    bool dynamics = FALSE;
    int k_max = 10000000;
    bool output = TRUE;
    bool stability = FALSE;

    check_help();

    int i = 0;
    while (++i < argc) if (argv[i][0] == '-')
	switch (argv[i][1]) {
	    case 'b': k_max = 2;
		      break;
	    case 'd': dynamics = 1 - dynamics;
		      break;
	    case 'D': i++;
		      if (i == argc || argv[i][0] == '-') {
			  debug = 1;
			  i--;
		      } else
			  debug = atoi(argv[i]);
		      break;
	    case 'k': k_max = atoi(argv[++i]);
		      break;
	    case 'o': output = 1 - output;
		      break;
	    case 's': stability = 1 - stability;
		      break;
            default:  cerr << USAGE << endl;
		      get_help();
	}

    // Note that we are simply assuming that the tree is flat initially...

    while (root = get_sdyn(cin)) {

//        sdyn* copy_root;	      // pointer to a copy of the N-body system
//        copy_tree(root, copy_root); // Should really save a copy
                                      // before proceeding...

	if (root->get_name() == NULL) root->set_name("root");

	make_tree(root, dynamics, stability, k_max, debug);

	if (debug) pp(root, cerr);
	if (output) put_sdyn(cout, *root);

	// Delete the newly constructed tree.

	delete_node(root);
    }
}
#endif
