
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Initiate and maintain unperturbed motion.
//
// Externally visible functions:
//
//	void hdyn::print_unperturbed_binary_params
//	void hdyn::update_dyn_from_kepler
//	bool hdyn::is_close_pair
//	bool hdyn::is_unperturbed_and_approaching
//	void hdyn::startup_unperturbed_motion
//	bool hdyn::integrate_unperturbed_motion
//	real hdyn::set_unperturbed_timestep
//	real hdyn::get_unperturbed_steps
//	bool hdyn::is_weakly_perturbed
//	bool hdyn::is_stable
//
// To do:  Mixing of many types of motion is too complex.  Especially
//	   problematic is the identification of pericenter reflection by
//	   checking for recession, which can fail in a variety of ways
//	   in marginal cases.
//
//	   Also, picking up unperturbed motion near pericenter biases the
//	   unperturbed timestep on an extremely short step, and can lead
//	   to scheduling problems.

#include "hdyn.h"
#include <star/dstar_to_kira.h>

#define USE_DT_PERTURBERS		true

//-------------------------------------------------------------------------
//
// Parameters now contained in hdyn::kira_options:
//
//#define OPTIMIZE_SCHEDULING		true
//
// Define what level of unperturbed motion is permitted by default:
//
//#define ALLOW_UNPERTURBED		true
//
//#define ALLOW_MULTIPLES		false	// "false" here turns off all
//						// unperturbed multiple motion
//
//-------------------------------------------------------------------------

// The following parameters define "unperturbed" motion.
//
// Note that gamma2 is used in determining membership in the perturber
// list.  Note also that d_min_sq does not appear in this file, so the
// criteria for unperturbed motion are independent of the criteria used
// in constructing and maintaining the tree (except that unperturbed
// motion can only occur in low-level nodes, so the binary must first
// be picked up by hdyn_tree).
//
// Where relevant, numbers are tailored to the default gamma = 1.e-7.
//
//-------------------------------------------------------------------------
//
// OUTSIDE_SEMI:  Comparison of orbital separation with semi-major axis for
// full merging.  Criterion is
//
//	SEP_FACTOR * separation > semi-major axis
//
// Macro is used three times, in is_unperturbed_and_approaching(),
// is_stable(), and set_unperturbed_timestep().
//
// SEP_FACTOR is slightly greater than 1 in order to accommodate nearly
// circular binaries

#define SEP_EPSILON			0.01	// arbitrary, but small
#define SEP_FACTOR			(1.0 + SEP_EPSILON)

#define OUTSIDE_SEMI(sep, sma)		((SEP_FACTOR)*(sep) > (sma))

// The comparison is often performed on elements of a kepler structure (k):

#define KEP_OUTSIDE_SEMI(k) \
    ((SEP_FACTOR)*((k).get_separation()) > ((k).get_semi_major_axis()))

//-------------------------------------------------------------------------
//
// The remaining parameters were originally #defined here, but are
// now part of the hdyn::kira_options class (Steve, 6/00).  Defaults
// are defined in kira_defaults.h.
//
//
// MIN_UNPERT_STEPS:  Minimum number of steps permitted in unperturbed motion.
//
//#define MIN_UNPERT_STEPS		5
//
// FULL_MERGE_TOLERANCE:  Limit on perturbation squared (independent of
// gamma2, as of 1/14/02) at which we permit merging of an entire binary
// orbit.  Used once, in function is_unperturbed_and_approaching().
//
//#define FULL_MERGE_TOLERANCE		1e-8
//#define RELAX_FACTOR			10	// relaxed factor to continue
//
// Allowing unperturbed motion at a greater value of gamma than is
// used in determining neighbors is reasonable since dominant
// perturbative terms are expected to be periodic, and so are
// ineffective over a full orbit period.
//
// Relaxed criterion to continue is probably reasonable (for true
// binaries only), as the first unperturbed step takes them to
// apocenter, where the external perturbation is greater.
//
// PARTIAL_MERGE_FACTOR:  Factor below gamma2 at which we permit
// merging of at least part of an orbit.  Used once, in function
// is_unperturbed_and_approaching().
//
//#define PARTIAL_MERGE_FACTOR		1e-2
//
// FULL_MERGE_TOL_FOR_CLOSE_BINARY:  Absolute limit on perturbation_squared
// (INDEPENDENT of gamma2, note) at which an unperturbed binary is regarded
// as "close."  Used once, in function is_close_pair().
//
//#define FULL_MERGE_TOL_FOR_CLOSE_BINARY	1e-4
//
// The following parameters are used in is_stable() for assessing
// multiple stability.
//
// MULTIPLE_MERGE_TOLERANCE:  Perturnation limit (independnt of gamma2,
// as of 1/14/01) at which we permit merging of an inner binary orbit.
// Defines "weakly perturbed outer orbit."
//
//#define MULTIPLE_MERGE_TOLERANCE		1e-8	// was 1.e-10
//
// UNCONDITIONAL_STABLE_FAC:  Simplest approach is to regard a multiple as
// stable if
//
//	outer periastron / inner semi-major axis > UNCONDITIONAL_STABLE_FAC
//
// for (both) inner binaries.
//
//#define UNCONDITIONAL_STABLE_FAC		5	// was 9
//
// Aarseth and Mardling's criterion is a little more elaborate than a simple
// distance ratio:
//
//	outer periastron / inner semi-major axis
//	    > AARSETH_STABLE_FAC * pow((1+q)*(1+e_outer)/sqrt(1-e_outer), 0.4)
//
// where e_outer is the eccentricity of the outer orbit and q is the mass
// ratio m_outer / m_binary.  Note that 1+q = m_total/m_binary.
//
// An empirical inclination factor should also be included in either case.
//
//#define USE_AARSETH_CRITERION		false
//#define AARSETH_STABLE_FAC		2.8
//
// PARTIAL_STABLE_FAC:  Allow inner binary to be unperturbed for part of
// the outer orbit (questionable!) if
//
//	outer periastron / inner semi-major axis > PARTIAL_STABLE_FAC
//
//#define PARTIAL_STABLE_FAC		30	// large to reduce tidal effect
//
// This should really also contain a mass factor...

//-------------------------------------------------------------------------

// STATIC flags to keep track of unperturbed motion:

static int init_binary_type;	// binary_type set by is_unperturbed_and...
static int binary_type;		// eventual binary_type

#define UNKNOWN_STATUS			0
#define STABLE_OUTER 			1
#define STABLE_INNER			2
#define NOT_APPROACHING			3
#define PERICENTER_REFLECTION		4
#define FULL_MERGER			5
#define CONTINUE_MERGER			6
#define PERTURBED			7
#define MULTIPLE_CM			8
#define UNPERTURBED_MULTIPLE_COMPONENT	9
#define UNKNOWN_PERTURBERS		10
#define DEFERRED			11

static char* bt[12] = {"unknown status",
		      "stable multiple outer binary",
		      "stable multiple inner binary",
		      "components not approaching",
		      "pericenter reflection",
		      "full binary merger",
		      "continue existing merger",
		      "perturbed binary",
		      "multiple system",
		      "unperturbed multiple component",
		      "no valid perturber list",
		      "unperturbed motion deferred"};

static int multiple_type;

#define NOT_MULTIPLE			0
#define FULL_MULTIPLE			1
#define APOCENTER_REFLECTION		2

static char* tt[3] = {"not multiple",
		      "full_multiple merger",
		      "apocenter reflection"};

//------------------------------------------------------------------------

// Local functions:

local inline real get_energy(hdyn* b, real& separation)
{
    // Compute the energy of b relative to its binary sister.

    if (b->is_top_level_node()) return 0;

    hdyn* s = b->get_binary_sister();
    separation = abs(b->get_pos() - s->get_pos());

    return 0.5*square(b->get_vel() - s->get_vel())
	      - b->get_parent()->get_mass() / separation;
}

local inline real get_semi(hdyn* bcm)
{
    // Compute the semi-major axis of the binary daughters of bcm.

    real semi = 0;

    if (bcm->is_parent()) {

	// Get semi-major axis of the orbit.
	// Use existing kep if orbit is already unperturbed.

	if (bcm->get_oldest_daughter()->get_kepler()) {
	    semi = bcm->get_oldest_daughter()->get_kepler()
					     ->get_semi_major_axis();
	} else {
	    real sep;
	    real energy = get_energy(bcm->get_oldest_daughter(), sep);
	    semi = -0.5 * bcm->get_mass() / energy;
	}
    }

    if (semi < 0) semi = -VERY_LARGE_NUMBER;

    return semi;
}

local inline real get_period(real mass, real sma)
{
    // Compute the period of a binary with specified mass and
    // semi-major axis.

    return TWO_PI * sqrt(sma*sma*sma / mass);
}

local inline void print_found_multiple(hdyn* b,
				       bool real_multiple,
				       kepler& outer,
				       real inner_semi1,
				       real inner_semi2,
				       real mass1,
				       real mass2)
{
    int p = cerr.precision(HIGH_PRECISION);
    cerr << endl;

    if (!real_multiple)
	cerr << endl << "impending multiple: ";

    cerr << b->get_parent()->format_label();

    if (real_multiple)
	cerr << " is a hierarchical stable multiple";
    else
	cerr << " is a multiple CM";

    cerr << endl << "    at system time " << b->get_system_time() << endl;
    cerr << "    components " << b->format_label();
    cerr << " and " << b->get_binary_sister()->format_label();

    if (!real_multiple)
	cerr << "; outer orbit is weakly perturbed:" << endl;

    cerr.precision(8);		// nonstandard precision

    cerr << "    perturbation = " << sqrt(b->get_perturbation_squared())
	 << endl
	 << "    outer peri   = " << outer.get_periastron()
	 << "    inner semi   = ";

    if (inner_semi1 > 0) cerr << inner_semi1 << "  ";
    if (inner_semi2 > 0) cerr << inner_semi2;

    cerr << endl
	 << "    outer period = " << outer.get_period()
	 << "    inner period = ";

    if (inner_semi1 > 0) cerr << get_period(mass1, inner_semi1) << "  ";
    if (inner_semi2 > 0) cerr << get_period(mass2, inner_semi2);

    cerr << endl
	 << "    parent dt = " << b->get_parent()->get_next_time()
	     				- b->get_system_time()
	 << endl;

    cerr.precision(p);
}

local inline void print_startup_message(hdyn * b, int type, bool new_kep)
{
    int p = cerr.precision(HIGH_PRECISION);
    cerr << endl;

    if (new_kep)
	cerr << "starting new";
    else
	cerr << "continuing existing";

    cerr << " unperturbed motion for "
	 << b->format_label();
    cerr << " and " << b->get_binary_sister()->format_label();

    if (b->get_slow())
	cerr << " (slowdown = " << b->get_kappa() << ")";   // shouldn't occur
//    else
//	cerr << " (no slowdown)";

    cerr << endl;

    cerr << "    at time " << b->get_time()
	 << "  " << " (system_time = " << b->get_system_time() << ")"
	 << endl;

    cerr << "    " << bt[type] << " (binary_type = " << type << ")";

    cerr.precision(8);	// nonstandard precision

    if (b->get_kepler())
	cerr << "  period = " << b->get_kepler()->get_period();

    cerr << endl;
    cerr.precision(p);

    cerr << "    perturbation = " << sqrt(b->get_perturbation_squared());

    if (b->get_parent()->get_nn() && b->get_parent()->get_nn()->is_valid())
	cerr << "  (" << b->get_parent()->get_nn()->format_label() << ")";

    cerr << endl;

#if 0
     if (streq(b->format_label(), "xxx"))
	 b->get_kepler()->print_all();
#endif

}

local inline void print_binary_data(hdyn* bi)
{
    PRL(bi->get_perturbation_squared());
    PRL(bi->get_unperturbed_timestep());
    PRL(bi->get_fully_unperturbed());
    PRL(bi->get_d_nn_sq());
    PRL(bi->get_perturbation_radius_factor());
    PRL(bi->find_perturber_node());
    if (bi->find_perturber_node())
	PRL(bi->find_perturber_node()->get_valid_perturbers());
    PRL(bi->get_d_coll_sq());
    PRL(bi->get_coll()->format_label());
}

local inline bool is_multiple(hdyn* b)   // ("multiple" means 3 or more stars)
{
    // Return true if low-level node b is one of the outer components
    // of a multiple system.

    if (b->is_top_level_node()) return false;
    return !( b->is_leaf() && b->get_binary_sister()->is_leaf() );
}

local inline real dt_overshoot(hdyn *b)
{
    // Return the amount by which a binary's unperturbed step may
    // overshoot its parent's next step.

    if (!b->get_parent())
	return 0;
    else
	return 0.9 * b->get_parent()->get_timestep();
}

//----------------------------------------------------------------------

#include "hdyn_inline.C"	// used only by check_perturbers()
				// and associated functions

local inline real dt_perturbers(hdyn *b)
{
    // Return the time for any perturber of multiple CM b to exceed
    // the threshhold for unperturbed multiple motion, or -1 if no
    // perturber list exists.

    // Effectively replace gamma2 by multiple_merge_tolerance in
    // the usual perturbation criterion.

    hdyn *top = b->get_top_level_node();
    real t_min = VERY_LARGE_NUMBER;

    if (top->get_valid_perturbers()) {

	real scale = binary_scale(top);
	real gamma = sqrt(b->get_kira_options()->multiple_merge_tolerance);

	// Loop over perturbers.

	int np = top->get_n_perturbers();
	hdyn **plist = top->get_perturber_list();

	for (int j = 0; j < np; j++) {

	    hdyn *p = plist[j];
	    if (!p->is_top_level_node()) {	      // inefficient -- covers
		if (p->get_elder_sister()) continue;  // some components twice
		p = p->get_top_level_node();
	    }

	    // Estimate the time needed for node p to perturb top
	    // at level gamma.  The critical distance for p to be
	    // such a perturber is rpert.

	    real rpert3 = crit_separation_cubed(top, p->get_mass(),
						scale, gamma);

	    real rpert = pow(rpert3, 1.0/3);

	    // Time estimate is a combination of the free-fall and
	    // crossing times.

	    vec dpos = top->get_pos() - p->get_pos();
	    vec dvel = top->get_vel() - p->get_vel();
	    vec dacc = top->get_acc() - p->get_acc();

	    real dr = abs(dpos);
	    real vr = dvel * dpos / dr;
	    real ar = dacc * dpos / dr;

	    real tr = time_to_radius(dr - rpert, vr, ar);

	    if (tr < t_min) t_min = tr;
	    if (tr <= 0) break;
	}

    } else

	t_min = -1;

    return t_min;
}

local void check_perturbers(hdyn *b)
{
    // Print the time for any perturber of multiple CM b to exceed
    // the threshhold for unperturbed multiple motion.

    real t_min = dt_perturbers(b);

    if (t_min >= 0)

	cerr << "    " << b->get_top_level_node()->get_n_perturbers()
	     << " perturbers, perturber crossing time = "
	     << t_min << endl;

    else

	cerr << "    no valid perturber list" << endl;
}

//----------------------------------------------------------------------



// Modifiers for unperturbed options:

void set_allow_unperturbed(hdyn *b,
			   bool value)	// default = true
{
    b->get_kira_options()->allow_unperturbed = value;
    if (!b->get_kira_options()->allow_unperturbed)
	b->get_kira_options()->allow_multiples = false;
}

void set_allow_multiples(hdyn *b,
			 bool value)	// defualt =true
{
    b->get_kira_options()->allow_multiples = value;
    if (b->get_kira_options()->allow_multiples)
	b->get_kira_options()->allow_unperturbed = true;
}

void toggle_unperturbed(hdyn *b, int level)
{
    if (level == 0)
	b->get_kira_options()->allow_unperturbed
	    = !b->get_kira_options()->allow_unperturbed;
    else
	b->get_kira_options()->allow_multiples
	    = !b->get_kira_options()->allow_multiples;

    // Consistency:

    if (!b->get_kira_options()->allow_unperturbed)
	b->get_kira_options()->allow_multiples = false;
    if (b->get_kira_options()->allow_multiples)
	b->get_kira_options()->allow_unperturbed = true;
}

void print_unperturbed_options(hdyn *b)
{
    cerr << "unperturbed options:  allow_unperturbed = "
	 << b->get_kira_options()->allow_unperturbed
	 << ",  allow_multiples = "
	 << b->get_kira_options()->allow_multiples
	 << endl;
}



void hdyn::print_unperturbed_binary_params()
{
    cerr << "    perturbed timestep for " << format_label()
	 << " is " << timestep << endl;

    cerr << "    pos*vel = " << pos*vel
	 << endl
	 << "    a = " << kep->get_semi_major_axis()
	 << "  e = " << kep->get_eccentricity()
	 << "  r = " << kep->get_separation()
	 << endl
	 << "    peri = " << kep->get_periastron()
	 << "  apo = " << kep->get_semi_major_axis()
				* (1 + kep->get_eccentricity())
	 << endl
	 << "    period = " << kep->get_period()
	 << "  perturbation = " << sqrt(perturbation_squared)
	 << endl;

    PRI(4); PRL(unperturbed_timestep);

    PRI(4); PRL(get_parent()->get_next_time());
    PRI(4); PRL(get_parent()->timestep);
    PRI(4); print_nn(get_parent(), 2);
}

void hdyn::update_dyn_from_kepler(bool need_acc_and_jerk)	// default = true
{
    if (diag->unpert_function_id) {
	cerr << ">> update_dyn_from_kepler for "
	     << format_label() << endl;
    }

    hdyn *sister = get_binary_sister();

    // Note: if this function is called from integrate_unperturbed_motion,
    // time and system_time are already the same.  However, if it is called
    // from dissociate_binary, the binary is currently at a retarded time
    // and we must predict the kepler to system_time before terminating it.

    if (diag->report_end_unperturbed) {
	if (time != system_time) {
	    cerr << "in update_dyn_from_kepler: ";
	    PRC(time), PRL(system_time);
	}
    }

    time = system_time;

    // Update kepler to new time.  Special treatment is needed in the case
    // of slow binary motion, as the the kepler structure really describes
    // the orbit in tau, not time.  (However, the "time" attached to the
    // kepler is the actual time at which the unperturbed motion started.)
    // Note that this only works as coded because we are doing pericenter
    // reflection, and the kepler structure is about to be deleted.

    if (!slow)
	kep->transform_to_time(time);
    else
	kep->transform_to_time(time - unperturbed_timestep
			       		* (1 - 1.0/get_kappa()));

    // Note that the latter expression implies the former, as get_kappa() = 1
    // for binaries with slow = NULL, but it seems cleaner this way.

    // Compute dyn components from updated kepler values.

    real factor = -sister->mass / (mass + sister->mass);
    pos = kep->get_rel_pos() * factor;
    vel = kep->get_rel_vel() * factor;
    pred_pos = pos;
    pred_vel = vel;

    prev_posvel = posvel;
    posvel = pos*vel;		    // used in is_unperturbed_and_approaching

    sister->time = time;

    real sfactor = 1 + factor;
    sister->pos = kep->get_rel_pos() * sfactor;
    sister->vel = kep->get_rel_vel() * sfactor;
    sister->pred_pos = sister->pos;
    sister->pred_vel = sister->vel;

    // Convenient, but possibly unnecessary, to update accelerations,
    // jerks, etc. for the binary components at the end of every
    // unperturbed step.  Doesn't seem to be expensive to do this.
    // However, should it be desirable to limit computation of acc and
    // jerk, use need_acc_and_jerk to do so.  We don't need acc and
    // jerk at the end of a normal (continuing) unperturbed step, but
    // we do need them if unperturbed motion is about to end.  Modify
    // integrate_unperturbed_motion as needed.
    //							(Steve, 5/02)

    // Continue to set perturbation_squared for use elsewhere.

    if (!need_acc_and_jerk) {

	hdyn* pnode = find_perturber_node();
	PRL(pnode);

	if (pnode && pnode->valid_perturbers && pnode->n_perturbers == 0)
	    perturbation_squared = 0;
	else
	    need_acc_and_jerk = true;
    }

    if (need_acc_and_jerk) {

	// Recompute perturbation_squared based on updated pos and vel.

	clear_interaction();

	// Doing these steps here allows us to bypass calculate_acc_and_jerk
	// and go directly to calculate_acc_and_jerk_on_low_level_node.

	d_coll_sq = VERY_LARGE_NUMBER;
	coll = NULL;
	sister->d_coll_sq = VERY_LARGE_NUMBER;
	calculate_acc_and_jerk_on_low_level_node();

	store_old_force();
    }

    update_binary_sister(this);		// seems to repeat some of the above...
}

bool hdyn::is_close_pair()
{
    if (kep == NULL) return false;
    if (kep->get_energy() >= 0) return false;

    real rp = kep->get_periastron();
    real sum_radius = radius + get_binary_sister()->radius;

    if (rp < 5 * sum_radius) {

	// This criterion is *very bad*: it should reflect the
	// Roche radius...

	// When this test is done, the system is already close
	// to apocenter, so perturbation is reasonable.

	if (perturbation_squared
	    	< options->full_merge_tol_for_close_binary) {

	    // cerr << "time = " << system_time << " ";

	    // pretty_print_node(cerr);

	    // cerr << " close pair criterion ";
	    // cerr << perturbation_squared << " " ;
	    // PRC(rp); PRL(sum_radius);

	    return true;
	}
    }
    return false;
}



// is_weakly_perturbed:  Return true iff 'this' system is lightly perturbed
//			 and satisfies some other basic acceptance criteria
//			 for the outer binary of an unperturbed multiple
//			 system.
//
//			 New function written by Steve, 8/98.

static char* wp[11] = {"unknown status",
		       "not low-level node",
		       "unbound orbit",
		       "no perturber node",
		       "no valid perturbers",
		       "perturbation too large",
		       "extended outer orbit",
		       "inside semi and perturbed",
		       "short parent time step",
		       "already uperturbed",
		       "weakly perturbed"};

bool hdyn::is_weakly_perturbed(int& status)
{
    if (diag->unpert_function_id) {
	cerr << ">> check is_weakly_perturbed for "
	     << format_label() << " at time " << system_time << endl;
    }

    status = 0;

    if (!is_low_level_node()) {
	status = 1;
	return false;
    }

    // 'This' is a component of a binary, possibly a multiple.
    // Return true if the binary defined by this and its binary
    // sister is weakly perturbed (i.e. not strongly perturbed).

    real separation, energy;

    if ((energy = get_energy(this, separation)) >= 0) {
	status = 2;
	return false;
    }

    hdyn* pnode = find_perturber_node();

#if 0

    // Old code: Don't proceed until perturber-list structure
    //		 is properly set.

    if (!pnode) {
	status = 3;
	return false;
    }

#else

    // New code mimics that for binaries in function
    // is_unperturbed_and_approaching().
    //					    (Steve, 4/99)

    if (!pnode) {

	// Possible that a newly formed center of mass node hasn't
	// yet taken a step, so no perturber list exists.  However,
	// we shouldn't break up an unperturbed system.

	if (kep) {

	    // Assume that some reorganization has just taken place in
	    // the unperturbed multiple system of which 'this' is a
	    // member, and allow the unperturbed motion to continue.

	    status = 9;
	    return true;

	} else {

	    status = 3;
	    return false;

	}
    }

#endif

    if (!pnode->valid_perturbers) {

	// Really defer multiple motion in this case.

	status = 4;
	return false;
    }

    // Definition of "not weakly perturbed" (see *** note below):

    if (perturbation_squared > options->multiple_merge_tolerance) {
	status = 5;
	return false;
    }

    // Avoid merging if the outer orbit is too extended.
    // Should perhaps make this condition consistent with the
    // relax_factor used elsewhere?

    real semi_major_axis = -0.5 * parent->get_mass() / energy;

    if (semi_major_axis > 5 * separation) {	  // 5 is ~arbitrary...
	status = 6;
	return false;
    }

    // Require the outer orbit to be outside its semi-major axis,
    // or completely unperturbed.

    // Note from Steve (9/00): the OUTSIDE_SEMI requirement here seems
    // unnecessary -- better to retain more compact configurations?
    // If it is changed, must modify logic in set_unperturbed_timestep
    // too.  *** However, use of perturbation_squared above is acceptable
    // only if system is outside peri, so maybe OK to keep this as is. ***
    // If we change the OUTSIDE_SEMI condition, then should use the
    // perturbation normalized to the sma...

    bool is_weakly_pert = false;

    if (OUTSIDE_SEMI(separation, semi_major_axis))
	is_weakly_pert = true;
    else if (pnode->n_perturbers == 0)
	is_weakly_pert = true;

    status = 7;

    if (is_weakly_pert) {

	// The following condition isn't directly related to the strength of
	// the perturbation, but may be checked if the binary turns out to
	// be stable.  May as well check these here and avoid the stability
	// check if it is not needed.

	// The minimum unperturbed step for a multiple is currently the
	// outer orbit period.  Return here if that will not be possible,
	// given the parent time step.

	real period = get_period(parent->get_mass(), semi_major_axis);
	real pdt2 = (real)(get_parent()->get_next_time() - time)
	    			+ dt_overshoot(this);

	if (period <= pdt2 || USE_DT_PERTURBERS) {	// latter condition
							// ==> consider later
	    status = 10;
	    return true;

	} else if (diag->report_multiple && diag->unpert_report_level > 0) {

	    cerr << endl << "Multiple " << get_parent()->format_label()
		 << " is weakly perturbed, but parent time step is too short"
		 << endl
		 << "    period = " << period << "  pdt2 = " << pdt2
		 << "  parent step = " << get_parent()->timestep
		 << endl;

	    check_perturbers(this);

#if 0
	    pp3(get_parent());

	    print_binary_from_dyn_pair(this, get_younger_sister());
	    cerr << endl;
	    if (is_parent()) {
		print_binary_from_dyn_pair(get_oldest_daughter(),
					   get_oldest_daughter()
					     ->get_younger_sister());
		cerr << endl;
	    }
	    if (get_younger_sister()->is_parent()) {
		print_binary_from_dyn_pair(get_younger_sister()
					     ->get_oldest_daughter(),
					   get_younger_sister()
					     ->get_oldest_daughter()
					     ->get_younger_sister());
		cerr << endl;
	    }
#endif

	}

	status = 8;

    }

    return false;
}



#define NEAR_MULTIPLE_FAC	1.2

// is_stable:  Return true iff this is a stable system.
//
//	       New (recursive) function written by Steve, 8/98.
//
// Note: Apocenter reflection is somewhat suspect... (Steve, 4/99)

bool hdyn::is_stable(int& status,
		     bool top_level)	// default = true
{
    if (diag->unpert_function_id) {
	cerr << ">> check is_stable for "
	     << format_label() << endl;
    }

    status = 0;

    if (!is_low_level_node()) {
	status = 1;
	return false;
    }

    // 'This' is a component of a binary, possibly a multiple.
    // Return true iff the binary defined by this and its binary
    // sister is stable.
    //
    // Unperturbed criteria for this to be a stable system:
    //
    // (a) the interior motion is stable, meaning that the ratio of
    //     the semimajor axis of each inner binary to the pericenter
    //	   of the outer binary is less than some critical value.
    //
    // (b) each component is stable: a single star, a binary, or a
    //	   stable multiple.

    hdyn* sister = get_binary_sister();
    if (!sister) {
	status = 2;
	return false;
    }

    // Binaries are stable:

    if (is_leaf() && sister->is_leaf()) {
	status = 3;
	return true;
    }

    // Set up a kepler structure describing the outer orbit.
    // Use existing kep if the outer orbit is already unperturbed.

    kepler outerkep;

    if (kep == NULL) {
	outerkep.set_time(time);
	outerkep.set_total_mass(parent->get_mass());
	outerkep.set_rel_pos(pos - sister->pos);
	outerkep.set_rel_vel(vel - sister->vel);
	outerkep.initialize_from_pos_and_vel();
    } else
	outerkep = *kep;

    real outer_peri = outerkep.get_periastron();

    real inner_semi1 = get_semi(this);
    real inner_semi2 = get_semi(sister);
    real inner_semi_sum = inner_semi1 + inner_semi2;

    if (inner_semi_sum < 0) {
	status = 4;
	return false;			// unbound inner orbit
    } else if (inner_semi_sum == 0) {
	status = 5;
	return true;			// binary is stable -- shouldn't happen
    }					// (already tested above)

    // Check criterion (a), applying the stability criterion to both
    // inner binaries, if appropriate.

    real peri_fac = 0;

    if (options->use_aarseth_criterion) {

	real e_outer = outerkep.get_eccentricity();
	real total_mass = outerkep.get_total_mass();

	// Apply Aarseth's criterion in the form
	//
	//    (outer_peri/AARSETH_STABLE_FAC)^5 * (1-e_outer)
	//					/ (total_mass*(1+e_outer))^2
	//
	//		> semi_major_axis^5 / binary_mass^2

	real binary_fac = 0;

	if (inner_semi1 > 0)
	    binary_fac = pow(inner_semi1, 5) / pow(mass, 2);

	if (inner_semi2 > 0)
	    binary_fac = Starlab::max(binary_fac, pow(inner_semi2, 5)
			     			/ pow(sister->mass, 2));

	if (binary_fac > 0)
	    peri_fac = pow(outer_peri/options->aarseth_stable_fac, 5)
				    * (1-e_outer)
		       / (binary_fac * pow(total_mass*(1+e_outer), 2));

    } else {

	// Use a simple distance criterion.  Note that, for binary-binary
	// systems, we use inner_semi_sum instead of the individual semi-major
	// axes.

	peri_fac = outer_peri / (inner_semi_sum
				   * options->unconditional_stable_fac);

    }

    // NOTE: The above definition (either version) of peri_fac omits an
    // inclination factor that effectively reduces the critical outer
    // periastron by a factor of ~ 2 for retrograde orbits relative to
    // prograde orbits.  We need an additional factor of something like
    //
    //		(4  -  2 * i / PI) / 4
    // or
    //		(3 + cos i) / 4
    //
    // where i is the inclination angle between the inner and outer
    // orbits (i = 0 for prograde, PI for retrograde).
    //
    // Not clear what to do for binary-binary systems.  Assume that the
    // inclination of the wider binary is the important quantity.

    hdyn *c = this;
    if (inner_semi2 > inner_semi1) c = sister;

    c = c->get_oldest_daughter();

    // Determine the normal vector of the inner orbit (component c).

    real cos_i = 0;

    real dx2 = square(c->get_pos());
    real dv2 = square(c->get_vel());

    if (dx2 > 0 && dv2 > 0) {
	vec n_inner = c->get_pos() ^ c->get_vel() / sqrt(dx2 * dv2);
	cos_i = n_inner * outerkep.get_normal_unit_vector();
    }

    if (options->use_aarseth_criterion)
	peri_fac *= pow(4 / (3 + cos_i), 5);
    else
	peri_fac *= 4 / (3 + cos_i);

    bool stable_a = (peri_fac > 1);

    if (stable_a) {

	// Unconditionally stable multiple.  Note that the perturbations
	// on the inner binaries are never checked explicitly.

	if (top_level)
	    multiple_type = FULL_MULTIPLE;

#if 0
	PRC(cos_i);
	if (options->use_aarseth_criterion)
	    PRL(pow(peri_fac, 0.2));
	else
	    PRL(peri_fac);
#endif

    } else {

	// Not fully stable.  Check for partial unperturbed motion.
	// NOTE: the validity of this procedure is questionable at best.

	if (top_level && options->partial_stable_fac*inner_semi_sum
					< outerkep.get_separation()) {

	    // Conditionally stable multiple (part of outer orbit only).
	    // (One of the few explicit xreal casts added by Steve, 5/00.)

	    real peri_time = (xreal)outerkep.pred_advance_to_periastron()
				- time;

	    if (peri_time > 0.8 * outerkep.get_period()) {   // 0.8 ~arbitrary

		// (Note that the above requirement actually means
		//  that the outer orbit is separating...)

		if (top_level)
		    multiple_type = APOCENTER_REFLECTION;

		stable_a = true;

	    } else {

		if (top_level && diag->report_impending_multiple_status)
		    cerr << "Outer orbit " << parent->format_label()
			 << " is partially unperturbed but not compact...\n";
	    }
	}
    }

    // Check criterion (b).

    bool stable_b = stable_a;

    if (stable_a) {

	int status_1;
	if (is_parent())
	    stable_b &= get_oldest_daughter()->is_stable(status_1, false);

	int status_2;
	if (sister->is_parent())
	    stable_b &= sister->get_oldest_daughter()
			      ->is_stable(status_2, false);

	status = 100 + 10*status_1 + status_2;
    }


    // Diagnostic output:

    if (top_level) {

	if (diag->report_impending_multiple_status
	    && !stable_a
	    && (options->unconditional_stable_fac*inner_semi_sum // non-Aarseth
			< NEAR_MULTIPLE_FAC*outer_peri		 // only
		|| options->partial_stable_fac*inner_semi_sum
			< NEAR_MULTIPLE_FAC*outerkep.get_separation())) {

	    // Getting close to the stable multiple criterion...

	    print_found_multiple(this, false, outerkep,
				 inner_semi1, inner_semi2,
				 mass, sister->mass);
	}

	if (stable_b && diag->report_multiple) {

	    if (diag->multiple_report_level > 0) {

		if (multiple_type == FULL_MULTIPLE)
		    print_found_multiple(this, true, outerkep,
					 inner_semi1, inner_semi2,
					 mass, sister->mass);
		else if (multiple_type == APOCENTER_REFLECTION) {

		    int p = cerr.precision(HIGH_PRECISION);
		    cerr << "\nfound partially stable multiple "
			 << parent->format_label()
			 << " at time " << time << endl;
		    cerr.precision(p);

		    hdyn* pnode = find_perturber_node();
		    if (!pnode)
			cerr << "pnode is NULL" << endl;
		    else {
			PRC(pnode->format_label());
			PRL(pnode->n_perturbers);
		    }

		    cerr << "outer period = " << outerkep.get_period() << endl;
		    cerr << "inner period = ";
		    if (inner_semi1 > 0)
			cerr << get_period(mass, inner_semi1) << "  ";
		    if (inner_semi2 > 0)
			cerr << get_period(sister->mass, inner_semi2);
		    cerr << endl;
		}
	    }
	}

	if (stable_a && !stable_b && diag->report_multiple) {
	    cerr << "\nmultiple " << parent->format_label()
		 << " stable, components not, at time " << time
		 << endl;
	}
    }

    return stable_b;
}



bool hdyn::is_unperturbed_and_approaching()

// Test unperturbed criterion for startup *or* continuation of unperturbed
// motion.  The two are distinguished by the existence of a kepler structure.
// Note that the notion of "unperturbed" is based on the magnitude of the
// perturbation, not on the number of perturbers on the perturber list.
// Note also that the validity of the perturber list is checked in the
// binary case, but not (yet) for multiples.

// Use of perturbation_squared is OK because the system must also be
// outside its semi-major axis to be accepted for merger.  If we relax this
// requirement, we should also normalize the perturbation to its value at
// a separation equal to the semi-major axis.

// By construction when this function is called, 'this' is already a
// low-level node, so its binary properties are well defined.

// **** Note from Steve, 7/99: this function is called at *every* perturbed
// **** step, which seems excessive -- really only need to check immediately
// **** after apocenter.  However, unperturbed periastron passages require
// **** a little more care.  Is there an efficient way to implement this?

// Function name refers to basic unperturbed criterion for simple binaries.
// For multiples, criterion is considerably more complicated...

{
    if (!options->allow_unperturbed) return false;	// fixes all binary
							// problems!

    if (diag->unpert_function_id) {
	cerr << ">> check is_unperturbed_and_approaching for "
	     << format_label() << " at time " << system_time << endl;
    }

    init_binary_type = binary_type = UNKNOWN_STATUS;	// should not occur

    // If this function returns true, then binary_type must be one of
    // the following:
    //
    //			PERICENTER_REFLECTION	// for binaries
    //			FULL_MERGER
    //			CONTINUE_MERGER
    //			UNPERTURBED_MULTIPLE_COMPONENT
    //
    //			STABLE_OUTER		// for multiples

    if (is_multiple(this)) {

	init_binary_type = binary_type = MULTIPLE_CM;

	if (!options->allow_multiples) return false;	// fixes all multiple
							// problems!!

	// Unperturbed criteria for this to be a stable system:
	//
	// (1) the outer binary is bound and stable, meaning that the
	//     perturbation is small; in addition, require that the
	//     current orbital phase be outside the semimajor axis.
	//
	// (2) the interior motion is stable, meaning that the ratio of
	//     the semimajor axis of the (or each) inner binary to the
	//     pericenter of the outer binary is less than some critical
	//     value.
	//
	// (3) each component is stable: a single star, a binary, or a
	//     stable multiple.

	multiple_type = NOT_MULTIPLE;
	int weak_stat;

	if (is_weakly_perturbed(weak_stat)) {

#if 0
	    cerr << endl
	     	 << format_label()
	    	 << " is weakly perturbed at time " << time
	    	 << "  perturbation = " << sqrt(perturbation_squared)
	    	 << endl;
#endif

	    int stable_stat;
	    if (is_stable(stable_stat)) {

		init_binary_type = binary_type = STABLE_OUTER;
		return true;

		// Note that multiples are detected only when the outer orbit
		// is advanced.  Up to the point of multiple recognition,
		// inner binaries are subject to the usual (more stringent)
		// unperturbed criterion.  With this convention, 'this' is
		// always one of the outer components of the multiple system.

	    }

	} else {

#if 0
	    cerr << endl
		 << format_label() << " is not weakly perturbed at time "
		 << time
		 << endl
		 << "perturbation = " << sqrt(perturbation_squared)
		 << "  status = " << weak_stat
		 << " (" << wp[weak_stat] << ")"
		 << endl;
#endif

	}

    } else {

	// Rest of function applies only to binaries.

	if (slow && slow->get_stop()) {

	    // Slow binary is scheduled for termination, presumably because
	    // we have already checked and expect full unperturbed motion to
	    // start soon.  No need to continue checking in this case.
	    // Return simply as a perturbed binary.

	    // if (streq(format_label(), "100a"))
	    //     cerr << "100a false slow stop..." << endl;

	    init_binary_type = binary_type = PERTURBED;
	    return false;
	}

	// The inner component of an unperturbed multiple must always
	// be regarded as unperturbed.  A side effect of this criterion
	// is that unperturbed multiples always become perturbed from
	// the top down.

	if (get_parent()->get_kepler()) {

	    // Parent node is unperturbed.

	    // if (streq(format_label(), "100a"))
	    //     cerr << "100a true unpert par..." << endl;

	    init_binary_type = binary_type = UNPERTURBED_MULTIPLE_COMPONENT;
	    return true;
	}

	bool approaching = (posvel < 0);	// posvel is recomputed at the
						// end of every perturbed step

	if (!approaching && !get_kepler() && !slow) {

	    // Accept nearly circular perturbed binaries even if they are
	    // separating slowly.  Unperturbed binaries will by construction
	    // be approaching for fully unperturbed orbits, and receding
	    // for pericenter reflection, so we shouldn't need to deal with
	    // ambiguous cases.

	    // Don't accept slow binaries in this case, as there is no
	    // possibility of promotion to full unperturbed motion and
	    // periastron reflection in a receding orbit may lead to
	    // problems.

	    // ***** Still do this at every outgoing step... (Steve, 7/99)
	    // Work with squares to avoid sqrt() hidden in vector abs().
	    // Possibly can use additional time criterion, as in slow binary
	    // motion...  However, this doesn't seem to be a significant time
	    // sink in the current code (Steve, 9/99).

	    approaching = (posvel*posvel < 1.e-6 * square(pos)*square(vel));
	}

	if (!approaching) {

	    // if (streq(format_label(), "100a"))
	    //     cerr << "100a false not appr at " << get_time() << endl;

	    init_binary_type = binary_type = NOT_APPROACHING;
	    return false;
	}

	if (!find_perturber_node()) {

	    // Possible that a newly formed center of mass node hasn't
	    // yet taken a step, so no perturber list exists.  However,
	    // we don't want to break up an unperturbed binary because
	    // of this.  Set binary_type appropriately and return true
	    // if this is already unperturbed.

	    init_binary_type = binary_type = UNKNOWN_PERTURBERS;

	    if (kep) {

		// Binary is already unperturbed; call presumably came
		// from integrate_unperturbed_motion().  Assume that some
		// reorganization has just taken place in the multiple
		// system of which the binary is a member, and allow the
		// unperturbed motion to continue for now.

		return true;

	    } else {

		// Binary is perturbed; call presumably came from kira
		// after the relative motion was computed.  In this
		// case, the perturbation should still be valid.
		// Continue checking for unperturbed motion after a
		// rudimentary consistency check.

		if (perturbation_squared < 0 || perturbation_squared > 1)
		    return false;

	    }
	}

	if (!get_kepler()
	    && options->optimize_scheduling
	    && !options->optimize_block) {

	    // Only consider new unperturbed motion at certain phases
	    // of the timestep cycle.

	    int it = (int) get_system_time();
	    real tt = get_system_time() - it;	// try to avoid overflow in it
	    real dtt = timestep/get_kappa();	// true timestep
	    it =(int)(tt/dtt + 0.1);		// should be a power of 2...

	    // Arbitrarily require that it be a multiple of 4 (so we are
	    // synchronized with steps two blocks up).  This may not really
	    // be necessary, given the synchronization strategy used below.

	    if (it%4 != 0) {
		init_binary_type = binary_type = DEFERRED;
		return false;
	    }
	}

 	// System is not a multiple, the components are "approaching," by
	// the above definition, and the perturbation is to be trusted.

	// Criterion is based on the real perturbation, not the slow-binary
	// version.  The value of perturbation_squared does *not* contain the
	// slowdown factor.

	// if (streq(format_label(), "100a"))
	//     cerr << "100a perturbation = " << perturbation_squared << endl;

	if ((perturbation_squared
	     	< gamma2 * options->partial_merge_factor)
	    && is_low_level_leaf()
	    && younger_sister->is_low_level_leaf()) {

	    // Sufficient to pick up the binary here as partly unperturbed;
	    // partial merger may be promoted to full merger later.

	    // if (streq(format_label(), "100a"))
	    //     cerr << "100a true peri refl..." << endl;

	    // PRC(perturbation_squared);
	    // PRL(gamma2 * options->partial_merge_factor);

	    init_binary_type = binary_type = PERICENTER_REFLECTION;
	    return true;

	} else {

	    real pert_fac = 1;

	    if (!kep) {

		// Use of perturbation_squared is OK here because we must
		// be outside semi for full merger...  (Alternatively,
		// we could always normalize the perturbation to separation
		// equal to semi.)

	        if ((perturbation_squared < options->full_merge_tolerance)
		    && is_low_level_leaf()
		    && younger_sister->is_low_level_leaf()) {

		    // if (streq(format_label(), "100a"))
		    //     cerr << "checking 100a..." << endl;

		    // Note that the relevant "timestep" is the step in the
		    // absence of slow motion (i.e. when unperturbed motion
		    // would actually start).

		    // if (streq(format_label(), "100a")) {
		    //     PRL(get_parent()->timestep);
		    //     PRL(10 * timestep/get_kappa());
		    // }

		    // Not entirely clear what good this timestep
		    // limit does... (Steve, 9/00)

		    if (get_parent()->timestep
			  > 10 * timestep/get_kappa()) {  // 10 is ~arbitrary

		        kepler kepl;
			hdyn *sister = get_binary_sister();

			kepl.set_time(time);
			kepl.set_total_mass(parent->get_mass());
			kepl.set_rel_pos(pos - sister->pos);
			kepl.set_rel_vel(vel - sister->vel);
			kepl.initialize_from_pos_and_vel();

			// if (streq(format_label(), "100a")) {
			//      cerr << "100a kepl..." << endl;
			//      PRL(kepl.get_energy());
			//      PRL(KEP_OUTSIDE_SEMI(kepl));
			//      PRL(get_parent()->get_next_time() - time);
			//      PRL(kepl.get_period());
			// }

			if (kepl.get_energy() < 0.0) {

			    if (KEP_OUTSIDE_SEMI(kepl)) {  // necessary if we
							   // use the current
							   // perturbation...

				// Include assumption that the parent timestep
				// will increase by at least a factor of two
				// once unperturbed motion starts (see note
				// below).

				real dtp = get_parent()->get_next_time()
				    		- time;
				if (!slow) dtp += get_parent()->get_timestep();

				// Require the period of a simple unperturbed
				// binary to fit into ~1 parent timestep.  For
				// multiples, this condition may be modified.

			        if (kepl.get_period() < dtp) {

				    // Eligible for full unperturbed motion.
				    // Defer only if slow set.

				    // Note from Steve (8/99).  This may lead to
				    // a "race" condition in the case of slow
				    // binary motion, as the parent time step in
				    // the slow case will generally be longer
				    // than in the normal case.  Thus, the slow
				    // motion may think there is enough time to
				    // allow unperturbed motion to start, but
				    // this condition may not be met once normal
				    // motion resumes.  The same issue exists
				    // when unperturbed motion starts, as the
				    // center of mass timestep will in general
				    // increase significantly.  May be time to
				    // reconsider the use of the parent time
				    // step as a limiting factor here and below.

				    // Don't permit a direct transition to
				    // unperturbed motion for slow binaries,
				    // as the slow motion has to end at the
				    // proper phase in the orbit.  Present
				    // strategy is to flag the slow motion
				    // for termination, after which unperturbed
				    // motion may begin normally.

				    if (slow) {

					if (!slow->get_stop()) {
					    cerr << endl
					    << "is_unperturbed_and_approaching:"
					    << " scheduling end of slow motion"
					    << endl
					    << "                               "
					    << " for " << format_label()
					    << " at time " << get_system_time()
					    << endl;

					    slow->set_stop();
					}

				    } else {

					// PRC(perturbation_squared);
					// PRL(options->full_merge_tolerance);

					init_binary_type = binary_type
					    		 = FULL_MERGER;

					// Ready to return true, but must still
					// check perturbation at apocenter, to
					// avoid immediate termination of
					// unperturbed motion after the first
					// step, and repetition of this cycle
					// (which leads to systematic errors).

					// Factor by which to increase
					// perturbation for apocenter
					// comparison.

					pert_fac = pow(kepl.get_apastron()
						       / kepl.get_separation(),
						       6);

					// return true;

				    }
				}
			    }
			}
		    }
		}

	    }

	    if (kep || binary_type == FULL_MERGER) {

	        // Note relaxed criterion for continuing unperturbed
	        // binary (but *not* multiple) motion.

	        real crit_pert2 = options->full_merge_tolerance
		    			 * options->relax_factor;

		// PRC(pert_fac*perturbation_squared); PRL(crit_pert2);

	        if (pert_fac*perturbation_squared < crit_pert2
		    || is_close_pair()) {

		    if (kep)
			init_binary_type = binary_type = CONTINUE_MERGER;

		    return true;

		}
	    }
	}
    }

    // if (streq(format_label(), "100a"))
    //     cerr << "100a false pert at time " << get_time() << endl;

    init_binary_type = binary_type = PERTURBED;
    return false;
}



// startup_unperturbed_motion: only invoked if is_unperturbed_and_approaching
//			       returns true...

void hdyn::startup_unperturbed_motion()
{
    // Begin treatment of unperturbed motion.  The main actions of this
    // function are (1) to create kepler structure(s) if necessary, and
    // (2) to determine the unperturbed time step(s).  May also modify
    // binary_type.

    // Note: In the case of unperturbed multiple motion, 'this' is one
    // of the outer components (see is_unperturbed_and_approaching).

    if (diag->unpert_function_id) {
	cerr << endl << ">> startup_unperturbed_motion for "
	     << format_label() << " at time " << system_time << endl;
    }

    bool new_unpert = true;
    if (get_kepler()) new_unpert = false;  // may happen if this is the inner
					   // component of a multiple system

    update_kepler_from_hdyn();	// creates a kepler structure if none exists...

    fully_unperturbed = false;	// may be set true in set_unperturbed_timestep()

    if (diag->report_start_unperturbed) {
	if (binary_type != PERICENTER_REFLECTION
	    || diag->report_pericenter_reflection) {
	    print_startup_message(this, binary_type, new_unpert);
	}
    }

    // Use of save_binary_type here may be somewhat redundant, since
    // init_binary_type seems to contain the same information...

    int save_binary_type = binary_type;	 // from is_unperturbed_and_approaching

    // Some care is needed with slow binary motion.  If the binary will
    // simply be reflected around pericenter, we will likely want to
    // continue the slow motion afterwards, so we don't want to discard
    // the slow data structures.  However, the reflection must be done
    // in tau, while set_unperturbed_timestep() works with the "real"
    // binary period, so we need to replace timestep by dtau when
    // computing steps.  If the binary would normally have been
    // "promoted" to fully unperturbed motion, we probably want to defer
    // the full startup, continue with pericenter reflection only, and
    // schedule the slow motion for termination at the proper point in
    // the orbit.

    // Temporarily replace timestep by tau in case of slow motion.

    if (slow) timestep /= get_kappa();

    // Set_unperturbed_timestep will *not* promote to fully unperturbed
    // motion if slow is set, but instead will schedule the slow motion
    // for termination.

    // Used to be an int, changed to real because integers have too
    // small a range (SPZ:02/1998):

    real steps = set_unperturbed_timestep(true);    // unperturbed timestep
						    // will be set equal to
						    // timestep * steps, in
						    // *all* cases

    // Restore timestep in case of slow motion.

    if (slow) timestep *= get_kappa();

    // Note that we may run into precision problems if we promote a very
    // eccentric binary from pericenter reflection to full merger, as
    // timestep may be very short, making 'steps' very large.

    if (diag->report_start_unperturbed && binary_type != save_binary_type) {

	if (save_binary_type == PERICENTER_REFLECTION
	    && !diag->report_pericenter_reflection)
	    print_startup_message(this, save_binary_type, new_unpert);

        cerr << "    binary_type changed to " << binary_type
	     << ":  " << bt[binary_type] << endl;

//  	if (name_is("11")) {
//  	    kep->print_all(cerr);
//  	    pp3(get_parent());
//  	}

    }

    // If steps = 0, there is no reason to start unperturbed motion.
    // Moreover, if steps is 0, program goes into an infinite loop, so
    // discard the kepler structure and return.

    // Avoid kepler if steps is less than 5 or so, in order to avoid
    // unnecessary roundoff due to kepler conversion.

    if (steps < options->min_unpert_steps) {

	if (diag->report_start_unperturbed
	    || (steps <= 0 && diag->report_zero_unpert_steps)) {
	    cerr << "\nstartup_unperturbed_motion:  calculated step size = "
		 << steps << endl
		 << "do not apply unperturbed motion to " << format_label()
		 << " at time " << system_time << "\n";
	}

	delete kep;
	kep = NULL;
	unperturbed_timestep = -VERY_LARGE_NUMBER;
	get_binary_sister()->unperturbed_timestep = -VERY_LARGE_NUMBER;
	get_binary_sister()->kep = NULL;

	if (slow)
	    slow->set_stop(false);

	return;
    }

    unperturbed_timestep = timestep * steps;
    get_binary_sister()->unperturbed_timestep = unperturbed_timestep;

    // More diagnostics.

    if (diag->report_start_unperturbed) {
	if (binary_type != PERICENTER_REFLECTION
	    || diag->report_pericenter_reflection) {

	    cerr << "    dt = " << timestep
		 << "  dt_unpert = " << unperturbed_timestep << endl;

	    if (diag->unpert_report_level > 0)
		print_unperturbed_binary_params();

#if 0

	    real predicted_mean_anomaly = kep->get_mean_anomaly()
	      + unperturbed_timestep * kep->get_mean_motion();
	    predicted_mean_anomaly = sym_angle(predicted_mean_anomaly);
	    PRI(4); PRL(predicted_mean_anomaly);

#endif

	}
    }

    // If this is the outer binary of a multiple, the inner component has
    // already been synchronized (in kira); make it unperturbed too.

    if (is_multiple(this)) {

        binary_type = STABLE_INNER;	// for startup of inner component

	if (is_parent())
	    get_oldest_daughter()->startup_unperturbed_motion();

	if (get_binary_sister()->is_parent())
	    get_binary_sister()->get_oldest_daughter()
			       ->startup_unperturbed_motion();
    }

    // Finally, replace fully unperturbed components by CM on perturber
    // lists, if necessary.

    if (!RESOLVE_UNPERTURBED_PERTURBERS && fully_unperturbed) {

	// Make a list of leaves or unperturbed nodes below parent.

	int nl = 0;
	for_all_nodes(hdyn, get_parent(), bb)
	    if (bb != get_parent()) nl++;

	hdynptr *cpt_list = new hdynptr[nl];

	// Don't know which nodes or leaves might be on others' lists,
	// so simply list all possibilities.

	nl = 0;
	for_all_nodes(hdyn, get_parent(), bb)
	    if (bb != get_parent()) cpt_list[nl++] = bb;

	correct_perturber_lists(get_root(), cpt_list, nl, get_parent());
	delete [] cpt_list;
    }

#if 0

    // Keep track of critical orbital data for unperturbed systems
    // in relatively wide triple orbits.		(Steve, 8/03)

    // *** NO LONGER NEEDED, but retain the code for now as an example.
    // *** See integrate_unperturbed_motion for an explanation of how
    // *** the tidal error arises and how to correct it in place.

    // *** NOTE that, if used, this approach would require writing the
    // *** extra data to the output file to avoid restart problems.

    if (fully_unperturbed) {

	hdyn *par = get_parent(), *pnn = par->get_nn();

	if (pnn) {

	    vec ppos = hdyn_something_relative_to_root(par,
						       &hdyn::get_pred_pos);
	    vec pvel = hdyn_something_relative_to_root(par,
						       &hdyn::get_pred_vel);
	    vec npos = hdyn_something_relative_to_root(pnn,
						       &hdyn::get_pred_pos);
	    vec nvel = hdyn_something_relative_to_root(pnn,
						       &hdyn::get_pred_vel);

	    real m12 = par->get_mass();
	    real m3 = pnn->get_mass();
	    real m123 = m12 + m3;
	    real mu = m12 * m3 / m123;
	    real R = abs(npos - ppos);
	    real Vsq = square(nvel - pvel);

	    real E = mu*(0.5*Vsq - m123/R);

	    if (E < 0) {

		// Need kT for comparison.  Anticipate doing this rarely, so
		// compute it the hard way.  (Could save kT in the story...)

		int ntop = 0;
		real kin = 0;
		for_all_daughters(hdyn, get_root(), bb) {
		    ntop++;
		    kin += bb->get_mass()*square(bb->get_pred_vel());
		}
		real kT = 2*kin/(3*ntop);

		if (E < -0.1*kT) {		// conservative -- may improve

		    // Save orbital data.

		    kep_config *k = kep->get_kep_init();
		    if (k) delete k;
		    k = new kep_config;
		    k->time = system_time;
		    k->nn = (void*)pnn;
		    k->energy = E;

		    vec pos1 = hdyn_something_relative_to_root(this,
						       &hdyn::get_pos);
		    hdyn *sis = get_younger_sister();
		    vec pos2 = hdyn_something_relative_to_root(sis,
						       &hdyn::get_pos);
		    real r13 = abs(npos - pos1);
		    real r23 = abs(npos - pos2);
		    real m1 = get_mass();
		    real m2 = sis->get_mass();

		    // Save the tidal potential energy as the energy error.

		    k->de = -m3 * (m1/r13 + m2/r23 - m12/R);

		    kep->set_kep_init(k);

		    if (diag->report_start_unperturbed) {
			int p = cerr.precision(HIGH_PRECISION);
			cerr << endl
			     << "startup_unperturbed_motion: "
			     << "saving neighbor data for "
			     << parent->format_label()
			     << endl
			     << "    at time " << system_time
			     << ",  outer ";
			cerr.precision(p);
			real r12 = abs(get_pos()-sis->get_pos());
			PRC(E/kT); PRL(k->de);
			PRI(4); PRC(r12); PRL(R);
		    }
		}
	    }
	}
    }
#endif

    // cerr << "Leaving startup_unperturbed_motion: "; PRL(fully_unperturbed);

    // On normal exit, a kepler structure exists, timestep is the last
    // perturbed timestep, suitable for restart if no orbital evolution
    // occurs, and unperturbed_timestep is an integral number of timesteps.
    // The hdyn is frozen at the instant at which unperturbed motion was
    // started or continued.  It is updated at each unperturbed step.

    // *** The only difference between pericenter reflection and full ***
    // *** unperturbed motion is the value of fully_unperturbed.      ***

    // For unperturbed multiple motion, each unperturbed component is
    // scheduled and updated separately.  Synchronization is maintained
    // by choice of time steps, but update times will in general *not*
    // be the same for all levels.

#if 0
    putrq(get_parent()->get_dyn_story(), "unpert_startup_time",
	  get_system_time());
    int ifull = fully_unperturbed;
    putiq(get_parent()->get_dyn_story(), "fully_unperturbed", ifull);
#endif

}



real hdyn::set_unperturbed_timestep(bool check_phase)	// no default

// Use of the argument check_phase:
//
// Argument check_phase = true indicates that we should check that the
// system is outside its semi-major axis (using KEP_OUTSIDE_SEMI)
// before applying full merging.
//
// This function is presently called only from
//
//	startup_unperturbed_motion()		(this file)
//	integrate_unperturbed_motion()		(this file)
//
// In the first two cases, the variable fully_unperturbed is false
// on entry, and may be set true here.  A call from anywhere else
// will generally be made to update an existing unperturbed binary
// (e.g. after binary evolution), and fully_unperturbed will be true.
// (It is important to distinguish these possibilities, as the binary
// will not necessarily be approaching in the latter cases.)
//
// Note that this function has accidentally been given the same
// name as the (real-argument) member function that actually sets
// the variable unperturbed_timestep for the hdyn class (oops!)...
//
// Argument check_phase is set true when this function is called from
// startup_unperturbed_motion().  It is set to !force_time when the
// call comes from integrate_unperturbed_motion().  The default value
// for force_time is false (kira call uses this default).  It is
// generally undesirable to set force_time to be true; however, it may
// be set true by a call from integrate_node() if its arguments
// integrate_unperturbed_motion and force_unperturbed_time are both
// true (default: true, false).  Integrate_node() seems only to be
// called from synchronize_node(), with values false and true.
//
// *** Thus, it appears that check_phase is *always* true when this
// *** function is called!  We thus promote pericenter reflection to
// *** full merger only when the orbit is "near" apocenter.  Hmmm...

{
    // Determine and set the unperturbed timestep.  Zero return value
    // means no unperturbed motion.
    //
    // *** May change the value of binary_type. ***
    //
    // On arrival in this function, binary_type must be one of the
    // "true" values set by is_unperturbed_and_approaching.  The
    // MULTIPLE values get special treatment.  The "binary" values,
    // in order of *increasing* external perturbation, are:
    //
    //		PERICENTER_REFLECTION
    //		FULL_MERGER
    //		CONTINUE_MERGER
    //		UNPERTURBED_MULTIPLE_COMPONENT
    //
    // Since the binary has already been deemed "unperturbed", an
    // unperturbed step should be taken if possible, even if it
    // extends beyond the parent step.

//    cerr << "in set_unperturbed_step for " << format_label()
//	 << " at " << time << "/" << system_time << endl;

    if (diag->unpert_function_id) {
	cerr << ">> set_unperturbed_timestep for "
	     << format_label() << endl;
    }

    if (!check_phase) {
        // cerr << "set_unpert_timestep called with no test phase\n";
        // PRL(kep->get_energy());
    }

    hdyn * sister = get_binary_sister();
    real steps = 0;

    // Establish "default" outcome.

    if (is_multiple(this)) {

	// No default -- unperturbed motion does not occur if the
	// criteria in get_unperturbed_steps() are not met.

    } else {

	// Default is periastron reflection.  Determine the number of
	// perturbed steps to the reflection point.

	binary_type = PERICENTER_REFLECTION;

	// real peri_time = kep->pred_advance_to_periastron() - time;

	// Old code used pred_advance_to_periastron(), but this may be quite
	// inaccurate.  Better to use mean motion.  (Steve, 4/99)

	real mean_anomaly = kep->get_mean_anomaly();
	if (kep->get_energy() < 0)
	    mean_anomaly = sym_angle(mean_anomaly);	// place in (-PI, PI]

	if (mean_anomaly >= 0) {

	    // Note that we always expect mean_anomaly < 0, so this
	    // should never happen:

	    steps = 0;

	} else {

	    real peri_time = -mean_anomaly / kep->get_mean_motion();
	    steps = ceil(2 * peri_time / timestep);

	    // Recall that ceil rounds to the next integer up, so we
	    // overshoot slightly.  We *don't want to do this for slow
	    // binaries, as we must terminate (components receding) at
	    // the end of the step.

	    if (slow) steps -= 1;

	}
    }

    if (kep->get_energy() < 0) {

	// Consider the possibility of complete merging in a bound orbit.

	if (is_multiple(this) && multiple_type == APOCENTER_REFLECTION) {

	    // Determine if we are close to apocenter in a multiple.

	    // Note: pred_advance_to_periastron() may not yield a very
	    // accurate time...
	    // (One of the few explicit xreal casts added by Steve, 5/00.)

	    real apo_time = (xreal)kep->pred_advance_to_periastron()
	                          - time - (xreal)(0.5*kep->get_period());
	    real steps = floor(((2 * apo_time) / timestep));

	    // Note: "floor" here means the separation at the end of
	    //	      the step is slightly greater than at the start.

	    if (diag->report_multiple) {
		cerr << "multiple (" << tt[multiple_type] << "): "
		     << format_label() << endl;
		PRC(apo_time), PRC(kep->get_period()); PRL(steps);
	    }

	    binary_type = STABLE_OUTER;

	} else if (!check_phase
		   || KEP_OUTSIDE_SEMI(*kep)	     // also picks up multiple
		   || (get_parent()->kep != NULL)) {

	    // Full merger criteria:
	    //     (1) not checking phase (discouraged), or
	    //     (2) system is outside its semi-major axis, or
	    //     (3) parent is already unperturbed.

	    // At this point, we would normally promote the binary to fully
	    // unperturbed motion.  However, if slow motion is in progress,
	    // don't allow that, but schedule the slow motion for termination
	    // (to catch the promotion next time around).

	    if (slow) {

		if (!slow->get_stop()) {
		    cerr << "set_unperturbed_timestep (#1): "
			 << "scheduling end of slow motion"
			 << endl
			 << "                               for "
			 << format_label()
			 << " at time " << get_system_time() << endl;
		    slow->set_stop();
		}

		// See note on the use of the parent time step in function
		// is_unperturbed_and_approaching().

	    } else {

		// The criteria in get_unperturbed_steps() will return zero
		// if it is not possible to fit an integral number of orbits
		// or an advance to apastron into the time allowed by the
		// parent step.  If usteps <= 0, this must mean that the timing
		// of the parent step is the problem.  In that case, we fall
		// back to pericenter reflection.  However, if that occurs
		// near apocenter (as it will if the orbit had previously been
		// advanced to apastron), we convert it to a complete orbit.

		// If usteps > 0, for a true binary, we "overshoot" slightly,
		// ending at the next perturbed step past apocenter.  For a
		// multiple, we don't proceed to apocenter.  Instead, we
		// require an integral number of orbits ending just before
		// (i.e. at separation outside) a complete period multiple.

		// Note that this is one of only two places in kira where
		// get_unperturbed_steps() is called.  (The other is in
		// util/hdyn_kepler.C.)

		real usteps = get_unperturbed_steps(!is_multiple(this));

		// Unperturbed timestep will be set equal to timestep * usteps
		// unless pericenter reflection is promoted to full merger.

		// PRL(usteps);

		if (usteps <= 0) {

		    if (diag->report_zero_unpert_steps) {

			cerr << endl << "    set_unperturbed_timestep: "
			     << "zero step for unperturbed binary\n";

			int p = cerr.precision(HIGH_PRECISION);
			PRI(4); PRL(time);
			cerr << "    parent: " << get_parent()->format_label()
			     << endl;
			PRI(4); PRL((get_parent()->get_next_time()));
			PRI(4); PRL(kep->get_period());
			PRI(4); PRL(usteps);
			cerr.precision(p);
		    }

		    // Note that steps and binary_type are left unchanged.

		} else {

		    steps = usteps;
		    binary_type = FULL_MERGER;	// may be changed from
		    				// CONTINUE_UNPERTURBED
		    fully_unperturbed = true;

		}
	    }
	}

#if 0

	real predicted_mean_anomaly = kep->get_mean_anomaly()
	  + steps * timestep * kep->get_mean_motion();
	predicted_mean_anomaly = sym_angle(predicted_mean_anomaly);
	PRI(4); PRL(predicted_mean_anomaly);

#endif

	// Now steps is the number of perturbed time steps to the reflection
	// point (rounded to the nearest integer), or the number of perturbed
	// steps corresponding to an integral number of orbits (plus mapping
	// to apocenter for true binaries) if full merging is applied.  Note
	// that, in either case, the resulting unperturbed step is *not*
	// necessarily a power of 2, but it does maintain (or improve) the
	// binary's place in the block scheduling scheme if and when perturbed
	// motion resumes.

	// ** Could improve scheduling by using time step near apocenter, not
	// ** promoting pericenter reflection to full unperturbed motion, and
	// ** restricting unperturbed motion to a power of two perturbed steps.

	// Special case occurs if the reflection is occurring near apocenter.
	// In this case, promote the step to a full merger for a single orbit.

	// Hmmm... Not altogether clear why we consider this separately here...
	//							(Steve, 8/99)

	if (binary_type == PERICENTER_REFLECTION && KEP_OUTSIDE_SEMI(*kep)) {

	    // Same treatment of slow binaries as above.

#if 0
	    // Flag this to see if it ever occurs (Steve, 9/00).

	    cerr << endl
		 << "**** binary_type == PERICENTER_REFLECTION "
		 << "&& KEP_OUTSIDE_SEMI(*kep) ****"
		 << endl << endl;
#endif

	    if (slow) {

		if (!slow->get_stop()) {
		    cerr << "set_unperturbed_timestep (#2): "
			 << "scheduling end of slow motion"
			 << endl
			 << "                               for "
			 << format_label()
			 << " at time " << get_system_time() << endl;
		    slow->set_stop();
		}

		// See note on the use of the parent time step in function
		// is_unperturbed_and_approaching().

	    } else {
		steps = ceil(kep->get_period()/timestep);
		binary_type = FULL_MERGER;
		fully_unperturbed = true;
	    }
	}
    }

    // Redundant comment: for slow binary motion, timestep is temporarily
    // dtau, so multiply unperturbed_timestep by kappa here (value is
    // actually overwritten in startup_unperturbed_motion anyway...).

    unperturbed_timestep = timestep * steps;
    if (slow) unperturbed_timestep *= get_kappa();

    sister->kep = kep;
    sister->unperturbed_timestep = unperturbed_timestep;
    sister->fully_unperturbed = fully_unperturbed;

    if (steps < 0) {
	cerr << endl
	     << "set_unperturbed_timestep: steps negative!"
	     << " -- the code will break soon...;-("
	     << endl;
	PRC(timestep); PRC(steps); PRL(unperturbed_timestep);
	kep->print_all();
    }

    return steps;
}



#define DEBUG_SCHEDULE
#undef DEBUG_SCHEDULE

#ifdef DEBUG_SCHEDULE
bool debug_flag = false;
#endif

// Careful!  Debugging causes problems with optimization under
// RH Linux 7.2 with g++ 2.96...  Adding/cutting lines can change
// the results and lead to incorrect output.	(Steve, 7/02)
//
// Dec(Compaq) OS4.0/gcc-whatever: Seems fine.
//
// RH 7.2/gcc-2.96: Works with inline below, fails without (-O2).
//		    Optimization -O1 seems to be OK...

local
//inline
xreal latest_time(xreal t_min, xreal t_max, real dtblock,
		  real ma_min, real ma_max,
		  real mean_anomaly, real mean_motion)
{
    // Determine the latest time (integer * dtblock) that lies between
    // t_min and t_max and has mean anomaly in the range (ma_min, ma_max).
    // Note: input mean_anomaly refers to time t_min.
    //
    // For use with optimized binary scheduling (Steve, 6/00).

    // On entry, should have -PI <= ma_min < ma_max < 0.

    if (ma_min < -PI) return t_min;
    if (ma_max >= 0) return t_min;
    if (ma_min >= ma_max) return t_min;

    real f = floor((real)t_max / dtblock);
    xreal t_last = dtblock * f;			// target time on this block

#ifdef DEBUG_SCHEDULE
    if (debug_flag) {
        PRC(t_last); PRL(t_min);
    }
#endif

    if (t_last <= t_min) return t_min;

    // t_last is the last block time in the allowed range.
    // Determine the mean anomaly at that time.

    real ma = sym_angle(mean_anomaly + mean_motion * (real)(t_last - t_min));

#ifdef DEBUG_SCHEDULE
    if (debug_flag) {
        PRC(ma); PRC(ma_min); PRL(ma_max);
    }
#endif

    if (ma > ma_min && ma < ma_max) return t_last;

    // Determine (modulo 2 PI) the change in ma per dtblock.

    real dma = sym_angle(dtblock * mean_motion);

#ifdef DEBUG_SCHEDULE
    if (debug_flag) {
        PRL(dma);
    }
#endif

    if (dma == 0) return t_min;

    // Find how many backward steps of dtblock we must take for
    // mean_anomaly to move into the desired range.  Simplify the
    // calculation by placing the lower limit of the range at ma_min
    // effectively moving all quantities into the range (0, 2 PI).

    ma -= ma_min;
    ma_max -= ma_min;
    ma_min = 0;

    // Now we want to find a time t_last corresponding to 0 < ma < ma_max.

#ifdef DEBUG_SCHEDULE
    if (debug_flag) {
	PRC(t_min); PRC(t_max); PRL(dtblock);
	PRC(ma_max); PRC(ma/(2*M_PI)); PRL(dma/(2*M_PI));
    }
#endif

    // Deal with the easy cases first...

    if (dma < 0 && dma > -ma_max) {

	// A single backward step will increase ma by 0 < -dma < ma_max.
	// By construction, we now have 0 < -dma < ma_max < ma < 2*M_PI.

	t_last -= dtblock * ceil((2*M_PI - ma) / (-dma));

    } else if (dma > 0 && dma < ma_max) {

	// A single backward step will decrease ma by 0 < dma < ma_max.
	// In this case, 0 < dma < ma_max < ma < 2*M_PI.

	t_last -= dtblock * ceil((ma - ma_max) / dma);

    } else {

	// Do it the hard way (step by step), for now.

#ifdef DEBUG_SCHEDULE
	if (debug_flag) {
	    cerr << "the hard way..." << endl << flush;#endif

	}
#endif

	while (t_last > t_min) {
	    t_last -= dtblock;
	    ma -= dma;
	    if (ma <= 0) ma += 2*M_PI;
	    if (ma > 2*M_PI) ma -= 2*M_PI;
//	    PRC(t_last); PRL(ma);
	    if (ma < ma_max) break;
	}
    }

#ifdef DEBUG_SCHEDULE
    if (debug_flag) {
	PRL(t_last);
    }
#endif
    return t_last;
}

real hdyn::get_unperturbed_steps(bool to_apo,	// default true (for binary)
						//   -- set false for multiple
				 bool predict)	// default false; *never* used

// Additional (first) argument to_apo added by Steve 8/5/98.
//
// Note from SPZ.  Previous implementation returned an int.
// Stellar evolution allows binaries with orbital periods as
// small as a few milliseconds (see Portegies Zwart 1998), so
// return type was changed to real.
//
// int hdyn::get_unperturbed_steps(bool predict, bool to_apo)
// real hdyn::get_unperturbed_steps(bool predict, bool to_apo)
// unsigned long hdyn::get_unperturbed_steps(bool predict, bool to_apo)
//
// This function is called only from
//
//	set_unperturbed_timestep()		(this file)
//	reinitialize_kepler_from_hdyn()		(hdyn_kepler.C)
//	evolve_stars()				(kira_stellar.C)
//
// It is sometimes more convenient to bypass the normal checks applied
// by set_unperturbed_timestep() when we already know the binary
// configuration.

{
    // Return the number of current time steps to advance
    // unperturbed binary motion up to (but not too far past)
    // the next CM time step.

    // If to_apo is true (default), the step will continue on
    // to the next apocenter (before the CM step).

    // First check if parent is younger binary sister, in which
    // case its time step may not be "definitive."

    // Note from Steve, 7/17/97:  For unknown reasons, this
    // function fails to compile properly under g++ version
    // cygnus-2.7-96q4 on Dec UNIX V4.0B, at optimization
    // levels higher than 0...

    // *** NOTE: Do not assume approaching components, or ***
    // ***	 time = system_time or parent time...	  ***


    if (!kep) return 0;


//    cerr << "in get_unperturbed_steps for " << format_label()
//	 << " at " << time << "/" << system_time << endl;

    hdyn * p = get_parent();
    if (p->is_low_level_node() && p->get_elder_sister() != NULL)
	p = p->get_elder_sister();

    // Set a limit on the next unperturbed step.

    real pdt = p->get_next_time() - time;	// time to end of parent time
						// step; note that time may
						// differ from system_time

    if (pdt > unpert_step_limit) {
	if (diag->report_continue_unperturbed)
	    cerr << "get_unperturbed_steps: unperturbed step for "
		 << format_label()
		 << " limited by unpert_step_limit" << endl;
	pdt = unpert_step_limit;
    }

    // Next step must end after current system time.  (Not relevant
    // during normal unperturbed step, but needed when recomputing step
    // after asynchronous binary evolution.)

    if (pdt < system_time - time) {
	if (diag->report_continue_unperturbed)
	    cerr << "get_unperturbed_steps: unpert step for "
		 << format_label()
		 << " increased to reach system_time" << endl;
	pdt = system_time - time;
    }
	
    if (pdt < 0) {
	cerr << endl << "get_unperturbed_steps: ";
	PRC(p->get_next_time()); PRL(time);
	pp3(p);
	cerr << "*no* corrective action taken\n";
	return 0;
    }

    real pdt2 = pdt + dt_overshoot(this);

    // Special treatment of multiple motion, since the cost of not
    // starting unperturbed motion is so high...

    if (USE_DT_PERTURBERS && is_multiple(this)) {

	// Include perturber crossing time in pdt...

	real pert_dt = dt_perturbers(this);
	if (pert_dt > 0) pdt2 = Starlab::max(pdt2, 0.25*pert_dt);	// conservative
    }

    // Goal: to advance the binary by as great a time as possible,
    //	     subject to the constraints that (a) we do not wish to
    //	     go too far past the parent node's step, and (b) we do
    //	     not wish to end up with a separation smaller than the
    //	     current separation.  In the event of a conflict, it is
    //	     better to break constraint (a) than constraint (b), as
    //	     we will continue with the current perturbed time step
    //	     on restart.  The actual unperturbed timestep will be
    //	     an integral multiple of the perturbed step, to retain
    //	     some semblance of block scheduling.
    //
    //	     **** We do, however, wish to end up after the parent
    //	     **** step, and not just before it (Steve, 9/99).
    //
    //  (1)  If to_apo is false (e.g. for unperturbed multiples, to
    //       minimize the tidal error), then we simply want to
    //       advance the system by an integral number of orbits.  If
    //       this is not possible within the parent step, or within
    //       half a parent time step after the present parent step,
    //       then we return zero steps (thus preventing unperturbed
    //       motion from starting or continuing).  The unperturbed
    //       step will *never* exceed the parent step by more than 1
    //       orbital period.
    //
    //       Since the unperturbed time step is constrained by the
    //       value of the perturbed step, it is not in general
    //       possible to advance a binary by an exact number of orbit
    //       periods.  Unperturbed motion is initiated only on the
    //       inward part of the binary orbit, so we take the time
    //       step to be slightly *less* than a whole number of
    //       periods, ensuring that the final separation is slightly
    //       greater than the initial one.  We also ensure that the
    //       step will always pass the next apocenter (assuming that
    //       a step is to be taken), even if that will violate either
    //       or both of constraints (a) and (b) above.
    //
    //  (2)  If to_apo is true (e.g. for true binaries, where we will
    //       neglect the tidal error associated with changing the
    //       binary phase), we will also try to map the orbit to (just
    //       after) the apocenter preceding the parent's next step.
    //       If necessary, we may advance to the apocenter beyond the
    //       parent time step, subject to the same restrictions as in
    //       (1) above.

    // Determine the integral number of orbits to take.

    // (Recall that floor rounds down to the next integer, so orb
    //  orbital periods end *before* pdt.)

    // real orb = floor(pdt / kep->get_period());
    real orb = ceil(pdt / kep->get_period());

    // (Recall that floor rounds down to the next integer, so orb
    //  orbital periods end *after* pdt.)

    if (orb <= 0) {

	// See if we can take an orbital step in less than half the
	// next parent step or the perturber crossing time (pdt2).

	if (kep->get_period() <= pdt2)
	    orb = 1;

	// Note that, if we can't fit into pdt but we can into pdt2,
	// we conservatively take a step of only one period, regardless
	// of how many periods would actually fit into pdt2.
    }

    // Get the normalized mean anomaly.  If this function is invoked as part
    // of a normal step, we should have kepler time = time = system_time, so
    // this is unambiguous.  However, if this function is called from
    // elsewhere, care must be taken to define the time.  The first use of
    // mean_anomaly below is in determining the time to apastron, so the
    // time in question should be the time for this particle.  Adopt this
    // convention, and check later when some other time is implied.

    real mean_anomaly = kep->get_mean_anomaly();
    real mean_motion =  kep->get_mean_motion();
    if (kep->get_time() != time)
	mean_anomaly += ((real)(time - kep->get_time())) * mean_motion;
    mean_anomaly = sym_angle(mean_anomaly);		// in (-PI, PI]

    // Determine number of perturbed steps in the next unperturbed step.

    real pert_steps = 0;

    if (to_apo) {

	// Get time to next apocenter.

	// Old code used pred_advance_to_periastron(), but this may be quite
	// inaccurate.  Better to use mean motion.  (Steve, 4/99)

	// Note that mean_anomaly here must correspond to the base time
	// of the unperturbed step.

	real apo_time = (M_PI - mean_anomaly) / mean_motion;

	// Note added by Steve, 7/99.  If the unperturbed motion has for
	// some reason (e.g. rounding error) managed to end just before
	// apocenter, i.e. with mean_anomaly > 0), then apo_time could be
	// very short, possibly as small as 1 perturbed step.  In this case,
	// go to the *next* apocenter, in order to avoid very short steps.
	// (This is what we would do here if the step ended properly, with
	// mean_anomaly < 0.)

	if (mean_anomaly > 0.9*M_PI
	    && kep->get_period() < pdt2) apo_time += kep->get_period();

	pert_steps = ceil( (Starlab::max(0.0, orb - 1) * kep->get_period() + apo_time)
		     						/ timestep
		     + 1);	// extra "1" to try to overcome
				// possible problems with roundoff

	// Use of ceil here ensures that the unperturbed step ends just after
	// the apocenter immediately preceding the last allowed full orbit.
	// If orb = 0, this step takes us to the next apocenter, regardless
	// of the parent step.

    } else {

	pert_steps = floor(orb * kep->get_period() / timestep);

	// Use of floor here ensures that the unperturbed step ends just
	// before an integral number of orbits.

	// However, this means that the binary separation at the end point
	// increases slowly (if we start with mean_anomaly < 0).  Must make
	// sure that the step really ends after apocenter (and not just
	// before).
    }

    // Recheck that step will advance beyond system_time.  By construction,
    // it should fall short by at most 1 period...

    if (time + pert_steps*timestep < system_time) {
	if (diag->report_continue_unperturbed)
	    cerr << "get_unperturbed_steps: unpert step for "
		 << format_label()
		 << " increased to reach system_time" << endl;
	pert_steps += ceil(kep->get_period() / timestep);
    }

    // Check that we really will pass apocenter at end of step.

    real d_mean_anomaly = timestep * mean_motion;
    real end_mean_anomaly = sym_angle(mean_anomaly
				       + pert_steps * d_mean_anomaly);
    int  mcount = 0;

    if (end_mean_anomaly > 0) {

	// Mean anomaly should be negative for incoming motion at end of step.
	// A while loop here could be dangerous.  Better to check iteration
	// count as we go.

	mcount = (int)(kep->get_period() / timestep);

	while (end_mean_anomaly <= M_PI && mcount > 0) {
	    end_mean_anomaly += d_mean_anomaly;
	    pert_steps += 1;
	    mcount--;
	}

	if (end_mean_anomaly >= M_PI) end_mean_anomaly -= 2*M_PI;

	if (false && mcount <= 0)
	    cerr << "get_unperturbed_steps: mcount = 0 for "
		 << format_label() << " at time " << time << endl;
    }

    //-----------------------------------------------------------------

    // We have now determined the best unconstrained value of pert_steps.
    // Optionally try to modify the choice to improve scheduling.

    if (pert_steps > 0 && options->optimize_scheduling) {

	// Best (maximum) value for the number of perturbed steps to take
	// is pert_steps.  Try to force the actual step into the block scheme
	// in the best possible location for scheduling purposes.
	//
	// ~Arbitrary criterion: step must end after apocenter, but not
	// too far after it -- require -PI < mean_anomaly < -0.9 PI, say.

	// Note that this procedure overrides and largely ignores
	// the previous timestep determination...

	if (!options->optimize_block) {

	    // Try to make pert_steps a multiple of the largest possible
	    // power of 2.
	    //
	    // *** Note that the "base" (perturbed) step may be quite short,
	    // *** so this may not actually help with overall synchronization
	    // *** (Steve, 6/00)...

	    // First determine the amount of "wiggle room."

	    int minus = (int)((end_mean_anomaly + M_PI) / d_mean_anomaly);
	    int plus = (int)((-0.9*M_PI - end_mean_anomaly) / d_mean_anomaly);

	    // Limits should be unnecessary, but...

	    if (minus > pert_steps/2) minus = (int)(pert_steps/2);
	    if (plus > pert_steps/2) plus = (int)(pert_steps/2);

	    // Desired range is  (pert_steps - minus)  to  (pert_steps + plus).

	    real new_steps = pert_steps - minus;
	    real p2 = 2;
	    while (new_steps <= pert_steps + plus) {

		if (fmod(new_steps, p2) != 0) new_steps += p2/2;

		// Now new_steps is a multiple of p2.

		p2 *= 2;
	    }

	    // On exit from the loop, new_steps is too big -- reduce it here.

	    p2 /= 2;
	    new_steps -= p2/2;

	    // PRC(new_steps); PRL(p2);

	    pert_steps = new_steps;

	} else {

	    // NEW STRATEGY (Steve, 9/99, 6/00).  Choose the number of
	    // perturbed steps in order to optimize the "block ranking" of
	    // the unperturbed motion, i.e. the motion will be synchronized
	    // with the highest block possible at the end of the next step.
	    // OK to take a step much smaller than the optimal pert_steps,
	    // so long as the overall synchronization is improved.  This
	    // strategy should allow the motion to migrate to and remain in
	    // a tolerably high block (in terms of "get_next_time") after a
	    // few unperturbed steps.
	    //
	    // Schematic procedure:
	    //
	    //   1. Start at top block (dt = 1, nb = 0) and work down (dt
	    //      /= 2; nb++) until dt = current (perturbed) timestep.
	    //
	    //   2. Within each block, determine the latest time (t =
	    //	    integer * 2^{-nb}) that lies within acceptable limits
	    //	    -- i.e. within the allowed time range and within the
	    //	    permitted phase window.
	    //
	    //   3. If such a time exists, accept it unconditionally as
	    //	    next_time and set pert_steps accordingly.

	    // Define range of acceptable end-times.  Allow the step to
	    // end up to half a parent timestep past the next parent step.
	    // Choice of t_min means that we require an "optimized" step
	    // of at least one orbit before we accept it.

	    // *** Careful with t_min choice in case of approaching binary ***
	    // *** TO DO...

	    xreal t_min = system_time + kep->get_period();	// note use of
								// system_time
	    real mean_anomaly_min = mean_anomaly
					+ ((real)(t_min-time))*mean_motion;

	    real pdt = p->get_timestep();
	    xreal t_max = p->get_next_time() + 0.5*pdt;

	    // Define range of acceptable final mean anomalies.

	    real ma_min = -0.999*M_PI;
	    real ma_max = -0.750*M_PI;		// was -0.900 (Steve, 5/02)

	    // Step will end at time t_next in the block defined by dtblock
	    // if t_next > t_min.

	    xreal t_next = t_min;
	    real dtblock = Starlab::min(unpert_step_limit, 2*pdt);
	    int kb = get_effective_block(dtblock);

#ifdef DEBUG_SCHEDULE
	    debug_flag = false;
	    if (name_is("11032") || name_is("13993")) debug_flag = true;
#endif

	    while (dtblock >= timestep) {
		t_next = latest_time(t_min, t_max, dtblock,
				     ma_min, ma_max,
				     mean_anomaly_min,
				     mean_motion);
		if (t_next > t_min) break;
		kb++;
		dtblock /= 2;
	    }

	    if (t_next <= t_min
		|| kb >= get_effective_block(time+pert_steps*timestep)) {

		// Nothing to be gained from optimizing.  Don't
		// change pert_steps.

#ifdef DEBUG_SCHEDULE
		cerr << "    retaining unoptimized pert_steps" << endl;
#endif

	    } else {

		real old_steps = pert_steps;
		pert_steps = floor(((real)(t_next - time)) / timestep + 0.1);

#ifdef DEBUG_SCHEDULE

		// (These debugging lines can be quite expensive...)

		if (debug_flag) {
//		if (system_time > 0.112) {

		    cerr << endl << "get_unperturbed_steps for "
			 << format_label() << " at time " << system_time << ":"
			 << endl;
		    int pp = cerr.precision(20);

		    PRI(4); PRL(system_time);
		    PRI(4); PRL(time);
		    PRI(4); PRL(timestep);
		    PRI(4); PRL(get_effective_block(time));
		    PRI(4); PRL(get_effective_block(timestep));

		    PRI(4); PRL(old_steps);
		    PRI(4); PRL(time+old_steps*timestep);
		    PRI(4); PRL(get_effective_block(old_steps*timestep));
		    PRI(4); PRL(get_effective_block(time+old_steps*timestep));
		    PRI(4); PRL(sym_angle(mean_anomaly
				+ mean_motion*old_steps*timestep));

		    PRI(4); PRC(ma_min); PRL(ma_max);
		    PRI(4); PRC(kb); PRC(dtblock); PRL(dtblock/timestep);
		    PRI(4); PRC(t_next); PRL(t_min);
		    PRI(4); PRC(pert_steps); PRL(pert_steps/old_steps);
		    PRI(4); PRL(get_effective_block(pert_steps*timestep));
		    PRI(4); PRL(get_effective_block(time+pert_steps*timestep));
		    PRI(4); PRL(sym_angle(mean_anomaly
				+ mean_motion*pert_steps*timestep));

		    if (get_effective_block(time+pert_steps*timestep) != kb)
			cerr << "    ????" << endl;

		    cerr.precision(pp);



	    real predicted_mean_anomaly = kep->get_mean_anomaly()
	      + pert_steps * timestep * kep->get_mean_motion();
	    predicted_mean_anomaly = sym_angle(predicted_mean_anomaly);
	    PRI(4); PRL(predicted_mean_anomaly);





		}
#endif

	    }

#if 0
	    cerr << endl; PRC(t_next); PRC(t_min); PRL(system_time);
	    PRC(pert_steps); PRL(pert_steps*timestep);
#endif
	}
    }

    return pert_steps;		// unit = current unperturbed timestep
}

void hdyn::recompute_unperturbed_step()
{
    //	    cerr << "Recomputing unperturbed step for "
    //		 << "(" << b->format_label() << ",";
    //	    cerr << b->get_binary_sister()->format_label() << ")\n";
    //	    PRL(b->get_fully_unperturbed());

    // Ordinarily, if the motion is flagged as unperturbed on restart,
    // the unperturbed timestep should be trustworthy.  However, it
    // is possible that the node time has been synchronized by some
    // program, and that time + unperturbed_timestep may not take us
    // to apastron, as is desirable (and assumed?) elsewhere in kira.
    // Recompute the step, using the standard criteria.  (SLWM, 4/02)

    // int usteps = get_unperturbed_steps(true);
    // unsigned long usteps = get_unperturbed_steps(true);
    real usteps = get_unperturbed_steps(true);

    if (usteps > 0) {
	unperturbed_timestep = timestep*usteps;
	get_binary_sister()->unperturbed_timestep
	    = unperturbed_timestep;
    } else {
	// Just live with this, and assume that it will be handled
	// properly once the step is over...
    }
}

void hdyn::recompute_unperturbed_steps()
{
    // Use time = system time as an indicator that the binary may
    // have beeen synchronized improperly...

    for_all_nodes(hdyn, this, b)
        if (b->is_low_level_node()
	    && !b->get_elder_sister()
	    && b->get_unperturbed_timestep() > 0
	    && b->time == system_time)			// new 07/03 (Steve)

	    b->recompute_unperturbed_step();
}



local vec rotate(vec v, vec r[3])
{
    return vec(v*r[0], v*r[1], v*r[2]);
}

local void rotate_kepler(kepler *k)
{
    // Apply a random rotation to all vectors of kepler k.

    real alpha = randinter(-M_PI, M_PI),
	 beta  = randinter(-M_PI, M_PI),
	 gamma = randinter(-M_PI, M_PI);		// Euler angles
    real ca = cos(alpha), sa = sin(alpha);
    real cb = cos(beta),  sb = sin(beta);
    real cg = cos(gamma), sg = sin(gamma);

    // Rotation matrix rows are rot[i]:

    vec rot[3];
    rot[0] = vec( cg*cb*ca-sg*sa,  cg*cb*sa+sg*ca, -cg*sb);
    rot[1] = vec(-sg*cb*ca-cg*sa, -sg*cb*sa+cg*ca,  sg*sb);
    rot[2] = vec(sb*ca, sb*sa, cb);

    // Rotate the longitudinal and transverse vectors, then
    // reconstruct the new normal vector from them.

    vec l = k->get_longitudinal_unit_vector();
    l = rotate(l, rot);
    vec t = k->get_transverse_unit_vector();
    t = rotate(t, rot);
    vec n = l^t;
    k->set_orientation(l, t, n);
}

bool hdyn::integrate_unperturbed_motion(bool& reinitialize,
					bool force_time)   // default = false
{
    // Update unperturbed binary motion and check for continuation.
    // Return true iff binary is still unperturbed at end of step.
    //
    // Flag reinitialize is set true if a reinitialization of the system
    // is needed on return (e.g. after binary mass loss).
    //
    // Flag force_time will force the binary time to the current system
    // time if set.  Normally, binaries should be allowed to remain
    // asynchronous to preserve the proper phase.

    // NOTE: Integrating a binary to a time other than the scheduled
    //	     end of its timestep may result in its being resolved at
    //	     an undesirable part of its orbit, causing large errors.

    if (diag->unpert_function_id) {
	cerr << endl << ">> integrate_unperturbed_motion for "
	     << format_label() << " at time " << system_time << endl;
    }

    bool unpert = true;
    reinitialize = false;

    if (!kep) return false;

    hdyn *sister = get_binary_sister();
    time += unperturbed_timestep;

    if (slow) slow->inc_tau(unperturbed_timestep/get_kappa());

    if (time != system_time) {

	// Note: this can happen in case of changing binary orbit
        //	 because of stellar evolution.  Probably undesirable
	//	 to terminate the motion in this case, and also need
	//	 special care in determining new step.  TO BE FIXED...

	// It is also possible for this to occur when perturbed steps
	// become very short and we begin to lose significance in the
	// last bit.
	//						(Steve, 7/99)

	real dt = time - system_time;
	if (dt != 0) {

	    if (!force_time && diag->report_continue_unperturbed) {

		cerr << endl << "integrate_unperturbed_motion for "
		     << format_label() << endl;

		int p = cerr.precision(HIGH_PRECISION);
		cerr << "time and system_time are different:" << endl;
		PRC(time); PRL(system_time);

		// fprintf(stderr,
		//	"       time = '%Lx'\nsystem_time = '%Lx'\n",
		//	time, system_time);

		cerr << "Forcing time to system time..." << endl;
		cerr.precision(p);
	    }
	    time = system_time;
	}
    }

    // Store initial separation, for use below.

    real initial_sep_sq = kep->get_separation() * kep->get_separation();

    // An unperturbed step consists of transforming the kepler
    // structure to the new time, then updating the hdyn time, pos,
    // vel, acc, jerk, and perturbation, then determining if the
    // unperturbed motion is to continue.

    // Note:  Because the period is not exactly an integral number
    // of time steps, the mean and true anomalies, pos, vel, etc.,
    // will change slightly from one step to another.  This drift
    // may lead to restart problems, as the kepler structure
    // obtained using the new pos and vel may not be precisely the
    // same as the original kepler.  (SLWM, 3/98)

    // (Basically, the operation kepler --> (pos, vel) --> kepler is
    // not quite cyclic, leading to small phase errors that subsequently
    // grow in time.)

    update_dyn_from_kepler();	    // update pos, vel, and perturbation;
				    // note special treatment of slow motion

    // This update will introduce a tidal error (especially in unperturbed
    // multiples) if the unperturbed time step is not an integral multiple
    // of orbit periods.  The error may be acceptably small for binaries,
    // which are picked up with a more stringent unperturbed criterion (they
    // are continued with a relaxed criterion, but this is probably OK if
    // they are advanced to apocenter on the first step).  However, both
    // the inner and the outer binaries of an unperturbed multiple may be
    // quite strongly perturbed, and the tidal errors significant.

    // Fix:  always advance the inner and outer binaries by an *integral*
    //	     number of periods in the unperturbed multiple case (and be
    //	     very wary of partially unperturbed multiple motion...).
    //
    //							(Steve, 5/8/98)

    // Possible alternatives are to allow the tidal error and either
    // simply keep a running total of such errors, or to absorb the
    // error into the motion (e.g. by correcting the outer orbit in
    // some way) -- neither is currently implemented (but see below).

    bool save_unpert = fully_unperturbed;
    fully_unperturbed = false;

    if (unperturbed_timestep < kep->get_period()) {
        kc->partial_unpert_orbit++;
	// PRL(kc->partial_unpert_orbit);
    }
    else {
        kc->full_unpert_step++;
        kc->full_unpert_orbit += (int) (unperturbed_timestep
					   	/ kep->get_period());
	// PRC(kc->full_unpert_step);
	// PRL(kc->full_unpert_orbit);
    }

#if 0
    if (streq(format_label(), "xxx")) {
	cerr << "in integrate_unperturbed motion for " << format_label()
	     << endl;
	get_kepler()->print_all();
    }
#endif

    // Note that we use the same unperturbed criterion as in the main code,
    // and apply the same unperturbed_timestep function.

    bool is_u = is_unperturbed_and_approaching();

#if 0
    PRC(format_label()); PRC(kep->get_rel_pos()*kep->get_rel_vel()); PRL(is_u);
#endif

#if 0
    if (streq(format_label(), "xxx")) {
	PRL(posvel);
	PRL(is_u);
    }
#endif

    // int set_u = 0;
    // unsigned long set_u = 0;
    real set_u = 0;

    // cerr << "check set_unperturbed_timestep time: ";
    // PRC(cpu_time());

    if (is_u)
	set_u = set_unperturbed_timestep(!force_time);

    // PRL(cpu_time());

    // End of step.  Check to see if unperturbed motion can be extended.

    if (is_u && set_u) {

	// Is still unperturbed.  Not much more to do here now...

        if (diag->report_continue_unperturbed
	    || (diag->report_multiple && is_multiple(this)
		&& (diag->multiple_report_level > 0
		    || diag->unpert_report_level > 0))) {

	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << "\ncontinuing unperturbed motion for "
		 << format_label();
	    cerr << " and " << get_binary_sister()->format_label()
		 << endl
		 << "    at time " << time << " (next: " << set_u << " steps)"
		 << endl;

	    cerr.precision(p);
	    cerr << "    dt = " << timestep
		 << "  dt_unpert = " << unperturbed_timestep
		 << "  period = " << get_kepler()->get_period()
		 << endl
		 << "    perturbation = " << sqrt(perturbation_squared)
		 << "  (" << bt[init_binary_type] << ")"
		 << endl;

	    if (diag->unpert_report_level > 0)
		print_unperturbed_binary_params();

	}

    } else {

	// Time to end the unperturbed motion.

	unpert = false;

	bool verbose = diag->report_end_unperturbed
			    || (is_multiple(this) && diag->report_multiple);

	// PRL(verbose);

#if 1
	// The following lengthy workarounds are fixes for unperturbed
	// systems that find themselves in relatively wide triple orbits.
	//						    (Steve, 8/03)
	//
	// The problem is that, since the orientation of an unperturbed
	// binary remains fixed and its interaction with its neighbors is
	// computed in the center-of-mass approximation, a tidal energy
	// error accumulates during the unperturbed motion.
	//
	// The computation of the energy always resolves the binary into
	// its (static) components while it is unperturbed.  This has the
	// conveniant effect of ensuring that the computed energy remains
	// continuous as unperturbed motion starts.  However, since the
	// outer orbit is computed in the center-of-mass approximation,
	// while the energy is "exact" (the binary is frozen, but resolved),
	// the computed orbit differs from the true one due to the neglect
	// of the next (quadrupole) moment of the acceleration.  Integrated
	// around the outer  orbit, this error is responsible for a tidal
	// (quadrupole) error in the potential.
	//
	// This is a small but systematic effect, and can lead to a
	// significant net error if the outer orbit survives for many
	// periods.
	//
	// Rather than attempting to include the quadrupole (and higher?)
	// terms in the relative motion of the binary CM and its neighbors,
	// for now at least we choose to correct the tidal error each time
	// the unperturbed segment ends.
	//
	// Fix 1: randomize the orientation of any unperturbed binary
	//	  satisfying the criteria -- effective at eliminating the
	//	  systematic growth in the error, but overkill.
	// Fix 2: (better) absorb the tidal error directly into the CM and
	//	  nn relative motion.
	//
	// Is there a corresponding effect if the binary stays unperturbed
	// during the entire outer orbit?  No, since the tidal error in that
	// case will be strictly periodic.  The problem arises because the
	// partially unperturbed orbit is handled in two distinctly different
	// ways.

	hdyn *par = get_parent(), *pnn = par->get_nn();

	// Last clause of the following if() is a bit restrictive, but very
	// likely to be true, and for now we don't want to deal with all the
	// possible special configurations that might conceivably turn up if
	// we exclude it.
	//						     (Steve, 8/03)

	if (pnn
	    && binary_type != NOT_APPROACHING
	    && par->get_parent() == pnn->get_parent()) {

	    // Note that the code here is very similar to that used
	    // in startup_unperturbed_motion().

	    // Relative to root is too general if par and pnn have the
	    // same parent, but keep as is...

	    vec ppos = hdyn_something_relative_to_root(par,
						       &hdyn::get_pred_pos);
	    vec pvel = hdyn_something_relative_to_root(par,
						       &hdyn::get_pred_vel);
	    vec npos = hdyn_something_relative_to_root(pnn,
						       &hdyn::get_pred_pos);
	    vec nvel = hdyn_something_relative_to_root(pnn,
						       &hdyn::get_pred_vel);

	    real m12 = par->get_mass();
	    real m3 = pnn->get_mass();
	    real m123 = m12 + m3;
	    real mu123 = m12 * m3 / m123;
	    real R = abs(npos - ppos);
	    real Vsq = square(nvel - pvel);

	    real E = mu123*(0.5*Vsq - m123/R);

	    if (E < 0) {

		// We need kT for comparison.  This should be a relatively
		// rare calculation, so just do it the hard way.  (Could
		// perhaps save kT in the root dyn story...?)

		int ntop = 0;
		real kin = 0;
		for_all_daughters(hdyn, get_root(), bb) {
		    ntop++;
		    kin += bb->get_mass()*square(bb->get_pred_vel());
		}
		real kT = 2*kin/(3*ntop);

		if (E < -0.1*kT) {		// same 0.1 as in startup_...

		    bool randomize = false;	// try Fix 2 first

		    if (!randomize) {

			// Fix 2: Compute the tidal error and absorb it.
			//
			// Note:  Numerical experiments indicate that the
			// 	      change in the total energy is indeed well
			//	      described by de_tidal computed below.

			vec pos1 = ppos + get_pos();
			hdyn *sis = get_younger_sister();
			vec pos2 = ppos + sis->get_pos();
			real r13 = abs(npos - pos1);
			real r23 = abs(npos - pos2);
			real m1 = get_mass();
			real m2 = sis->get_mass();

			// Tidal error:

			real de_tidal = -m3 * (m1/r13 + m2/r23 - m12/R);

			if (verbose) {
			    int p = cerr.precision(HIGH_PRECISION);
			    cerr << endl
				 << "integrate_unperturbed_motion: "
				 << "absorbing tidal error for "
				 << parent->format_label()
				 << endl
				 << "    at time " << system_time
				 << ",  ";
			    cerr.precision(p);
			    PRC(E/kT); PRL(de_tidal);
#if 0
			    // Aside: Quadrupole approximation -m3*dphi seems to
			    // be a good approximation to the exact de_tidal:

			    real r12 = abs(get_pos()-sis->get_pos());
			    real dph = m1/r13 + m2/r23 - m12/R;

			    real costh = (pos2-pos1)*(npos-ppos)
						/(R*r12);
			    real mu12 = m1*m2/m12;
			    real dphi = 0.5*mu12*r12*r12*(3*costh*costh-1)
				    			 /pow(R,3);

			    PRI(4); PRC(m1); PRC(m2); PRC(m3); PRL(m12);
			    PRI(4); PRC(r12); PRC(r13); PRL(r23);
			    PRI(4); PRC(R); PRL(m1/r13 + m2/r23 - m12/R);
			    PRI(4); PRL(pos1);
			    PRI(4); PRL(pos2);
			    PRI(4); PRL(ppos);
			    PRI(4); PRL(npos);
			    PRI(4); PRL(pos2-pos1);
			    PRI(4); PRL(npos-ppos);
			    PRI(4); PRL((pos2-pos1)*(npos-ppos));
			    PRI(4); PRC(mu12); PRL(costh);
			    PRI(4); PRC(r12/R); PRL(r12*r12/pow(R,3));
			    PRI(4); PRC(dphi), PRL(-m3*dphi/de_tidal);
#endif
			}

			if (abs(de_tidal) > 0.25*mu123*Vsq) {

			    // Error de_tidal is suspiciously large.  Could
			    // correct exactly by rotating the binary, but
			    // for now just resort to randomizing it.

			    if (verbose)
				cerr << "    de_tidal too large: "
				     << "randomizing binary orientation"
				     << endl;

			    randomize = true;

			} else {

			    // Absorb the error.  Too complicated to adjust
			    // velocities in non-synchronous nodes, so start by 
			    // advancing par and pnn to the current time, if
			    // necessary.  (Do this before terminating the
			    // unperturbed motion, for consistency.)

			    par->synchronize_node();
			    pnn->synchronize_node();

			    // PRI(4); PRL(par->get_timestep());
			    // PRI(4); PRL(pnn->get_timestep());

			    // Note that we should really repeat the entire
			    // previous calculation, since all quantities will
			    // change slightly once we move from predicted to
			    // corrected quantities.  Ignore this detail and
			    // use the computed de_tidal, but recompute the
			    // velocities, Vsq, etc.

			    pvel = hdyn_something_relative_to_root(par,
							      &hdyn::get_vel);
			    nvel = hdyn_something_relative_to_root(pnn,
							      &hdyn::get_vel);
			    vec Vcm = (m12*pvel + m3*nvel) / m123;
			    vec V = nvel - pvel;
			    Vsq = square(V);

			    // Adjust velocities and correct jerks for the
			    // mutual interaction between par and pnn.  Want
			    // to change the relative velocity of par and pnn
			    // to increase 0.5*mu123*Vsq by -de_tidal, without
			    // changing their center-of-mass velocity.

			    real Vfac = sqrt(1 - de_tidal/(0.5*mu123*Vsq));
			    real facp = -Vfac*m3/m123;
			    real facn = Vfac*m12/m123;

			    if (verbose) {
				PRI(4); PRC(Vfac); PRC(facp); PRL(facn);
			    }

			    vec dpvel = Vcm + facp*V - pvel;
			    vec dnvel = Vcm + facn*V - nvel;

			    par->inc_vel(dpvel);
			    pnn->inc_vel(dnvel);

			    // *Neglect* jerk corrections for now...

			}
		    }

		    if (randomize) {		// value may have changed

			// Fix 1: Randomize the binary orientation.

			rotate_kepler(kep);
			update_dyn_from_kepler();

			if (verbose) {
			    int p = cerr.precision(HIGH_PRECISION);
			    cerr << endl
				 << "integrate_unperturbed_motion: "
				 << "applied random rotation to "
				 << parent->format_label()
				 << endl
				 << "    at time " << system_time
				 << ",  outer ";
			    cerr.precision(p);
			    PRL(E/kT);
			}
		    }
		}
	    }
	}
#endif

	// Correct any perturber lists containing the CM.  (This would be
	// repaired anyway at the next perturbee CM step, but that may be
	// too late...)

	// PRL(save_unpert);
	// PRL(binary_type);

	// Expand the binary into components before removing the kepler
	// in order to avoid complaints from expand_perturber_lists().

	if (save_unpert && !RESOLVE_UNPERTURBED_PERTURBERS)
	    expand_perturber_lists(get_root(), get_parent(), verbose);

	// NOTE: Ending the inner component of an unperturbed multiple
	//       should really end the outer component too.  Shouldn't
	//	 happen, given the current "unperturbed" criteria.

	if (verbose) {

	    // Pericenter reflection will be separating on termination.
	    // A transition from fully unperturbed to not approaching
	    // probably indicates internal binary evolution and should
	    // be flagged.

	    if (binary_type != NOT_APPROACHING	// (= end of peri reflection)
		|| diag->report_pericenter_reflection
		|| fully_unperturbed || save_unpert) {

		int p = cerr.precision(HIGH_PRECISION);
		cerr << "\nending unperturbed motion for "
		     << format_label();
		cerr << " and " << get_binary_sister()->format_label()
		     << " at time " << time << endl;
		cerr << bt[binary_type] << " (binary_type = " << binary_type;

		cerr.precision(p);
		cerr << ")  perturbation = " << sqrt(perturbation_squared);

		hdyn *pnn = get_parent()->get_nn();
		if (pnn && pnn->is_valid())
		    cerr << "  (" << pnn->format_label() << ")";

		cerr << endl;

#if 1
//		PRI(4); PRL(get_effective_block(time));
		PRI(4); PRL(sym_angle(kep->get_mean_anomaly()));
//		int pp = cerr.precision(20);
//		PRI(4); PRL(time);
//		PRI(4); PRL(system_time);
//		cerr.precision(pp);
#endif

		if (diag->unpert_report_level > 0
		    || diag->end_unpert_report_level > 0) {

		    print_unperturbed_binary_params();

		    if (binary_type != NOT_APPROACHING) {

			hdyn* p = get_parent();
			hdyn* pnn = p->get_nn();

			print_binary_from_dyn_pair(this, get_binary_sister(),
						   0, 0, true);
			cerr << endl;
			if (pnn) {
			    print_binary_from_dyn_pair(p, pnn, 0, 0, true);
			    cerr << endl;
			} else
			    cerr << "nn is NULL" << endl;

			if (pnn
			    && (diag->unpert_report_level > 1
				|| diag->end_unpert_report_level > 1)) {

			    pp3_minimal(get_top_level_node(), cerr);

			    // Attempt to estimate the work required to get
			    // past this encounter.

			    kepler outerkep;

			    outerkep.set_time(time);
			    outerkep.set_total_mass(p->get_mass());
			    outerkep.set_rel_pos(p->pos - pnn->pos);
			    outerkep.set_rel_vel(p->vel - pnn->vel);
			    outerkep.initialize_from_pos_and_vel();

			    real r_end = outerkep.get_separation();
			    if (perturbation_squared
				  > options->full_merge_tolerance)
				r_end
				  *= sqrt(perturbation_squared
					   / options->full_merge_tolerance);
			    if (r_end > outerkep.get_apastron())
				r_end = 0.9999*outerkep.get_apastron();
			    real transit_time
				= 2 * (outerkep.pred_advance_to_radius(r_end)
				       - (real)time);

			    real ave_step
				= timestep * sqrt(kep->get_periastron()
						  /kep->get_separation());

			    PRI(4); PRL(transit_time);
			    PRI(4); PRL(ave_step);

			    hdyn *pnode = find_perturber_node();
			    if (pnode) {
				cerr << "    estimate "
				     << 2 * pnode->n_perturbers
					  * (transit_time/ave_step) / 1e6
				     << " million force calculations to"
				     << " end of encounter"
				     << endl;
				cerr << "    perturber node: "
				     << pnode->format_label()
				     << endl;
				cerr << "    n_pert = "
				     << pnode->n_perturbers
				     << endl;
				PRI(4); print_nn(find_perturber_node(), 2);
			    }

			    PRI(4); PRC(is_u); PRL(set_u);
			}
		    }
		}
	    }

#if 0
	    if (streq(format_label(), "xxx"))
		get_kepler()->print_all();
#endif

	}

	bool update_dynamics[2] = {false, false};

	// Check for mass loss by binary evolution during the unperturbed
	// step.
	//
	//----------------------------------------------------------------
	//
	// Code now (4/99) somewhat redundant, as create_or_delete_binary
	// currently does *not* perform binary evolution, but simply
	// removes the dstar part of the CM node.
	//
	//----------------------------------------------------------------

	if (get_use_dstar() && has_dstar(this)) {

	    create_or_delete_binary(get_parent(),
				    update_dynamics);

	    //cerr << "integrate_unperturbed_motion: "
	    //     << "back from create_or_delete_binary" << endl << flush;

	    // If a merger occurred in create_or_delete_binary, both
	    // this and its coll particle have already been deleted!
	    // Detect this by checking if this is still a valid node.

	    if (!is_valid()) {
		reinitialize = true;
		return false;
	    }

	    // If the final orbit is radically different from the
	    // initial one, we have to recompute the time step.
	    // However, no need to do this if mass loss is forcing
	    // a complete reinitialization of the N-body system.

	    // update_dynamics[0] = true ==> full reinitialization needed
	    // update_dynamics[1] = true ==> new binary timestep needed

	    if (!update_dynamics[0] && update_dynamics[1]) {

		if (diag->report_start_unperturbed    ||
		    diag->report_continue_unperturbed ||
		    diag->report_end_unperturbed)
		    cerr << "\nCorrected timestep for change "
			 << "in binary elements.\n";

		update_dyn_from_kepler();   // may introduce a tidal error...
					    // -- should include in de_total

		set_first_timestep(0.5*kep->get_period());
		get_binary_sister()->timestep = timestep;

		get_parent()->mass = mass + get_binary_sister()->get_mass();
	    }

	    reinitialize = update_dynamics[0];
	    if (reinitialize)
		cerr << "integrate_unperturbed_motion: reinitialization "
		     << "forced by binary evolution" << endl
		     << "parent = " << get_parent()->format_label()
		     << " time = " << get_system_time() << endl;
	}

	// Finish up the termination of unperturbed motion (cf hdyn
	// constructor).

	delete kep;
	kep = NULL;
	get_binary_sister()->kep = NULL;
	unperturbed_timestep = -VERY_LARGE_NUMBER;
	get_binary_sister()->unperturbed_timestep = -VERY_LARGE_NUMBER;
        fully_unperturbed = false;
        get_binary_sister()->fully_unperturbed = false;

	// Check that the timestep is reasonable.

	if (!(update_dynamics[0] || update_dynamics[1])) {

	    // Separation may have changed if phase was adjusted during step...

	    real sep_sq = square(pos - get_binary_sister()->pos);

	    if (sep_sq < 0.75*initial_sep_sq) {

		// Time step is probably too long.  Reduce it appropriately.

		real rfac = sep_sq/initial_sep_sq;
		real target_timestep = timestep * pow(rfac, 0.75);
		while (timestep > target_timestep) timestep /= 2;

		if (diag->report_continue_unperturbed)
		    cerr << endl
			 << "integrate_unperturbed_motion: timestep for "
		         << format_label() << " reduced to " << timestep
			 << endl << "                              "
			 << "on restart at time " << time
			 << ", rfac = " << rfac
			 << endl;

		get_binary_sister()->timestep = timestep;
	    }
	}
    }

#if 0
    if (binary_type != NOT_APPROACHING) {
	PRC(timestep);
	PRL(get_binary_sister()->timestep);
	PRL(get_parent()->timestep);
    }
#endif

    return unpert;
}
