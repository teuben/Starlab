
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//
//  _dyn_slow.C: functions related to slow binary motion.  Manipulations
//		 are defined for the _dyn_ class, but are presently used
//		 only in kira, for hdyn objects.
//
//.............................................................................
//    version 1:  Jul 1999   Steve McMillan
//    version 2:
//.............................................................................
//
// Externally visible functions:
//
//	void _dyn_::create_slow
//	void _dyn_::delete_slow
//	void _dyn_::extend_slow
//
//	void _dyn_::check_slow_perturbed
//	slow_perturbed *_dyn_::find_slow_perturbed
//	slow_perturbed *_dyn_::add_slow_perturbed
//	void _dyn_::remove_slow_perturbed
//	int _dyn_::count_slow_perturbed
//
// Note that there is a great deal of overlap in coding between the functions
// in each group.

#include "_dyn_.h"

void _dyn_::create_slow(int k)		// default = 1
{
    if (!is_low_level_node()) {
	warning("create_slow: not a low-level node");
	return;
    }

    slow = new slow_binary(k);

    slow->set_t_init(get_time());
    slow->set_t_apo(get_time());

    // New dtau is the old timestep.

    slow->set_dtau(timestep);
    timestep *= k;

    // Correct the non-2-body portion of the acceleration.
    // (Assume binaries, not multiples, for now...)

    _dyn_ *s = get_younger_sister();

    if (!s) {

	warning("create_slow: no younger sister");

    } else {

	real m2 = get_binary_sister()->mass;
	vector sep = pos * (1 + mass/m2);
	real r2 = sep*sep;
	vector a2 = -m2*sep / (r2*sqrt(r2));

	acc = a2 + (acc - a2) * k;
	old_acc = acc;			// OK because we do this between steps

	// Update the younger sister.

	s->slow = slow;
	real factor = -s->mass / mass;
	s->acc = factor * acc;
	s->old_acc = s->acc;

    }
}

void _dyn_::delete_slow()
{
    if (!slow) {

	cerr << "warning: delete_slow: " << format_label()
	     << " is not a slow binary" << endl;
	return;

    } else {

	int k = get_kappa();

	timestep /= k;

	// Correct the non-2-body portion of the acceleration.
	// (Assume binaries, not multiples, for now...)

	// Assume that this function is invoked from the elder sister
	// of a binary.

	_dyn_ *s = get_younger_sister();

	if (!s) {

	    warning("delete_slow: no younger sister...");

	} else {

	    real m2 = s->mass;
	    vector sep = pos * (1 + mass/m2);
	    real r2 = sep*sep;
	    vector a2 = -m2*sep / (r2*sqrt(r2));

	    acc = a2 + (acc - a2) / k;
	    old_acc = acc;		// OK because we do this between steps

	    // Update the sister.

	    real factor = -s->mass / mass;
	    s->acc = factor * acc;
	    s->old_acc = s->acc;

	}

	delete slow;

	slow = NULL;
	if (s) s->slow = NULL;
    }
}

void _dyn_::extend_slow(int k)		// no default
{
    if (!slow) {

	cerr << "warning: extend_slow: " << format_label()
	     << " is not a slow binary" << endl;
	return;

    } else {

	real k_fac = ((real)k) / get_kappa();
	slow->set_kappa(k);

	slow->set_t_init(get_time());
	slow->set_t_apo(get_time());

	// Must reset tau to zero; dtau is unchanged...

	slow->set_tau(0);
	slow->init_tau_pred();

	timestep *= k_fac;

	// Correct the non-2-body portion of the acceleration.
	// (Assume binaries, not multiples, for now...)

	_dyn_ *s = get_younger_sister();

	if (!s) {

	    warning("extend_slow: no younger sister");

	} else {

	    real m2 = get_binary_sister()->mass;
	    vector sep = pos * (1 + mass/m2);
	    real r2 = sep*sep;
	    vector a2 = -m2*sep / (r2*sqrt(r2));

	    acc = a2 + (acc - a2) * k_fac;
	    old_acc = acc;		// OK because we do this between steps

	    // Update the younger sister.

	    real factor = -s->mass / mass;
	    s->acc = factor * acc;
	    s->old_acc = s->acc;

	}
    }
}

// Member functions relating to the slow_perturbed lists.

bool is_valid_slow(_dyn_ *pert_node)
{
    if (!pert_node
	|| !pert_node->is_valid()
	|| !pert_node->get_oldest_daughter()
	|| !pert_node->get_oldest_daughter()->get_slow())
	return false;
    else
	return true;
}

typedef slow_perturbed* sp_ptr;

local void delete_slow_perturbed(_dyn_ *b,
				 sp_ptr &s,
				 slow_perturbed *prev,
				 sp_ptr &sp,
				 bool verbose)
{
    // Remove entry s from the slow_perturbed list of b.

    slow_perturbed *next = s->get_next();

    if (prev)
	prev->set_next(next);
    else
	sp = next;

    if (verbose) {
	cerr << "deleted " << s->get_node()->format_label();
	cerr << " from slow_perturbed list of " << b->format_label()
	     << " at time " << b->get_system_time()
	     << endl;
    }

    s->set_next(NULL);		// so we don't delete the entire chain
    delete s;
    s = next;
}

void _dyn_::check_slow_perturbed(bool verbose)		// default = false
{
    slow_perturbed *s = sp, *prev = NULL;
    while (s) {
	_dyn_ *pert_node = s->get_node();

	// Apply basic checks:

	if (!is_valid_slow(pert_node))

	    delete_slow_perturbed(this, s, prev, sp, verbose);

	else {

	    prev = s;
	    s = s->get_next();

	}
    }
}

slow_perturbed* _dyn_::find_slow_perturbed(_dyn_ *n,	  // checks, too...
					   bool verbose)  // default = false
{
    slow_perturbed *s = sp, *prev = NULL;
    while (s) {
	_dyn_ *pert_node = s->get_node();

	// Apply basic checks:

	if (!is_valid_slow(pert_node))

	    delete_slow_perturbed(this, s, prev, sp, verbose);

	else {

	    if (pert_node == n) return s;

	    prev = s;
	    s = s->get_next();

	}
    }
    return NULL;
}

slow_perturbed* _dyn_::add_slow_perturbed(_dyn_ *n,
					  bool verbose)	  // default = false
{
    if (!is_valid_slow(n)) return NULL;

    slow_perturbed *s = sp, *prev = NULL;
    while (s) {
	if (s->get_node() == n) return s;
	prev = s;
	s = s->get_next();
    }

    s = new slow_perturbed();
    s->set_node(n);
    s->set_kappa(n->get_oldest_daughter()->get_kappa());

    if (prev)
	prev->set_next(s);
    else
	sp = s;

    if (verbose) {
	cerr << "added " << n->format_label();
	cerr << " to slow_perturbed list of " << format_label()
	     << " at time " << get_system_time()
	     << endl;
    }

    return s;
}
	    
void _dyn_::remove_slow_perturbed(_dyn_ *n,
				  bool verbose)		// default = false
{
    slow_perturbed *s = sp, *prev = NULL;
    while (s) {

	if (s->get_node() == n) {
	    delete_slow_perturbed(this, s, prev, sp, verbose);
	    return;
	}

	prev = s;
	s = s->get_next();
    }
}

bool _dyn_::copy_slow_perturbed(_dyn_ *to,
				bool overwrite,		// default = false
				bool verbose)		// default = false
{
    // Copy the slow_perturbed list from this node to another.
    // Optionally delete any existing slow_perturbed list before overwriting.

    if (!sp) return false;

    if (to->get_sp()) {
	if (overwrite) {
	    if (verbose) cerr << "copy_slow_perturbed: "
			      << "deleting existing slow_perturbed list"
			      << endl;
	    delete to->get_sp();
	} else {
	    if (verbose) cerr << "copy_slow_perturbed: "
			      << "adding to existing slow_perturbed list"
			      << endl;
	}
    }

    slow_perturbed *s = sp;
    while (s) {
	to->add_slow_perturbed(s->get_node(), verbose);
	s = s->get_next();
    }

    return true;
}

void _dyn_::dump_slow_perturbed(char *string)		// default = ""
{
    // Dump all slow_perturber data.

    if (!sp) return;

    cerr << string
	 << "dump slow_perturbed list for "
	 << format_label() << " at time " << get_system_time() << ":"
	 << endl;

    slow_perturbed *s = sp;
    int n = 0;
    while (s) {
	n++;

	cerr << string;	PRC(n); PRL(s);
	cerr << string;	PRI(4); PRC(s->get_node());
			PRL(s->get_node()->format_label());
	cerr << string;	PRI(4); PRL(s->get_kappa());
	cerr << string;	PRI(4); PRL(s->get_next());

	s = s->get_next();
    }
}

void _dyn_::print_slow_perturbed()	// print the slow_perturbed node list
{
    if (!sp) return;

    cerr << "slow_perturbed list of " << format_label() << ":"
	 << endl;

    slow_perturbed *s = sp;
    while (s) {
	cerr << " " << s->get_node()->format_label();
	s = s->get_next();
    }
    cerr << endl;
}

int _dyn_::count_slow_perturbed()
{
    slow_perturbed *s = sp;
    int n = 0;
    while (s) {
	n++;
	s = s->get_next();
    }
    return n;
}
	    
