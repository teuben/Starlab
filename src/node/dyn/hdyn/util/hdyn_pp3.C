
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// hdyn_pp3:  Starlab hdyn-specific debugging tools.  Perform pp3_tree
////            on all nodes, descending from the root node.
////
//// Steve McMillan, 11/98.

// Externally visible functions:
//
//	void pp3_maximal
//	void pp3
//	void pp3_minimal
//	void pp3_tree

#include "hdyn.h"

#ifndef TOOLBOX

//========================================================================
//
// Helper functions:

#define MAX_INDENT	10
#define PP3_PRECISION	16

local bool pp3_check(hdyn* b)
{
    if (!b) return false;
    if (!b->is_valid()) return false;

    return true;
}

local void skip(int n, ostream & s)
{
    for (int i = 0; i < n; i++) s << " ";
}

local void print_node(char* label, hdyn * b, ostream & s)
{
    skip(MAX_INDENT, s); s << " \t" << label;

    if (b) {

	s << "(" << b << ")";
	if (b->is_valid())
	    s << "  " << b->format_label();
	else
	    s << "  invalid.";

    } else

	s << "(NULL)";

    s << endl;
}

local void print_id_and_time(hdyn * b, ostream & s, int level)
{
    char *string = b->format_label();
    int len = strlen(string);

    skip(2*level, s); s << string;

    skip(MAX_INDENT - len - 2*level, s);
    s << " \taddr: " << b << "  system_time: " << b->get_system_time() << endl;

#ifdef USE_XREAL
    skip(MAX_INDENT, s);
    s << " \txreal system_time = "; xprint(b->get_system_time(), s);
#endif

    skip(MAX_INDENT, s); s << " \tt: " << b->get_time()
				       << "  dt: " << b->get_timestep();
    if (b->get_kepler() || b->get_unperturbed_timestep() > 0)
	s << "  dt_u: " << b->get_unperturbed_timestep();
    else
	s << "  tp: " << b->get_t_pred();
    s << endl;

#ifdef USE_XREAL
    skip(MAX_INDENT, s);
    s << " \txreal time = "; xprint(b->get_time(), s);
#endif

}

local void print_pos_and_vel(hdyn * b, ostream & s, bool print_abs = false)
{
    skip(MAX_INDENT, s); s << " \tm: " << b->get_mass()
			   << "  |x|: " << abs(b->get_pos())
			   << endl;

    skip(MAX_INDENT, s); s << " \tx: " << b->get_pos()  << endl;
    if (!b->get_kepler()) {
	skip(MAX_INDENT, s); s << " \t-> " << b->get_nopred_pos()  << endl;
    }
    if (print_abs && b->is_low_level_node()) {
	skip(MAX_INDENT, s);
	s << " \tX: "
	  << hdyn_something_relative_to_root(b, &hdyn::get_pos)
	  << endl;
    }

    skip(MAX_INDENT, s); s << " \tv: " << b->get_vel()  << endl;
    if (!b->get_kepler()) {
	skip(MAX_INDENT, s); s << " \t-> " << b->get_nopred_vel()  << endl;
    }
    if (print_abs && b->is_low_level_node()) {
	skip(MAX_INDENT, s);
	s << " \tV: "
	  << hdyn_something_relative_to_root(b, &hdyn::get_vel)
	  << endl;
    }
}

//========================================================================
//
// Functions to print out information on a single node:

local void print_tree(hdyn * b, ostream & s, int level)
{
    if (!pp3_check(b)) return;

    print_node("parent:          ", b->get_parent(), s);
    print_node("elder_sister:    ", b->get_elder_sister(), s);
    print_node("younger_sister:  ", b->get_younger_sister(), s);
    print_node("oldest_daughter: ", b->get_oldest_daughter(), s);
}

local void print_minimal(hdyn * b, ostream & s, int level)
{
    if (!pp3_check(b)) return;

    // ID, address, times:

    print_id_and_time(b, s, level);

    // mass, pos, vel:

    print_pos_and_vel(b, s);

    // neighbor information:

    skip(MAX_INDENT, s); s << " \tnn: "; print_nn(b, 0, s);
    if (b->get_nn()) s << "  (" << sqrt(abs(b->get_d_nn_sq())) << ")";
    s << endl;
}

local void print_maximal(hdyn * b, ostream & s, int level)
{
    // Print out vital information in compact but readable form.

    if (!pp3_check(b)) return;

    //------------------------------------------------------------
    //
    // ID, address, times:

    print_id_and_time(b, s, level);

    //------------------------------------------------------------
    //
    // mass, pos, vel, acc, jerk:

    print_pos_and_vel(b, s, true);

    skip(MAX_INDENT, s); s << " \ta: " << b->get_acc()  << endl;
    skip(MAX_INDENT, s); s << " \tj: " << b->get_jerk() << endl;
    skip(MAX_INDENT, s); s << " \tk: " << 18*b->get_k_over_18() << endl;

    if (b->is_root()) {
	s << endl;
	return;
    }

    //------------------------------------------------------------
    //
    // nn, kepler, binary properties:

    skip(MAX_INDENT, s); s << " \tnn: "; print_nn(b, 0, s);
    if (b->get_nn()) s << "  (" << sqrt(abs(b->get_d_nn_sq())) << ")";
    if (b->is_low_level_node()) s << "  kep: " << b->get_kepler();
    s << endl;

    // Print binary properties only for elder binary component.

    if (b->is_low_level_node() && b->get_elder_sister() == NULL) {

	// Don't touch b's kepler structure, even if one exists.

	kepler k;
	initialize_kepler_from_dyn_pair(k, b, b->get_younger_sister(),
					true);		// minimal kepler

	int prec = s.precision(INT_PRECISION);
	skip(MAX_INDENT, s); s << " \tsma: " << k.get_semi_major_axis()
	    		       << "  ecc: "  << k.get_eccentricity()
			       << "  P: "    << k.get_period()
			       << endl;

	s.precision(prec);
    }

    //------------------------------------------------------------
    //
    // Perturbation:

    // Print perturber info only for elder binary component.

    if (!b->is_top_level_leaf()) {
	if (b->is_low_level_node() && b->get_elder_sister() == NULL) {
	    skip(MAX_INDENT, s);
	    if (b->get_perturbation_squared() >= 0
		&& b->get_perturbation_squared() < VERY_LARGE_NUMBER)
		s << " \tpert_sq: " << b->get_perturbation_squared() << endl;
	    else
		s << " \tpert_sq: unknown" << endl;
	}

	char pre[MAX_INDENT+3];
	for (int i = 0; i <= MAX_INDENT; i++) pre[i] = ' ';
	pre[MAX_INDENT+1] = '\t';
	pre[MAX_INDENT+2] = '\0';
#if 1
	if (!b->is_leaf())
	    b->print_perturber_list(s, pre);
	if (b->is_low_level_node() && b->get_elder_sister() == NULL)
	    b->find_print_perturber_list(s, pre);
#endif
    }

    //------------------------------------------------------------

    s << endl;
}

//========================================================================
//
// Externally visible functions:

void pp3_maximal(hdyn * b,				// maximal, recursive
		 ostream & s,	// default = cerr
		 int level)	// default = 0		// -1 ==> no recursion
{
    if (!pp3_check(b)) return;

    if (b->is_root() && level >= 0) {
	s << endl << "Static data:" << endl << endl;
	b->print_static(s);
    }

    int p = s.precision(PP3_PRECISION);
    print_maximal(b, s, level);
    s.precision(p);
    print_tree(b, s, level);

    if (level >= 0) {
	for_all_daughters(hdyn, b, daughter)
	    pp3_maximal(daughter, s, level + 1);
    }
}

void pp3(hdyn * b,					// standard, recursive
	 ostream & s,	// default = cerr
	 int level)	// default = 0			// -1 ==> no recursion
{
    if (!pp3_check(b)) return;

    if (b->is_root() && level >= 0) {
	s << endl << "Static data:" << endl << endl << flush;
	b->print_static(s);
    }

    int p = s.precision(PP3_PRECISION);
    print_maximal(b, s, level);
    s.precision(p);

    if (level >= 0) {
	for_all_daughters(hdyn, b, daughter)
	    pp3(daughter, s, level + 1);
    }
}

void pp3_minimal(hdyn * b,				// minimal, recursive
		 ostream & s,	// default = cerr
		 int level)	// default = 0		// -1 ==> no recursion
{
    if (!pp3_check(b)) return;

    int p = s.precision(PP3_PRECISION);
    print_minimal(b, s, level);
    s.precision(p);

    if (level >= 0) {
	for_all_daughters(hdyn, b, daughter)
	    pp3_minimal(daughter, s, level + 1);
    }
}

void pp3_tree(hdyn * b,					// data+links, recursive
	      ostream & s,	// default = cerr
	      int level)	// default = 0		// -1 ==> no recursion
{
    if (!pp3_check(b)) return;

    int p = s.precision(PP3_PRECISION);
    print_minimal(b, s, level);
    print_tree(b, s, level);
    s.precision(p);

    if (level >= 0) {
	for_all_daughters(hdyn, b, daughter)
	    pp3_tree(daughter, s, level + 1);
    }
}

void pp3_maximal(char *n,
		 ostream & s,	// default = cerr
		 int level)	// default = 0		// -1 ==> no recursion
{
    hdyn tmp;
    pp3_maximal((hdyn*)node_with_name(n, tmp.get_root()), s, level);
}

void pp3_minimal(char *n,
		 ostream & s,	// default = cerr
		 int level)	// default = 0		// -1 ==> no recursion
{
    hdyn tmp;
    pp3_minimal((hdyn*)node_with_name(n, tmp.get_root()), s, level);
}

void pp3(char *n,
	 ostream & s,		// default = cerr
	 int level)		// default = 0		// -1 ==> no recursion
{
    hdyn tmp;
    pp3((hdyn*)node_with_name(n, tmp.get_root()), s, level);
}

void pp3_tree(char *n,
	      ostream & s,	// default = cerr
	      int level)	// default = 0		// -1 ==> no recursion
{
    hdyn tmp;
    pp3_tree((hdyn*)node_with_name(n, tmp.get_root()), s, level);
}

#else

#include <sstream>

main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    hdyn *b;

    while (b = get_hdyn()) {

	//pp3_tree(b);

	ostringstream s;

	pp3(b->get_oldest_daughter(), s);
	cerr << s.str();

	delete b;
    }
}

#endif
