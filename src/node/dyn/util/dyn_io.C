
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// dyn_io:  test Starlab dyn class I/O functions.
////
////          Define scan_dyn_story and print_dyn_story for
////          the dyn class.
////
//// Options: none

#include "dyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize all static dyn data here...

xreal dyn::system_time          = 0.0;
real dyn::real_system_time      = 0.0;
bool dyn::use_sstar	        = false;

unsigned int  dyn::external_field	= 0;

int  dyn::tidal_type	= 0;
real dyn::omega	= 0;
real dyn::omega_sq	= 0;
real dyn::alpha1	= 0;
real dyn::alpha3	= 0;
vector dyn::tidal_center = vector(0,0,0);

real dyn::p_mass = 0;
real dyn::p_scale_sq = 0;
vector dyn::p_center = vector(0,0,0);

real dyn::pl_coeff = 0;
real dyn::pl_scale_sq = 0;
real dyn::pl_exponent = 0;
vector dyn::pl_center = vector(0,0,0);

static bool read_xreal = false;

istream & dyn::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s,input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

    	if (!strcmp("real_system_time", keyword)) {

	    read_xreal = true;
	    last_real = true;

	    // We don't know anything about parent nodes yet, so it is
	    // not easy to know if we are the root node.  Rule: if we
	    // find real_system_time, assume that we should read an
	    // xreal as system_time.  Otherwise, read it as real.
	    // The tortuous logic is to keep the determination of
	    // which choice we should make completely local.
	    //
	    // Unfortunately, this logic must be duplicated in all
	    // other *dyn::scan_dyn_story functions (see _dyn_io.C,
	    // hdyn_io.C, sdyn3_io.C)...

	} else if (!strcmp("system_time", keyword)) {

	    // Check input format before reading.

	    if (!last_real) read_xreal = false;

	    if (read_xreal) {

		//cerr << "dyn::scan_dyn_story: input "
		//     << "time data type is xreal"
		//     << endl;

		// The following should set real_system_time too...

		set_system_time(get_xreal_from_input_line(input_line));

	    } else {

		//cerr << "dyn::scan_dyn_story: input "
		//     << "time data type is real"
		//     << endl;

		real_system_time = system_time = strtod(val, NULL);
	    }
	} else {

	    last_real = false;

	    if (!strcmp("m", keyword))
		mass = strtod(val, NULL);
	    else if (!strcmp("r", keyword))
		set_vector_from_input_line(pos, input_line);
	    else if (!strcmp("v", keyword))
		set_vector_from_input_line(vel, input_line);
	    else
		add_story_line(dyn_story, input_line);
	}
    }

    return s;
}

ostream& dyn::print_dyn_story(ostream& s,
			      bool print_xreal,		// default = true
			      int short_output)	// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    // Print system time first (root node only).

    if (!parent) {

#ifdef USE_XREAL
	if (print_xreal) {

	    // Note (Steve 5/00): system_time is now xreal and hard to read.
	    // For convenience, also print out a real version of the time.
	    // By printing out real_system_time first, we set a flag that
	    // allows scan_dyn_story to read xreal, rather than real, input.

	    put_real_number(s, "  real_system_time  =  ", (real)system_time);
	    put_real_number(s, "  system_time  =  ", system_time);

	} else

	    put_real_number(s, "  system_time  =  ", (real)system_time);
#else

	put_real_number(s, "  system_time  =  ", system_time);

#endif
    }

    // Mass is now printed by node::print_dyn_story().

    node::print_dyn_story(s, print_xreal, short_output);

    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);

    return s;
}

#else
main(int argc, char** argv)
{
    dyn  * b;
    check_help();

    while (b = (dyn *) get_node(cin, new_dyn)) {
	cout << "TESTING put_dyn:" << endl;
        put_node(cout,*b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	delete b;
    }
    cerr << "Normal exit\n";
}
#endif
