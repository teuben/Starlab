
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// tdyn_io:  Starlab tdyn class I/O functions.
////
////           Define scan_dyn_story and print_dyn_story for
////           the tdyn class.
////
//// Options: none

#include "tdyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize any static tdyn data here...
//
// ...

static bool read_xreal = false;

istream & tdyn::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while (get_line(s, input_line), !matchbracket(END_DYNAMICS, input_line)) {

	// Special case:

	if (streq("tmpv =", input_line)) {

	    // Read time, mass, pos, and vel as unformatted data.
	    // Input time will be real, independent of USE_XREAL.

	    // *** Must coordinate with hdyn_io.C. ***

	    time = read_unformatted_real( s );

	    mass = read_unformatted_real( s );
	    read_unformatted_vector( s, pos );
	    read_unformatted_vector( s, vel );

	} else if(streq("t64mpv32 =", input_line)) {

	    time = read_unformatted_real( s );
	    mass = read_unformatted32_real( s );
	    read_unformatted32_vector( s, pos );
	    read_unformatted32_vector( s, vel );

	} else {

	    // Usual formatted input.

	    real last_real = false;

	    char keyword[MAX_INPUT_LINE_LENGTH];
	    char *val = getequals(input_line, keyword);

	    // See xreal notes in dyn_io.C...

	    if (!strcmp("real_system_time", keyword)) {

		read_xreal = true;
		last_real = true;

	    } else if (!strcmp("system_time", keyword)) {

		// Check input format before reading.

		if (!last_real) read_xreal = false;

		if (read_xreal) {

		    set_system_time(get_xreal_from_input_line(input_line));

		} else {

		    real_system_time = system_time = strtod(val, &val);

		}

	    } else {

		last_real = false;

		if (!strcmp("t", keyword)) {

		    if (read_xreal)
			time = get_xreal_from_input_line(input_line);
		    else {
			time = strtod(val, &val);
		    }

		} else if (!strcmp("m", keyword))
		    mass = strtod(val, &val);
		else if (!strcmp("r", keyword))
		    set_vector_from_input_line(pos, input_line);
		else if (!strcmp("v", keyword))
		    set_vector_from_input_line(vel, input_line);
		else if (!strcmp("a", keyword))
		    set_vector_from_input_line(acc, input_line);
		else if (!strcmp("j", keyword))
		    set_vector_from_input_line(jerk, input_line);
		else if (!strcmp("kep", keyword)) {
		    int i;
		    i = strtol(val, &val, 10);
		    kep = (kepler*)1;		// just use as a flag for now
		}
		else if (!strcmp("defunct", keyword))
		    defunct = true;
		else
		    add_story_line(dyn_story, input_line);
	    }
	}
    }
    return s;
}

ostream & tdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0
{
    put_story_header(s, DYNAMICS_ID);

    if (!parent) {

	// See xreal notes in dyn_io.C...

#ifdef USE_XREAL
	if (print_xreal) {

	    put_real_number(s, "  real_system_time  =  ", (real)system_time);
	    put_real_number(s, "  system_time  =  ", system_time);

	} else

	    put_real_number(s, "  system_time  =  ", (real)system_time);
#else

	put_real_number(s, "  system_time  =  ", system_time);

#endif
    }

    if (print_xreal)
	put_real_number(s, "  t  =  ", time);		// OK for real or xreal
    else
	put_real_number(s, "  t  =  ", (real)time);

    put_real_number(s, "  m  =  ", mass);
    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);

    // Not needed:

    // put_real_vector(s, "  a  =  ", acc);
    // put_real_vector(s, "  j  =  ", acc);

    put_story_footer(s, DYNAMICS_ID);

    return s;
}

#else

main(int argc, char** argv)
{
    tdyn *b;
    check_help();

#if 0

    while (b = get_tdyn(cin)) {
	cout << "TESTING put_tdyn:" << endl;
        put_node(cout, *b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	delete b;
    }

#else

    ifstream s("xxx");
    if (s) {
	cerr << "pause..." << flush;
	char tmp;
	cin >> tmp;

	b = get_tdyn(s, NULL, NULL, false);//, NULL, NULL, false);

	cerr << "pause..." << flush;
	cin >> tmp;
    } else
	cerr << "no file xxx" << endl;

#endif
}

#endif
