
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// tdyn_io:  Starlab pdyn and tdyn class I/O functions.
////
////           Define scan_dyn_story and print_dyn_story for
////           the pdyn and tdyn classes.
////
//// Options: none

#include "tdyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize any static tdyn data here...
//
// ...

static bool read_xreal = false;

istream & tdyn::scan_star_story(istream & s, int level)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    // In this case, simply read the Star info into the dyn variables.
    // This code replicates part of scan_dyn_story() below, and allows
    // us to read a snapshot with the stellar information encoded in
    // a (Star...)Star clause instead of in Dyn.

    while (get_line(s, input_line), !matchbracket(END_STAR, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (!strcmp("S", keyword)) {

//	    char cptr[MAX_INPUT_LINE_LENGTH];
//	    sscanf(val,"%s",cptr);
//	    set_stellar_type(cptr);

	    stellar_type = strtol(val, NULL, 10);

	} else if (!strcmp("T", keyword))
	    temperature = strtod(val, NULL);
	else if (!strcmp("L", keyword))
	    luminosity = strtod(val, NULL);

	// Ignore everything else -- no stories!
    }
}

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
	    const char *val = getequals(input_line, keyword);

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

		    real_system_time = system_time = strtod(val, NULL);

		}

	    } else {

		last_real = false;

		// Dynamical data:

		if (!strcmp("t", keyword)) {

		    if (read_xreal)
			time = get_xreal_from_input_line(input_line);
		    else {
			time = strtod(val, NULL);
		    }

		} else if (!strcmp("m", keyword))
		    mass = strtod(val, NULL);
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
		    i = strtol(val, NULL, 10);
		    kep = (kepler*)1;		// just use as a flag for now
		}

		// Stellar data:

		else if (!strcmp("S", keyword)) {

//		    char cptr[MAX_INPUT_LINE_LENGTH];
//		    sscanf(val,"%s",cptr);
//		    set_stellar_type(cptr);

		    stellar_type = strtol(val, NULL, 10);

		} else if (!strcmp("T", keyword))
		    temperature = strtod(val, NULL);
		else if (!strcmp("L", keyword))
		    luminosity = strtod(val, NULL);
		else if(streq("TL =", input_line)) {

		    // Short output always uses floats for T and L.

		    temperature = read_unformatted32_real( s );
		    luminosity = read_unformatted32_real( s );
		}

		// Bookkeeping:

		else if (!strcmp("defunct", keyword))
		    defunct = true;

		// else
		//    add_story_line(dyn_story, input_line);	// no stories

	    }
	}
    }
    return s;
}

local void print_local_time(xreal time,
			    ostream & s,
			    bool print_xreal,
			    int short_output)
{
    adjust_starlab_precision(HIGH_PRECISION);

    if (print_xreal)
	put_real_number(s, "  t  =  ", time);		// OK for real or xreal
    else
	put_real_number(s, "  t  =  ", (real)time);

    adjust_starlab_precision(-1);
}

ostream & pdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0, ignored
{
    // Use dyn::print_dyn_story() to print the dyn stuff...

    dyn::print_dyn_story(s, print_xreal, short_output);

//    if (stellar_type) put_string(s,      "  S  =  ", stellar_type);

    put_integer(s, "  S  =  ", stellar_type);
    put_real_number(s, "  T  =  ", temperature);
    put_real_number(s, "  L  =  ", luminosity);

    return s;
}

ostream & tdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0, ignored
{
    // Modifications by Steve (5/01) to streamline output.

    // Ordinarily we want to print time before dyn stuff, but it is
    // necessary to print system time first for the root node...

    if (parent) print_local_time(time, s, print_xreal, short_output);

    // Use pdyn::print_dyn_story() to print the pdyn stuff...

    pdyn::print_dyn_story(s, print_xreal, short_output);

    if (!parent) print_local_time(time, s, print_xreal, short_output);

    return s;
}

#else

main(int argc, char** argv)
{
    tdyn *b;
    check_help();

#if 1

    while (b = get_tdyn(cin)) {
	cout << "TESTING put_tdyn:" << endl;
        put_node(cout, *b, true, 1);
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

	b = get_tdyn(s, NULL, NULL, false);

	cerr << "pause..." << flush;
	cin >> tmp;
    } else
	cerr << "no file xxx" << endl;

#endif
}

#endif
