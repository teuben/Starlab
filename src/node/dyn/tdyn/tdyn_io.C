
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// tdyn_io:  Starlab pdyn and tdyn class I/O functions.
////
////           Define scan_star_story, scan_dyn_story and print_dyn_story
////           for the pdyn and tdyn classes.
////
//// Options: none

#include "tdyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize any static p/tdyn data here...
//
// ...

static bool read_xreal = false;

istream & pdyn::scan_star_story(istream & s, int level)
{
    // (Same function is inherited by tdyn.)

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
    return s;
}

bool tdyn::check_and_correct_node(bool verbose)		// default = false
{
    // Undo inherited dyn version; revert to node (just check masses).

    return node::check_and_correct_node(verbose);
}

istream & pdyn::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    bool last_real = false;
    bool reading_root = false;
    vec tmp;

    while (get_line(s, input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	// See xreal notes in dyn_io.C...

	// PRL(keyword);

	switch(keyword[0]) {
	    case 'a':

		// Acceleration:

		if (!strcmp("a", keyword)) {
		    set_vector_from_input_line(acc, input_line);
		    break;
		}
		goto other;
	
	    case 'L':

		// Luminosity:

		if (!strcmp("L", keyword)) {
		    luminosity = strtod(val, NULL);
		    break;
		}
		goto other;

	    case 'm':

		// Mass:

		if (!strcmp("m", keyword)) {
		    mass = strtod(val, NULL);
		    break;
		}
		goto other;
	
	    case 'r':

		// Position:

		if (!strcmp("r", keyword)) {
		    if (!reading_root)
			set_vector_from_input_line(pos, input_line);
		    else
			set_vector_from_input_line(tmp, input_line);
		    break;
		}

		// Real system time:

		if (!strcmp("real_system_time", keyword)) {

		    // "real_system_time" (a) gives the root-node time in
		    // real format and (b) serves as a flag that all other
		    // times are given in xreal form.

		    read_xreal = true;
		    last_real = true;
		    reading_root = true;
		    break;
		}
		goto other;

	    case 's':

		// System time:

		if (!strcmp("system_time", keyword)) {

		    // Check input format before reading.
		    // If we see "system_time" and haven't ever encountered
		    // "real_system_time", then all our times are plain
		    // "real"s.

		    if (!last_real) read_xreal = false;
		    reading_root = true;

		    if (read_xreal) {

			set_system_time(get_xreal_from_input_line(input_line));

		    } else {

			real_system_time = system_time = strtod(val, NULL);

		    }
		    break;
		}
		goto other;
	
	    case 'S':

		// Stellar type:

		if (!strcmp("S", keyword)) {

//		    char cptr[MAX_INPUT_LINE_LENGTH];
//		    sscanf(val,"%s",cptr);
//		    set_stellar_type(cptr);

		    stellar_type = strtol(val, NULL, 10);
		    break;
		}
		goto other;
	
	    case 't':

		if (!strcmp("tmpv", keyword)) {

		    // Read time, mass, pos, and vel as unformatted data.
		    // Input time will be real, independent of USE_XREAL.

		    // *** Must coordinate with hdyn_io.C. ***

		    real time = read_unformatted_real( s );
		    mass = read_unformatted_real( s );

		    if (!reading_root) {

			read_unformatted_vector( s, pos );
			read_unformatted_vector( s, vel );

		    } else {

			// Root pos and vel are used for center tracking.

			read_unformatted_vector( s, tmp );
			read_unformatted_vector( s, tmp );

		    }
		    break;
		}

		if (!strcmp("t64mpv32", keyword)) {
		    real time = read_unformatted_real( s );
		    mass = read_unformatted32_real( s );

		    if (!reading_root) {

			read_unformatted32_vector( s, pos );
			read_unformatted32_vector( s, vel );

		    } else {

			// Root pos and vel are used for center tracking.

			read_unformatted_vector( s, tmp );
			read_unformatted_vector( s, tmp );

		    }
		    break;
		}
		goto other;
	
	    case 'T':

		if (!strcmp("TL", keyword)) {

		    // Short output always uses floats for T and L.

		    temperature = read_unformatted32_real( s );
		    luminosity = read_unformatted32_real( s );
		    break;
		}

		// Temperature:

		if (!strcmp("T", keyword)) {
		    temperature = strtod(val, NULL);
		    break;
		}
		goto other;

	    case 'v':

		// Velocity:

		if (!strcmp("v", keyword)) {
		    if (!reading_root)
			set_vector_from_input_line(vel, input_line);
		    else
			set_vector_from_input_line(tmp, input_line);
		    break;
		}
		goto other;
	
	    default:
	      other:
		// else
		//    add_story_line(dyn_story, input_line);	// no stories!
		;
	}
    }
    return s;
}

istream & tdyn::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    bool last_real = false;
    bool reading_root = false;
    vec tmp;

    while (get_line(s, input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	// See xreal notes in dyn_io.C...

	// PRL(keyword);

	switch(keyword[0]) {
	    case 'a':

		// Acceleration:

		if (!strcmp("a", keyword)) {
		    set_vector_from_input_line(acc, input_line);
		    break;
		}
		goto other;
	
	    case 'b':

		// Alternate center tracking (root node only).  Save all
		// "center" quantities as stories for possible future use.

		if (reading_root) {

		    vec tmp;

		    if (!strcmp("bound_center_pos", keyword)) {
			set_vector_from_input_line(tmp, input_line);

			// Create a dyn story if none exists...

			if (!dyn_story)
			    dyn_story = mk_story_chapter(DYNAMICS_ID);

			putvq(dyn_story, "bound_center_pos", tmp);
			break;
		    }

		    if (!strcmp("bound_center_vel", keyword)) {
			set_vector_from_input_line(tmp, input_line);

			if (!dyn_story)
			    dyn_story = mk_story_chapter(DYNAMICS_ID);

			putvq(dyn_story, "bound_center_vel", tmp);
			break;
		    }
		}
		goto other;

	    case 'c':

		// Center tracking (root node only).  We will save any
		// "center" quantities to the root node pos and vel, and
		// save all such quantities as stories too.

		if (reading_root) {

		    if (!strcmp("center_pos", keyword)) {
			set_vector_from_input_line(pos, input_line);

			// Create a dyn story if none exists...

			if (!dyn_story)
			    dyn_story = mk_story_chapter(DYNAMICS_ID);

			putvq(dyn_story, "center_pos", pos);
			break;
		    }

		    if (!strcmp("center_vel", keyword)) {
			set_vector_from_input_line(vel, input_line);

			if (!dyn_story)
			    dyn_story = mk_story_chapter(DYNAMICS_ID);

			putvq(dyn_story, "center_vel", vel);
			break;
		    }

		    if (!strcmp("center_type", keyword)) {
			int type = strtol(val, NULL, 10);

			if (!dyn_story)
			    dyn_story = mk_story_chapter(DYNAMICS_ID);

			putiq(dyn_story, "center_type", type);
			break;
		    }
		}
		goto other;

	    case 'd':

		// Bookkeeping:

		if (!strcmp("defunct", keyword)) {
		    defunct = true;
		    break;
		}
		goto other;

	    case 'e':

		// Cluster escaper flag:

		if (!strcmp("esc", keyword)) {

		    // Use the prev pointer for temporary storage.  Careful!!
		    // Note: NULL means that esc is false.

		    // The esc flag is largely redundant, as all necessary
		    // information is contained in t_esc (below).

		    int esc = strtol(val, NULL, 10);
		    if (esc == 1) prev = (tdyn*) 42;
		    break;
		}
		goto other;
	
	    case 'j':

		// Jerk:

		if (!strcmp("j", keyword)) {
		    set_vector_from_input_line(jerk, input_line);
		    break;
		}
		goto other;
	
	    case 'k':

		// Kepler flag (1 = unperturbed, 2 = lightly perturbed):

		if (!strcmp("kep", keyword)) {
		    int i;
		    i = strtol(val, NULL, 10);
		    kep = (kepler*)i;
		    break;
		}
		goto other;

	    case 'L':

		// Luminosity:

		if (!strcmp("L", keyword)) {
		    luminosity = strtod(val, NULL);
		    break;
		}
		goto other;

	    case 'm':

		// Mass:

		if (!strcmp("m", keyword)) {
		    mass = strtod(val, NULL);
		    break;
		}
		goto other;
	
	    case 'r':

		// Position:

		if (!strcmp("r", keyword)) {
		    if (!reading_root)
			set_vector_from_input_line(pos, input_line);
		    else
			set_vector_from_input_line(tmp, input_line);
		    break;
		}

		// Real system time:

		if (!strcmp("real_system_time", keyword)) {

		    // "real_system_time" (a) gives the root-node time in
		    // real format and (b) serves as a flag that all other
		    // times are given in xreal form.

		    read_xreal = true;
		    last_real = true;
		    reading_root = true;
		    break;
		}
		goto other;

	    case 's':

		// System time:

		if (!strcmp("system_time", keyword)) {

		    // Check input format before reading.
		    // If we see "system_time" and haven't ever encountered
		    // "real_system_time", then all our times are plain
		    // "real"s.

		    if (!last_real) read_xreal = false;
		    reading_root = true;

		    if (read_xreal) {

			set_system_time(get_xreal_from_input_line(input_line));

		    } else {

			real_system_time = system_time = strtod(val, NULL);

		    }
		    break;
		}
		goto other;
	
	    case 'S':

		// Stellar type:

		if (!strcmp("S", keyword)) {

//		    char cptr[MAX_INPUT_LINE_LENGTH];
//		    sscanf(val,"%s",cptr);
//		    set_stellar_type(cptr);

		    stellar_type = strtol(val, NULL, 10);
		    break;
		}
		goto other;
	
	    case 't':

		// Time:

		if (!strcmp("t", keyword)) {
		    if (read_xreal)
			time = get_xreal_from_input_line(input_line);
		    else {
			time = strtod(val, NULL);
		    }
		    break;
		}

		if (!strcmp("tmpv", keyword)) {

		    // Read time, mass, pos, and vel as unformatted data.
		    // Input time will be real, independent of USE_XREAL.

		    // *** Must coordinate with hdyn_io.C. ***

		    time = read_unformatted_real( s );
		    mass = read_unformatted_real( s );

		    if (!reading_root) {

			read_unformatted_vector( s, pos );
			read_unformatted_vector( s, vel );

		    } else {

			// Root pos and vel are used for center tracking.

			read_unformatted_vector( s, tmp );
			read_unformatted_vector( s, tmp );

		    }
		    break;
		}

		if (!strcmp("t64mpv32", keyword)) {
		    time = read_unformatted_real( s );
		    mass = read_unformatted32_real( s );

		    if (!reading_root) {

			read_unformatted32_vector( s, pos );
			read_unformatted32_vector( s, vel );

		    } else {

			// Root pos and vel are used for center tracking.

			read_unformatted_vector( s, tmp );
			read_unformatted_vector( s, tmp );

		    }
		    break;
		}

		// Escape time.

		if (!strcmp("t_esc", keyword)) {

		    // Another kludge...  Careful again!!
		    // Return t_esc pointed to by t_next.
		    // MUST be sure to delete it on return.

		    real *t_esc = new real;
		    *t_esc = strtod(val, NULL);

		    next = (tdyn *)t_esc;
		}		
		goto other;
	
	    case 'T':

		if (!strcmp("TL", keyword)) {

		    // Short output always uses floats for T and L.

		    temperature = read_unformatted32_real( s );
		    luminosity = read_unformatted32_real( s );
		    break;
		}

		// Temperature:

		if (!strcmp("T", keyword)) {
		    temperature = strtod(val, NULL);
		    break;
		}
		goto other;

	    case 'v':

		// Velocity:

		if (!strcmp("v", keyword)) {
		    if (!reading_root)
			set_vector_from_input_line(vel, input_line);
		    else
			set_vector_from_input_line(tmp, input_line);
		    break;
		}
		goto other;
	
	    default:
	      other:
		// else
		//    add_story_line(dyn_story, input_line);	// no stories
		;
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

ostream & _pdyn_::print_dyn_story(ostream & s,
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

// Not needed (should be inherited?):

//  ostream & pdyn::print_dyn_story(ostream & s,
//  				bool print_xreal,	// default = true
//  				int short_output)	// default = 0, ignored
//  {
//      _pdyn_::print_dyn_story(s, print_xreal, short_output);
//      return s;
//  }

ostream & tdyn::print_dyn_story(ostream & s,
				bool print_xreal,	// default = true
				int short_output)	// default = 0, ignored
{
    // Modifications by Steve (5/01) to streamline output.

    // Ordinarily we want to print time before dyn stuff, but it is
    // necessary to print system time first for the root node...

    if (parent) print_local_time(time, s, print_xreal, short_output);

    // Use pdyn::print_dyn_story() to print the pdyn stuff...

    _pdyn_::print_dyn_story(s, print_xreal, short_output);

    if (!parent) print_local_time(time, s, print_xreal, short_output);

    return s;
}

#else

main(int argc, char** argv)
{
    tdyn *b;
    check_help();

#if 1

    while (b = get_tdyn()) {
	cout << "TESTING put_tdyn:" << endl;
        put_node(b, cout, true, 1);
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
