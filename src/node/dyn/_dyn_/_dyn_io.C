
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// _dyn_io: Starlab _dyn_ class I/O functions.
////
////          Define scan_dyn_story and print_dyn_story for
////          the _dyn_ class.
////
//// Options: none

#include "_dyn_.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize any static _dyn_ data here...
//
// ...

static bool read_xreal = false;

istream & _dyn_::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s, input_line), strcmp(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];

	sscanf(input_line, "%s%s", keyword, should_be_equal_sign);
	if (strcmp("=", should_be_equal_sign)) {
	    cerr << "Expected '=', but got '" << should_be_equal_sign
		 << endl;
	    exit(1);
	}

	// See xreal notes in dyn_io.C...

    	if (!strcmp("real_system_time", keyword)) {

	    read_xreal = true;
	    last_real = true;

	} else if (!strcmp("system_time", keyword)) {

	    // Check input format before reading.

	    if (!last_real) read_xreal = false;

	    if (read_xreal) {
		cerr << "_dyn_::scan_dyn_story: reading xreal time data"
		     << endl;

		set_system_time(get_xreal_from_input_line(input_line));
		cerr << "system_time = " << system_time << " ";
	    } else {
		cerr << "_dyn_::scan_dyn_story: reading real time data"
		     << endl; 
		real tmp;
		sscanf(input_line, "%*s%*s%lf", &tmp);
		real_system_time = system_time = tmp;
		cerr << "system_time = " << system_time << " ";
	    }

	} else {

	    last_real = false;

	    if (!strcmp("t", keyword)) {

		if (read_xreal)
		    time = get_xreal_from_input_line(input_line);
		else {
		    real tmp;
		    sscanf(input_line, "%*s%*s%lf", &tmp);
		    time = tmp;
		}

	    } else if (!strcmp("dt", keyword))
		sscanf(input_line, "%*s%*s%lf", &timestep);
	    else if (!strcmp("m", keyword))
		sscanf(input_line, "%*s%*s%lf", &mass);
	    else if (!strcmp("r", keyword))
		set_vector_from_input_line(pos, input_line);
	    else if (!strcmp("v", keyword))
		set_vector_from_input_line(vel, input_line);
	    else if (!strcmp("a", keyword))
		set_vector_from_input_line(acc, input_line);
	    else if (!strcmp("pot", keyword))
		sscanf(input_line, "%*s%*s%lf", &pot);
	    else if (!strcmp("R_eff", keyword))
		sscanf(input_line, "%*s%*s%lf", &radius);
	    else
		add_story_line(dyn_story, input_line);
	}
    }

    return s;
}

ostream & _dyn_::print_dyn_story(ostream & s,
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

    if (!short_output)
	put_real_number(s, "  dt =  ", timestep);

    put_real_number(s, "  m  =  ", mass);
    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);
    put_real_vector(s, "  a  =  ", acc);

    if (!short_output) {
	put_real_number(s, "  pot  =  ", pot);
	put_real_number(s, "  R_eff  =  ", radius);

	if (dyn_story)
	    put_story_contents(s, *dyn_story);
    }

    put_story_footer(s, DYNAMICS_ID);

    return s;
}

#else

main(int argc, char** argv)
{
    _dyn_  * b;
    check_help();

    while (b = get__dyn_(cin)) {
	cout << "TESTING put__dyn_:" << endl;
        put_node(cout, *b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	delete b;
    }
    cerr << "Normal exit" << endl;
}

#endif
