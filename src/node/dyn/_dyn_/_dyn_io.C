
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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

void _dyn_::print_static(ostream& s)		// default = cerr
{
    dyn::print_static(s);
}

void _dyn_::null_pointers()
{
    // Clear all pointers (don't touch what they point to)
    // -- for use in cleaning up temporary nodes...  Careful!!

    slow = NULL;
    sp =NULL;

    dyn::null_pointers();
}

istream & _dyn_::scan_dyn_story(istream & s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s, input_line), strcmp(END_DYNAMICS, input_line)) {

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
		cerr << "_dyn_::scan_dyn_story: reading xreal time data"
		     << endl;

		set_system_time(get_xreal_from_input_line(input_line));
		cerr << "system_time = " << system_time << " ";
	    } else {
		cerr << "_dyn_::scan_dyn_story: reading real time data"
		     << endl;
		real_system_time = system_time = strtod(val, NULL);
		cerr << "system_time = " << system_time << " ";
	    }

	} else {

	    last_real = false;

	    if (!strcmp("t", keyword)) {

		if (read_xreal)
		    time = get_xreal_from_input_line(input_line);
		else {
		    time = strtod(val, NULL);
		}

	    } else if (!strcmp("dt", keyword))
		timestep = strtod(val, NULL);
	    else if (!strcmp("m", keyword))
		mass = strtod(val, NULL);
	    else if (!strcmp("r", keyword))
		set_vector_from_input_line(pos, input_line);
	    else if (!strcmp("v", keyword))
		set_vector_from_input_line(vel, input_line);
	    else if (!strcmp("a", keyword))
		set_vector_from_input_line(acc, input_line);
	    else if (!strcmp("pot", keyword))
		pot = strtod(val, NULL);
	    else if (!strcmp("R_eff", keyword))
		radius = strtod(val, NULL);
	    else
		add_story_line(dyn_story, input_line);
	}
    }

    return s;
}

local void print_local_time(xreal time,
			    xreal system_time,
			    ostream & s,
			    bool print_xreal,
			    int short_output)
{
    // Always print high-precision real system_time for t in the
    // short_output case.  (In general, the output precision is
    // controlled with the STARLAB_PRECISION environment variable.)

    if (short_output) {

	adjust_starlab_precision(HIGH_PRECISION);
	put_real_number(s, "  t  =  ", (real)system_time);
	if (short_output) adjust_starlab_precision(-1);

    } else {

	if (print_xreal)
	    put_real_number(s, "  t  =  ", time);	// OK for real or xreal
	else
	    put_real_number(s, "  t  =  ", (real)time);

    }
}

ostream & _dyn_::print_dyn_story(ostream & s,
				 bool print_xreal,	// default = true
				 int short_output)	// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    // Ordinarily we want to print time before dyn stuff, but it is
    // necessary to print system time first for the root node...

    if (parent) print_local_time(time, system_time,
				 s, print_xreal, short_output);

    // Awkward (dyn output prints pos), but...

    vector tmp_pos, tmp_vel;
    if (short_output > 1) {
	tmp_pos = pos;
	tmp_vel = vel;
	pos = pred_pos;
	vel = pred_vel;
    }

    // Use dyn::print_dyn_story() to print dyn stuff...

    dyn::print_dyn_story(s, print_xreal, short_output);

    if (short_output > 1) {
	pos = tmp_pos;
	vel = tmp_vel;
    }

    if (!parent) print_local_time(time, system_time,
				  s, print_xreal, short_output);

    if (!short_output) {
	put_real_vector(s, "  a  =  ", acc);
	put_real_number(s, "  pot  =  ", pot);
	put_real_number(s, "  dt =  ", timestep);
	put_real_number(s, "  R_eff  =  ", radius);
    }

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
