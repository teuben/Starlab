
//// sdyn3_io:  Starlab sdyn3-specific I/O functions
////
//// Options:  none

#include "sdyn3.h"
#include "util_io.h"

#ifndef TOOLBOX

static bool read_xreal = false;

istream & sdyn3::scan_dyn_story(istream & s)
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

	    if (read_xreal)
		system_time = get_xreal_from_input_line(input_line);
	    else
		sscanf(input_line, "%*s%*s%lf", &system_time);

	} else {

	    last_real = false;

	    if (!strcmp("t", keyword)) {

		if (read_xreal)
		    time = get_xreal_from_input_line(input_line);
		else
		    sscanf(input_line, "%*s%*s%lf", &time);

	    } else if (!strcmp("m", keyword))
		sscanf(input_line, "%*s%*s%lf", &mass);
	    else if (!strcmp("r", keyword))
		set_vector_from_input_line(pos, input_line);
	    else if (!strcmp("v", keyword))
		set_vector_from_input_line(vel, input_line);
	    else if (!strcmp("dt", keyword))
		sscanf(input_line, "%*s%*s%lf", &timestep);
	    else if (!strcmp("a", keyword))
		set_vector_from_input_line(acc, input_line);
	    else if (!strcmp("j", keyword))
		set_vector_from_input_line(jerk, input_line);
	    else if (!strcmp("pot", keyword))
		sscanf(input_line, "%*s%*s%lf", &pot);
	    else
		add_story_line(dyn_story, input_line);
	}
    }

    return s;
}

ostream& sdyn3::print_dyn_story(ostream& s,
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
    put_real_vector(s, "  a  =  ", acc);
    put_real_vector(s, "  j  =  ", jerk);

    if (!short_output) {

	put_real_number(s, "  pot  =  ", pot);

	put_real_number(s, "  min_nn_dr2  =  ", min_nn_dr2);
	put_integer(s, "  min_nn_label  =  ", min_nn_label);
	put_real_number(s, "  min_min_ssd  =  ", min_min_ssd);
	put_integer(s, "  n_ssd_osc  =  ", n_ssd_osc);
	put_real_number(s, "  n_steps  =  ", n_steps);
	put_real_number(s, "  e_tot_init  =  ", e_tot_init);
	put_real_number(s, "  de_tot_abs_max  =  ", de_tot_abs_max);

	if (dyn_story)
	    put_story_contents(s, *dyn_story);
    }

    put_story_footer(s, DYNAMICS_ID);
    
    return s;
}

#else

main(int argc, char** argv)
{
    check_help();

    sdyn3 * b;
    while (b = (sdyn3 *) get_node(cin, new_sdyn3)){
	put_node(cout,*b);
	pp2(b);
    }
    cerr << "Normal exit\n";
}

#endif
