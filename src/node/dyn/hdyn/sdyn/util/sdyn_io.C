
//// sdyn_io: Starlab sdyn-specific I/O functions.
////
//// Options: none

#include "sdyn.h"
#include "util_io.h"

#ifndef TOOLBOX

istream & sdyn::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), strcmp(END_DYNAMICS, input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	char *val = getequals(input_line, keyword);

    	if(0){   // trick to keep the else if() statements homogeneous
    	    }else if(!strcmp("m",keyword)){
		mass = strtod(val, &val);
	    }else if(!strcmp("r",keyword)){
		set_vector_from_input_line(pos,input_line);
	    }else if(!strcmp("v",keyword)){
		set_vector_from_input_line(vel,input_line);
	    }else if(!strcmp("t",keyword)){
		time = strtod(val, &val);
	    }else if(!strcmp("dt",keyword)){
		timestep = strtod(val, &val);
	    }else if(!strcmp("a",keyword)){
		set_vector_from_input_line(acc,input_line);
	    }else if(!strcmp("j",keyword)){
		set_vector_from_input_line(jerk,input_line);
	    }else if(!strcmp("pot",keyword)){
		pot = strtod(val, &val);
	    }else{
		add_story_line(dyn_story, input_line);
	    }
	}
    return s;
}

ostream& sdyn::print_dyn_story(ostream& s,
			       bool print_xreal,	// default = true
			       int short_output)	// default = 0
							// (not implemented)
{
    put_story_header(s, DYNAMICS_ID);

    put_real_number(s, "  t  =  ", time);
    put_real_number(s, "  m  =  ", mass);
    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);
    put_real_vector(s, "  a  =  ", acc);

    if (!short_output) {
	put_real_vector(s, "  j  =  ", jerk);
	put_real_number(s, "  pot  =  ", pot);

	put_real_number(s, "  min_nn_dr2  =  ", min_nn_dr2);
	put_integer(s, "  min_nn_label  =  ", min_nn_label);
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

    sdyn * b;
    while (b = (sdyn *) get_node(cin, new_sdyn)){
        put_node(cout,*b);
	pp2(b);
    }
    cerr << "Normal exit\n";
}

#endif
