
//// node_io:  Starlab node I/O functions.
////
//// Options:  none

#include "node.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize all static node data here (?)

node * node::root           = NULL;

istream& node::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while (get_line(s,input_line), strcmp(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];

	sscanf(input_line, "%s%s", keyword,should_be_equal_sign);

	if (strcmp("=", should_be_equal_sign)) {
	    cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
	    exit(1);
	}

    	if (0) {   // trick to keep the else if() statements homogeneous

    	    }else if(!strcmp("m",keyword)){
		sscanf(input_line,"%*s%*s%lf",&mass);
	    }else{
		add_story_line(dyn_story, input_line);
	    }
    }
    return s;
}

ostream& node::print_dyn_story(ostream& s,
			       bool print_xreal,	// default = true
			       int short_output)	// default = 0
{
    put_story_header(s, DYNAMICS_ID);

    put_real_number(s, "  m  =  ", mass);

    if (!short_output && dyn_story)
        put_story_contents(s, *dyn_story);		// HP OK -- story.C

    put_story_footer(s, DYNAMICS_ID);
    
    return s;
}

#else

main(int argc, char** argv)
{
    check_help();

    node * b;
    while (b = get_node(cin)){
        put_node(cout,*b);
	pp2(b);
    }
    cerr << "Normal exit\n";
}
#endif
