
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

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

    // On entry, we have just read the "(Dynamics" [or "(D"] line
    // signifying the start of dyn_story input.  Keep reading and
    // storing until we encounter the corresponding closing
    // [END_DYNAMICS] line.

    while (get_line(s,input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

    	if (0) {   // trick to keep the else if() statements homogeneous

    	    }else if(!strcmp("m",keyword)){
		mass = strtod(val, NULL);
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
    // Modifications by Steve (5/01) to streamline output.

    put_real_number(s, "  m  =  ", mass);

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
