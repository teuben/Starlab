
//// set_com:  Modify positions and velocities to set the center-of-mass
////           position and velocity.  Writes to the root dyn story.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -r    specify center of mass position [(0,0,0)]
////              -v    specify center of mass velocity [(0,0,0)]

//   version 1:  Aug 2001   Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void dyn::set_com(vector pos, vector vel)	// defaults = 0
{
    vector com_pos;
    vector com_vel;

    compute_com(this, com_pos, com_vel); 

    vector dpos = pos - com_pos;
    vector dvel = vel - com_vel;

    for_all_daughters(dyn, this, bb) {
	bb->inc_pos(dpos);
	bb->inc_vel(dvel);
    }

    // Correct entries in dyn story.

    putvq(get_dyn_story(), "com_pos", pos);
    putvq(get_dyn_story(), "com_vel", vel);

    // If a "center" for an external field has been specified,
    // we should probably adjust that too.

    if (external_field) {

	if (find_qmatch(get_log_story(), "kira_plummer_center")) {

	    // Plummer field:

	    vector c = getvq(get_log_story(), "kira_plummer_center");
	    if (abs(c) < VERY_LARGE_NUMBER) {
		c += dpos;
		putvq(get_log_story(), "kira_plummer_center", c);
	    }
	}

    }
}

#else

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    vector r = 0;
    vector v = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "c:r:::v:::";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
	    case 'r':	r = vector(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;
	    case 'v':	v = vector(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;
            case '?':	params_to_usage(cerr, argv[0], param_string);
	    		get_help();
	    		exit(1);
        }            

    dyn *b;

    while (b = get_dyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        b->set_com(r, v);
	cerr << "set_com:  "; PRC(r); PRL(v);

	put_dyn(cout, *b);	
	rmtree(b);
    }
}

#endif

// endof: set_com

