
//// get_std_center: Determine the "standard" center of the input N-body
////                 system, determined as follows:
////
////                 (1) if densities are computed and up to date, use
////                     the density center (but *don't* recompute the
////                     densities), as determined by compute_max_cod(),
////
////                 (2) otherwise, use the modified center of mass,
////                     as returned by compute_mcom().
////
////               Center position and velocity are written to the dyn
////               story of the top-level node; they are also optionally
////               returned as function arguments in the library version.
////
//// Options:     -c    add a comment to the output snapshot [false]

#include "dyn.h"

#ifndef TOOLBOX

int get_std_center(dyn *b, vector& pos, vector& vel)
{
    int which = 0;

    // See if densities are available and up to date.

    if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {

	if (getrq(b->get_dyn_story(), "density_center_time")
		!= b->get_system_time())

	    warning("get_std_center: neglecting out-of-date density center");

	else {

	    // Assume that density_center_vel exists if density_center_pos
	    // is OK.

	    pos = getvq(b->get_dyn_story(), "density_center_pos");
	    vel = getvq(b->get_dyn_story(), "density_center_vel");

	    which = 1;
	}
    }

    if (which == 0) {

	if (find_qmatch(b->get_dyn_story(), "com_mpos")) {

	    // Try modified center of mass instead.

	    if (getrq(b->get_dyn_story(), "mcom_time")
		    != b->get_system_time())

		warning("get_std_center: neglecting out-of-date mcom");

	    else {

		pos = getvq(b->get_dyn_story(), "mcom_pos");
		vel = getvq(b->get_dyn_story(), "mcom_vel");
	    }

	} else

	    // Compute the modified center of mass.

	    compute_mcom(b, pos, vel);

	which = 2;
    }

    putrq(b->get_dyn_story(), "std_center_time", b->get_system_time());
    putvq(b->get_dyn_story(), "std_center_pos", pos);
    putvq(b->get_dyn_story(), "std_center_vel", vel);

    return which;
}

int get_std_center(dyn *b)
{
    vector pos, vel;
    return get_std_center(b, pos, vel);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use get_std_center() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    dyn * b;
    bool  c_flag = FALSE;       // if TRUE: a comment given on command line

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn(cin)) == NULL)
       err_exit("compute_com: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        get_std_center(b);

	// Write system to stdout and get next system (if any).

        put_dyn(cout, *b);
	rmtree(b);
	b = get_dyn(cin);
    }
}

#endif

// endof: get_std_center.C
