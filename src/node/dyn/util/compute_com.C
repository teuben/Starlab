
//// compute_com:  Determine the center of mass position and velocity of
////               the input N-body system.
////
////               Center of mass position and velocity are written to the
////               dyn story of the top-level node; they are also optionally
////               returned as function arguments in the library version.
////
//// Options:     -c    add a comment to the output snapshot [false]

#include "dyn.h"

#ifndef TOOLBOX

void compute_com(dyn *b, vector& pos, vector& vel)
{
    real total_mass = 0;

    pos = 0;
    vel = 0;
    
//  for_all_leaves(dyn, b, d)
    for_all_daughters(dyn, b, d) {
	total_mass += d->get_mass();
	pos += d->get_mass() * d->get_pos();
	vel += d->get_mass() * d->get_vel();
    }	

    pos /= total_mass;
    vel /= total_mass;

    putrq(b->get_dyn_story(), "com_time", b->get_system_time());
    putvq(b->get_dyn_story(), "com_pos", pos);
    putvq(b->get_dyn_story(), "com_vel", vel);
}

void compute_com(dyn *b)
{
    vector pos, vel;
    compute_com(b, pos, vel);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use compute_com() as a tool
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

        compute_com(b);

	// Write system to stdout and get next system (if any).

        put_dyn(cout, *b);
	rmtree(b);
	b = get_dyn(cin);
    }
}

#endif

// endof: compute_com.C
