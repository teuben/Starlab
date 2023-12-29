
//// Sample Starlab tool:  Bring all positions and velocities
////                       to the center-of-mass frame.
////
//// Options:     -c    add a comment to the output snapshot [false]

#include "dyn.h"

local void shift_to_com(dyn *b)
{
    // Compute the center of mass.

    real total_mass = 0;
    vector com_pos  = 0;
    vector com_vel  = 0;

    for_all_daughters(dyn, b, d) {

	total_mass += d->get_mass();
	com_pos    += d->get_mass() * d->get_pos();
	com_vel    += d->get_mass() * d->get_vel();

    }	

    com_pos /= total_mass;
    com_vel /= total_mass;

    // Shift all positions and velocities.

    for_all_daughters(dyn, b, d) {
        d->inc_pos(-com_pos);
        d->inc_vel(-com_vel);
    }

    // Add an entry to the root log story.

    char tmp[128];
    sprintf(tmp, "  modified system center of mass at time %f",
	    b->get_system_time());
    b->log_comment(tmp);

    putvq(b->get_log_story(), "old_com_pos", com_pos);
}

main(int argc, char *argv[])
{
    check_help();

    // Parse the command line.

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    bool c_flag = FALSE;
    char *comment;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
        switch(c) {

            case 'c': c_flag = TRUE;
	    	      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
        }            

    dyn *b;

    while (b = get_dyn(cin)) {

        if (c_flag) b->log_comment(comment);
        b->log_history(argc, argv);

        shift_to_com(b);

        put_dyn(cout, *b);        
        rmtree(b);
    }
}
