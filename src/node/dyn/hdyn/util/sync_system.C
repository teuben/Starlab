
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// sync_system:  Synchronize system to the system time.
////               This rather dumb function forces the times of all nodes
////               to the system time.
////               It solves the problem of non-restartability ffter having
////               run with slow KS (-K option in kira).
////               The problem with this solution is that it randomly shifts
////               the orbital phase of KS binaries.
////               On the other hand, it allows us to restarts kira.
////               It is recommended to use this function sparsely.

//   version 1:  March 2002   Simon Portegies Zwart

#include "hdyn.h"

#ifdef TOOLBOX

void sync_system(hdyn *b, xreal time) {

    real dtmax = 1;
    real ttime = (real)time;
    b->set_system_time(ttime);

    PRC(dtmax);PRC(ttime);PRL(b->get_system_time());
    cerr << "checking timestep consistency...";
    while(fmod(ttime, dtmax) != 0) {
	dtmax /= 2;
	PRC(dtmax);
	if(dtmax<VERY_SMALL_NUMBER) 
	    err_exit("sync_system: timestep too small");
    }
    cerr << "OK,  "; PRL(dtmax);

    xreal ct;
    for_all_nodes(hdyn, b, bb) {
	ct = bb->get_time();
	if(ct != time) {
	    cerr << "Time of particle " << bb->format_label()
		 << " not equal to system time (" << ct << " != "
		 << time << ")" << endl;
	    bb->set_time(time);
	    bb->set_timestep(Starlab::min(bb->get_timestep(), dtmax));
	}
    }
}

main(int argc, char ** argv)
{
    bool c_flag = false;
    char  *comment;
    xreal time = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "c:t:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
  	    case 't':	atof(poptarg);
			break;
            case '?':	params_to_usage(cerr, argv[0], param_string);
	    		get_help();
	    		exit(1);
        }            

    hdyn *b;

    while (b = get_hdyn(cin)) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	time = b->get_system_time();

	cerr << "Synchronizing system to time " << time << endl;
	sync_system(b, time);

	put_hdyn(b);
	rmtree(b);
    }
}

#endif
