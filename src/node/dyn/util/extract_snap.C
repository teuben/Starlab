
//// extract_snap:  Find and print the last snapshot, or the first
////                snapshot following a specified time.
////
//// Options:    -c    add a comment to the output snapshot [false]
////             -t    specify time [last snapshot]
////             -n    specify number of snapshots to extract [1]
////             -v    set verbose mode [off]
////
//// If -n is set to 0 and -v is set, then all we do is count snapshots.

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = FALSE;	// if TRUE, a comment given on command line

    real  t_extract;
    bool  t_flag = FALSE;	// if TRUE, a time was specified

    bool  v_flag = FALSE;	// if TRUE, print snap times as read

    int   n = 1;
    bool  n_flag = false;	// if TRUE, a number was specified

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:n:t:v";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c': c_flag = TRUE;
	    	      comment = poptarg;
		      break;
	    case 't': t_flag = TRUE;
		      t_extract = atof(poptarg);
		      break;
	    case 'n': n_flag = TRUE;
		      n = atoi(poptarg);
		      break;
	    case 'v': v_flag = TRUE;
		      break;
	    case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
		      exit(1);
	}            

    if (n < 0) err_exit("n < 0 specified.");
    if (n > 1 && !t_flag) {
	warning("n > 1 but no time specified -- 0 assumed.");
	t_extract = 0;
	t_flag = true;
    }

    dyn *b = NULL, *bp = NULL;
    int i = 0;

    while (b = get_dyn(cin)) {

	real time = b->get_system_time();
	i++;

	if (v_flag) cerr << "Snap time #" << i << " = " << time << endl;

	if (t_flag && time >= t_extract) {

	    if (n > 0) {

		if (c_flag == TRUE)
		    b->log_comment(comment);

		b->log_history(argc, argv);
		put_dyn(cout, *b);
	    }

	    if (--n <= 0) exit(0);
	}

	if (bp != NULL)		// Hmmm... if (!bp) fails here on merlot!
	    rmtree(bp);

	bp = b;
    }

    if (n > 0 && !t_flag) {
	bp->log_history(argc, argv);
	put_dyn(cout, *bp);
    }

}

#endif
