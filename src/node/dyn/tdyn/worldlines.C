
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

////    worldlines:  Create output snapshots from input worldbundles.
////                 Produce a snapshot or diagnostic starting at time t,
////                 optionally with time step dt.
////
////    options:     -d    output interval (negative integer ==> 2^[-d])   [0]
////                 -D    diagnostic: no snapshsot                       [no]
////                 -F    specify input file                          [stdin]
////                       *** file input is much faster (why?) ***
////                 -t    initial output time                 [start of data]
////                 -v    verbosity level                                 [1]
////
////    Snapshots are generated as soon as the data are available
////    -- i.e. we don't read in the entire dataset first.
////    Only one snapshot is produced if the time step is 0 (default).
////
////    Note: Output is in pdyn format, which currently cannot be used
////    directly as input to kira.
////
////    Created by  SPZ at MIT in Dec 2000
////    Modified by SLWM at DU in Dec 2003
//
//----------------------------------------------------------------------

#include "worldline.h"

#ifdef TOOLBOX

worldbundle *get_next_worldbundle(bool file, char *filename, int verbose)
{
    static int nh = 0;
    static ifstream *s = NULL;
    worldbundle *wb;

    if (file) {
	if (!s) {
	    s = new ifstream(filename);
	    if (!s) {
		cerr << "Data file " << filename << " not found." << endl;
		exit(1);
	    }
	}
	wb = read_bundle(*s, verbose>2);

    } else
	wb = read_bundle(cin, verbose>2);

    nh++;

    if (verbose > 1) {

	// Brief info on the new worldline.

	real t = wb->get_t_min();
	int nw = wb->get_nw(), ns = count_segments(wb),
	ne = count_events(wb);
	cerr << "worldbundle " << nh << ": "
	     << nw << " worldlines, "
	     << ns << " segments, "
	     << ne << " events, t = "
	     << wb->get_t_min() << " to " << wb->get_t_max()
	     << endl;
    }

    return wb;
}

local void convert_relative_to_absolute(dyn* b)
{
    if (b->get_parent()) b->inc_pos(b->get_parent()->get_pos());
    for_all_daughters(dyn, b, bb) convert_relative_to_absolute(bb);
}

main(int argc, char *argv[])
{
    bool t_flag = false;
    real t;			// first output time

    real dt = 0;

    bool file = false;
    char infile[128];

    bool diag = false;
    int verbose = 1;

    extern char *poptarg;
    int c;
    char* param_string = "d:DF:t:v:";
  
    check_help();
  
    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
      
	    case 'd': dt = atof(poptarg);
	    	      if (dt < 0) dt = pow(2.0, dt);
		      break;
	    case 'D': diag = true;
	    	      break;
	    case 'F': file = true;
	    	      strcpy(infile, poptarg);
	    	      break;
	    case 't': t_flag = true;
		      t = atof(poptarg);
		      break;
	    case 'v': verbose = atoi(poptarg);
		      break;
	    case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
		      exit(1);
	}

    // Verbosity levels:	0	no output
    //				1	print time when writing snapshot
    //				2	1+brief info on worldlines
    //				3	1+long info on worldlines


    worldbundle *wb = get_next_worldbundle(file, infile, verbose);
    if (!wb) {
	cerr << "No data!" << endl;
	exit(0);
    }

    if (!t_flag) t = wb->get_t_min();

    while (1) {
	if (t >= wb->get_t_min() && t <= wb->get_t_max()) {

	    pdyn *b = create_interpolated_tree2(wb, t, true);

	    if (b) {
		if (diag) {
		    vec com_pos, com_vel;
		    compute_com(b, com_pos, com_vel);
		    PRC(t); PRL(com_pos);
		    PRL(b->get_pos());
		} else {
		    if (verbose) {
			cerr << "Writing snapshot at "; PRL(t);
		    }
		    put_pdyn(b);
		}
	    }

	    if (dt == 0) break;
	    t += dt;

	} else
	    if (!(wb = get_next_worldbundle(file, infile, verbose)))
		break;
    }
}

#endif
