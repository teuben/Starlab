
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// loop_tdyn: loop through a tdyn structure, mimicking the action of
////            xstarplot22, but displaying nothing.
////
////    options:
////             -D    time interval between frames [0.015625 = 1/64]
////
////    input: worldbundle output file from kira
////
////    Author: SLWM, October 2003
//
//----------------------------------------------------------------------

#include "worldline.h"
#include "dyn_util.h"

#define DYN dyn
#define DYNPTR dynptr

main(int argc, char** argv)
{
    int verbose = 0;
    real dt = 0.015625;		// powers of 2 are preferred, but not essential

    char infile[128];
    strcpy(infile, "run.out");

    check_help();

    extern char *poptarg;
    char* params = "D:";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'D': dt = atof(poptarg);	// delay between frames [0.01]
		      if (dt < 0)		// usual convention
			  dt = pow(2.0, dt);
		      break;
            case '?': params_to_usage(cerr, argv[0], params);
	              get_help();
		      exit(0);
	}

    ifstream s(infile);
    if (!s) {
	cerr << "Data file " << infile << " not found." << endl;
	exit(1);
    }

    // Read in an array of worldbundles (a "worldhistory", when this
    // becomes the next level of structure in the world hierarchy).
    // Eventually, the display can start as soon as the first
    // worldbundle is ready, and the remainder can be read in
    // asynchronously.

    worldbundleptr wb, wh[1024];

    int nh = 0;
    while (nh < 1024 && (wb = read_bundle(s, verbose))) wh[nh++] = wb;

    cerr << endl << "statistics on " << nh << " worldbundle";
    if (nh != 1) cerr << "s";
    cerr << ":" << endl;

    int nwtot = 0, nstot = 0, netot = 0;
    for (int ih = 0; ih < nh; ih++) {
	wb = wh[ih];
	real t = wb->get_t_min();
	int nw = wb->get_nw(), ns = count_segments(wb), ne = count_events(wb);
	cerr << "worldbundle " << ih << ": "
	     << nw << " worldlines, "
	     << ns << " segments, "
	     << ne << " events, t = "
	     << wb->get_t_min() << " to " << wb->get_t_max()
	     << endl;
	nwtot += nw;
	nstot += ns;
	netot += ne;
    }
    cerr << "totals: " << nwtot << " worldlines, "
	 << nstot << " segments, " << netot << " events"
	 << endl << endl;

    // Now display the data.

    bool eod = false;
    int step_mode = 0;
    int init = 0;

    int ih = 0;
    wb = wh[ih];
    real t = wb->get_t_min();

    bool max_cm_set = false;

    while (1) {

	DYN *root = create_interpolated_tree2(wb, t);

	if (root) {
	    PRC(t); PRL(root->get_mass());
	}

	t += dt;

	// Move to the next worldbundle, or loop, as necessary.

#define EPS 1.e-12

	if (dt > 0 && t > wb->get_t_max() + EPS) {
	    if (++ih >= nh) {
		ih = 0;
		wb = wh[ih];
		t = wb->get_t_min();
	    } else
		wb = wh[ih];
	}

	// May be moving backwards, so check that too.

	if (dt < 0 && t < wb->get_t_min() - EPS) {
	    if (--ih < 0) {
		ih = nh-1;
		wb = wh[ih];
		t = wb->get_t_max();
	    } else
		wb = wh[ih];
	}

    }
}
