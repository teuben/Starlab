
//// dump_worldlines:  print out details of worldbundle structure.
////
//// Options:     -F    input file [world.dat]

//.............................................................................
//    version 1:  Aug 2001   Steve McMillan

#ifdef TOOLBOX
#include "worldline.h"

local void print_worldlines(worldbundleptr wh[], int nh)
{
    for (int ih = 0; ih < nh; ih++) {
	cerr << "worldbundle " << ih << ":" << endl;
	wh[ih]->dump(4);
    }
}

main(int argc, char** argv)
{
    char infile[128];
    strcpy(infile, "world.dat");

    check_help();

    extern char *poptarg;
    char* params = "F:";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'F': strcpy(infile, poptarg);
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
    while (nh < 1024 && (wb = read_bundle(s, false))) wh[nh++] = wb;

    // PRL(mass_scale_factor());

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

    // Details:

    print_worldlines(wh, nh);
}
#endif
