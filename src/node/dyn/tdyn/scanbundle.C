//// scanbundle:  read and print statistics on a series of worbundles.
////
//// Options:    -F    input file [run.out]

#include "worldline.h"

local void print_worldline_stats(worldbundleptr wh[], int nh)
{
    cerr << endl << "statistics on " << nh << " worldbundle";
    if (nh != 1) cerr << "s";
    cerr << ":" << endl;

    int nwtot = 0, nstot = 0, netot = 0;
    for (int ih = 0; ih < nh; ih++) {
	worldbundleptr wb = wh[ih];
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

    for (int ih = 0; ih < nh; ih++)
	wh[ih]->check();
}

main(int argc, char** argv)
{
    check_help();

    char infile[128];
    strcpy(infile, "run.out");

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

    worldbundleptr wb, wh[1024];

    int nh = 0;
    while (nh < 1024 && (wb = read_bundle(s, 1))) wh[nh++] = wb;

    for (int i = 0; i < 5; i++) cerr << endl;

    print_worldline_stats(wh, nh);
}
