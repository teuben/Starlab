#include "worldline.h"
#include "dyn_util.h"

#define DYN pdyn
#define DYNPTR pdynptr

main(int argc, char *argv[])
{
    char infile[128];
    strcpy(infile, "run.out");

    bool verbose = false;

    check_help();

    extern char *poptarg;
    char* params = "F:v";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'F': strcpy(infile, poptarg);
	    	      break;
	    case 'v': verbose = !verbose;
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

    // Read in an array of worldbundles.

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

    real t_min = wh[0]->get_t_min();
    real t_max = wh[nh-1]->get_t_max();

    while (1) {

	real t;
	cerr << "time = " << flush;
	cin >> t;
	if (cin.eof()) exit(0);

	if (t < t_min || t > t_max) continue;

	int ih;
	for (ih = 0; ih < nh; ih++)
	    if (wh[ih]->get_t_max() >= t) break;

	DYN *root = create_interpolated_tree2(wh[ih], t);

	if (root) {
	    cerr << "worldbundle " << ih << ", "
		 << root->n_daughters() << " top-level nodes, "
		 << root->n_leaves() << " leaves, total mass = "
		 << root->get_mass()
		 << endl;
	}

    }
}
