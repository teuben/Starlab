//// readbundle:  read and print statistics on a series of worbundles,
////              and print out data on the interpolated tree at fixed
////              time intervals
////
//// Options:    -d    time step [0.0625]
////             -F    input file [run.out]

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
}

main(int argc, char** argv)
{
    check_help();

    real dt = 0.0625;
    char infile[128];
    strcpy(infile, "run.out");

    extern char *poptarg;
    char* params = "d:F:";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'd': dt = atof(poptarg);
		      if (dt < 0)		// usual convention
			  dt = pow(2.0, dt);
		      break;
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
    while (nh < 1024 && (wb = read_bundle(s, false))) wh[nh++] = wb;

    print_worldline_stats(wh, nh);

    int ih = 0;
    wb = wh[ih];
    real t = wb->get_t_min();

    while (1) {

	real cpu = cpu_time();
	pdyn *root = create_interpolated_tree2(wb, t);
	cpu = cpu_time() - cpu;
	if (abs(cpu) < 1.e-6) cpu = 0;

	// Nominal output:

	PRC(t); PRC(root); PRL(root->format_label());
	PRI(4); PRC(wb->get_nw()); PRC(root->n_daughters()); PRL(cpu);

	for_all_nodes(pdyn, root, p)
	    if (p->get_worldline_index() < 0) {
		PRC(p); PRC(p->format_label()); PRL(p->get_worldline_index());
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
