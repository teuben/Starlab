//// scanbinary:  read and print binary statistics on a series of worbundles.
////
//// Options:    -F    input file                                    [stdin]
////                   *** file input is much faster (why?) ***
////
//// Program reads a collection of worldbundles, then repeatedly prompts
//// for the name of a particle or node and prints out the interaction
//// history of that object, in lines of the form
////
////                 time: top-level node name
////
//// End with exit, quit, or ^D.
////
//// Created by Steve, 9/2001.

#include "worldline.h"

local void print_interactions(ostream& str,
			      worldbundleptr wh[], int nh,
			      char *id)
{
    unique_id_t curr_id = -1;		// attempt to deal with repeats
    unique_id_t prev_id = -1;		// (1-cycles) and 2-cycles
    unique_id_t prev_prev_id = -1;

    real trep = 0;
    int  nrep = 0;
    int order = 0;

    for (int ih = 0; ih < nh; ih++) {
	worldbundle *wb = wh[ih];
	worldline *w = wb->find_worldline(id);
	// PRC(ih); PRL(w);
	if (w) {
	    int is = 0;
	    segment *s = w->get_first_segment();
	    // PRC(is); PRL(s);

	    while (s) {
		tdyn *b = s->get_first_event();
		// PRC(b); PRL(b->format_label());

		if (b) {

		    // Find the top-level node.

		    tdyn *p = b;
		    while (p && p->get_parent()
			   && strcmp(p->get_parent()->format_label(),
				     "root")) {
			// PRC(p); PRL(p->format_label());
			p = p->get_parent();
		    }

		    // See whether to print anything, and keep track of cycles.

		    unique_id_t this_id = unique_id(p);

		    // Note: 1-cycle has this_id = curr_id; 2-cycle has
		    // this_id = prev_id and curr_id = prev_prev_id...

		    bool print = (ih == nh-1 && s->get_next() == NULL);
		    
		    if (this_id == curr_id) {
			order = 1;
			nrep++;
			trep = p->get_time();
		    } else {
			if (this_id == prev_id && curr_id == prev_prev_id) {
			    order = 2;
			    nrep++;
			    trep = p->get_time();
			} else
			    print = true;
		    }

		    if (print) {
			int prec = str.precision(10);
			if (nrep > 0 && order > 1)
			    str << trep << ": [" << order
				<< "-cycle, rep = " << nrep << "]"
				<< endl;
			str << p->get_time() << ": "
			    << p->format_label()
			    << endl;
			nrep = 0;
			str.precision(prec);
		    }

		    prev_prev_id = prev_id;
		    prev_id = curr_id;
		    curr_id = this_id;
		}
		s = s->get_next();
		is++;
	    }
	}
    }
}

main(int argc, char** argv)
{
    check_help();

    bool file = false;
    char infile[128];
    int verbose = 0;

    extern char *poptarg;
    char* params = "F:v";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {

	    case 'F': file = true;
	    	      strcpy(infile, poptarg);
	    	      break;
	    case 'v': verbose = 1;
	    	      break;
            case '?': params_to_usage(cerr, argv[0], params);
	              get_help();
		      exit(0);
	}

    int nh = 1024;
    worldbundleptr *wh = new worldbundleptr[nh];

    if (file) {
	ifstream s(infile);
	if (!s) {
	    cerr << "Data file " << infile << " not found." << endl;
	    exit(1);
	}
	read_bundles(s, wh, nh, verbose);

    } else
	read_bundles(cin, wh, nh, verbose);


    preload_pdyn(wh, nh, false, false);

    while (1) {
	char id[64];
	cout << "name: " << flush; cin >> id;
	if (strlen(id) <= 0
	    || strstr(id, "exit")
	    || strstr(id, "quit")) {
	    if (strlen(id) <= 0) cout << endl;
	    break;
	}
	print_interactions(cout, wh, nh, id);
    }
}
