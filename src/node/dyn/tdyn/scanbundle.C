//// scanbundle:  read and print statistics on a series of worldbundles.
////
//// Options:    -F    specify input file                        [stdin]
////                   *** file input is much faster (why?) ***
////             -v    more (too?) detailed output                  [no]

#include "worldline.h"

main(int argc, char** argv)
{
    check_help();

    bool file = false;
    char infile[128];
    int verbose = 0;

    extern char *poptarg;
    char* params = "F:v";
    int   c;

    while ((c = pgetopt(argc, argv, params,
		    "$Revision$", _SRC_)) != -1)
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
	read_bundles(s, wh, nh, verbose+1);

    } else
	read_bundles(cin, wh, nh, verbose+1);


    for (int ih = 0; ih < nh; ih++) {
	if (verbose) {
	    cerr << "worldbundle " << ih << ":" << endl;
	    wh[ih]->dump(4);
	} else
	    wh[ih]->check();
    }
}
