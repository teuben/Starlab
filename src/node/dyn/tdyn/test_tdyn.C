#include "tdyn.h"

#define WHICH tdyn

main(int argc, char** argv)
{
    int n = 1000;

    extern char *poptarg;
    char* params = "n:";
    int   c;

    while ((c = pgetopt(argc, argv, params)) != -1)
	switch(c) {
	    case 'n': n = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], params);
	              get_help();
		      exit(0);
	}

    cerr << "pause..." << flush;
    char tmp;
    cin >> tmp;

    for (int i = 0; i < n; i++) {
	WHICH *b = new WHICH(NULL, NULL, false);
    }

    cerr << "pause..." << flush;
    cin >> tmp;
}
