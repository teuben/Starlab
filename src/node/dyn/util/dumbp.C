
//// dumbp:  dump out an N-body system in a dumb format suitable
////         for digestion by NBODY1-5
////
//// Options:     -p    specify precision of output [6 sig. fig.]
////              -t    also add the time

#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char** argv)
{
    dyn *root, *ni;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "p:t";

    int p = STD_PRECISION;
    int t = 0;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'p': p = atoi(poptarg);
		      break;
	    case 't': t = 1;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }

#ifdef BAD_GNU_IO
    char prec[32];
    sprintf(prec, "%%.%df ", p);
#else
    cout.precision(p);
#endif

    while (root = get_dyn(cin)) {

	for (ni = root->get_oldest_daughter(); ni != NULL;
	     ni = ni->get_younger_sister()) {

#ifdef BAD_GNU_IO

            if (t) printf(prec,ni->get_system_time());
	    printf(prec, ni->get_mass());
	    vector temp;
	    temp = ni->get_pos();
	    int k;
	    for (k = 0; k < 3; k++) printf(prec,temp[k]);
	    temp = ni->get_vel();
	    for (k = 0; k < 3; k++) printf(prec,temp[k]);
	    printf("\n");
#else
            if (t) cout << ni->get_system_time() << " ";
	    cout << ni->get_mass() << " "
		 << ni->get_pos()  << " "
		 << ni->get_vel()
		 << endl;
#endif

	}

	rmtree(root);
    }
}

#endif

// end of: dumbp.C
