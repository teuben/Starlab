
//// dumbp:  dump out an N-body system in a dumb format suitable
////         for digestion by NBODY1-5
////
//// Options:     -p    specify precision of output [6 sig. fig.]
////              -t    include time in output [no]

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
    bool time = false;

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'p': p = atoi(poptarg);
		      break;
	    case 't': time = !time;
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

    while (root = get_dyn()) {

	for_all_daughters(dyn, root, ni) {

	    real t = getrq(ni->get_dyn_story(), "t");

#ifdef BAD_GNU_IO

	    if (time) printf(prec, t);

	    printf("%d ", ni->get_index());
	    printf(prec, ni->get_mass());
	    vector temp;
	    temp = ni->get_pos();
	    int k;
	    for (k = 0; k < 3; k++) printf(prec,temp[k]);
	    temp = ni->get_vel();
	    for (k = 0; k < 3; k++) printf(prec,temp[k]);
	    printf("\n");

#else

	    if (time) cout << t << " ";
	    cout << ni->get_index() << " "
		 << ni->get_mass() << " "
		 << ni->get_pos()  << " "
		 << ni->get_vel()
		 << endl;

#endif

	}

#ifdef BAD_GNU_IO
	    printf("\n");
#else
	    cout << endl;
#endif

	rmtree(root);
    }
}

#endif

// end of: dumbp.C
