#include "dyn.h"

#ifdef TOOLBOX

main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case '?':	params_to_usage(cerr, argv[0], param_string);
	    		get_help();
	    		exit(1);
        }            

    dyn *b;
    real prev_time = -VERY_LARGE_NUMBER;
    cerr.precision(HIGH_PRECISION);

    while (b = get_dyn(cin)) {

	real time = b->get_system_time();

	if (prev_time > -VERY_LARGE_NUMBER) {
	    real dt = time - prev_time;
	    PRC(time); PRL(dt);
	} else
	    PRL(time);
	prev_time = time;

	rmtree(b);
    }
}

#endif
