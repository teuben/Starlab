/*
 *  setup.C:	Basic setup for starlab programs.
 *
 *		currently, only sets output precision...
 *
 *  Name shortened to satisfy some library managers.
 *
 *.............................................................................
 *    version 1:  Nov 1994   Steve McMillan
 *    version 2:  
 *.............................................................................
 *  non-local functions: 
 *    set_starlab_precision
 *.............................................................................
 */

//// setup:  check Starlab function to set output precision.
////
//// No options

#include "stdinc.h"

#ifndef TOOLBOX

static int precision = -2;

int get_starlab_precision()
{
    if (precision == -2) {
	char* s = getenv("STARLAB_PRECISION");
	if (s)
	    precision = atoi(s);
	else
	    //	precision = -1;
	    //  Note by J.M. 96-Aug-5 I changed this part
	    //  so that the default precision is 16, not 
	    //  the system default.
	    //  precision = 16;
	    //  Note by J.M. 97-Jan-5 I changed 16 to 18,
	    //  to guarantee the correctness of the last
	    //  bit. (not tested though).
	    precision = 18;
    }

    return precision;
}

void set_starlab_precision(ostream& o)
{
    get_starlab_precision();
    if (precision > 0) o.precision(precision);
}

#else

main(int argc, char** argv)
{
    check_help();
    set_starlab_precision(cout);
    
    real x = sqrt(2.0);
    cout << "x = " << x << endl;
}

#endif
