/*
 *  precision.C:  Set output precision for starlab programs.
 *
 *.............................................................................
 *    version 1:  Nov 1994   Steve McMillan
 *    version 2:  
 *.............................................................................
 *  non-local functions: 
 *    set_starlab_precision
 *.............................................................................
 */

//// precision:  Set Starlab output precision using environment
////             variable STARLAB_PRECISION, if available.
////
//// No options

// Probably should be implemented as a class (like kira_options)...

#include "stdinc.h"

#ifndef TOOLBOX

static int precision = -1;

// Accessor:

int get_starlab_precision() {return precision;}

local void get_default_precision()
{
    // Establish the default Starlab precision setting.
    // First try the environment, then use a default value.

    precision = -1;

    char* s = getenv("STARLAB_PRECISION");

    if (s) {
	int p = atoi(s);
	if (p > 0) precision = p;
    }

    if (precision < 0)

	// Note by J.M. 97-Jan-5:  changed 16 to 18,
	// to try to guarantee the correctness of the
	// last bit (not thoroughly tested).

	precision = 18;
}

int set_starlab_precision(ostream &s)
{
    if (precision < 0) get_default_precision();
    return s.precision(precision);	// convenient to return the
					// old "precision" setting
}

int adjust_starlab_precision(int p)
{
    if (p >= 0)
	precision = p;
    else
	get_default_precision();	// reread the environment

    return precision;
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
