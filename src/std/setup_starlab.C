/*
 *  setup.C:	Basic setup for starlab programs.
 *
 *		currently, only sets output precision...
 *
 *.............................................................................
 *    version 1:  Nov 1994   Steve McMillan
 *    version 2:  
 *.............................................................................
 *  non-local functions: 
 *    set_starlab_precision
 *.............................................................................
 */

#include "stdinc.h"

static int precision = 0;

#ifndef TOOLBOX

void set_starlab_precision(ostream& o)
{
    if (!precision) {
	char* s = getenv("STARLAB_PRECISION");
	if (s)
	    precision = atoi(s);
	else
	    precision = -1;
    }

    if (precision > 0) o.precision(precision);
}

#else

main()
{
   set_starlab_precision(cout);

   real x = sqrt(2.0);
   cout << "x = " << x << endl;
}

#endif
