#include "stdinc.h"

main(int argc, char **argv)
{
    extern char *poptarg;
    int c;
    real s = 0;
    char* param_string = "s:";

    while ((c = pgetopt(argc, argv, param_string,
			"$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case 's':	s = atof(poptarg);

	}

    }

    cerr << "sqrt(" << s << ") = ";
    cerr << sqrt(s);
}
