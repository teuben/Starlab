
//
// error.C
//

#include "stdinc.h"

void print_message(char* line)
{
    if (line)
	cerr << line;
    else
	cerr << " (no text)";
    cerr << endl;
}

void  err_exit(char * line)
{
    cerr << "error: ";
    print_message(line);
    exit(1);
}

void  warning(char * line)
{
    cerr << "warning: ";
    print_message(line);
}

