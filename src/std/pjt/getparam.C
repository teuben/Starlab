//
// this is a small library of classless access functions to the
// command line interface. The Program class uses a static
// pointer to itself to find itself back here! 
// Each routine MUST call init_main() to double check if Program
// has been loaded.
//
// The alternative is to define your own main() and instantiate a
// Program and do all the work. This set of routines is merely provided
// for those 99% of the programs which suffice with the simple default
// fixed set of program parameters and don't need any dynamic behavior.
// (Note: we have not quite defined how dynamic we allow the CLUI to be)
//
// History:
//	summer 1996  	Created in Amsterdam		Peter Teuben
//	10-apr-97	minor doc fixes			PJT
//	5-jun-97	fixed up Error/Debug/Dprintf	PJT
//
// Todo:
//      - expression parser for getrparam, getiparam
//
#include "program.h"
#include "stdlib.h"
#include "stdarg.h"

static Program dummy, *main = 0;     // Program

static void init_main(void) {        // check if everything initialized is here
    if (main == 0) {
        main = dummy.main();
        if (main == 0) 
            Error("getparam.C: main could not be initialized");
    }
}

int hasvalue(string key) 
{
    init_main();
    string val = main->get(key);
    if (val == 0 || *val == 0) return FALSE;
    return TRUE;
}

int isdefault(string key)
{
    init_main();
    return main->count(key);
}
    
int getflag(string key) 
{
    init_main();
    char *sp = main->get(key);
    if (strchr("tTyY1",*sp)) return 1;
    if (strchr("fFnN0",*sp)) return 0;
    Error("Syntax error for flag %s=%s",key,sp);
    return 0;
}

string getparam(string key)
{
    init_main();
    return main->get(key);
}

int getiparam(string key)
{
    init_main();
    // Note: no syntax checking yet
    return atoi(main->get(key));
}

real getrparam(string key)
{
    init_main();
    // Note: no syntax checking yet
    return atof(main->get(key));
}

int get_argc(void)
{
    init_main();
    return main->get_argc();
}

char **get_argv(void)
{
    init_main();
    return main->get_argv();
}

//  using printf-style output 

void Error(char *fmt ...)
{
    va_list ap;

    init_main();

    fprintf(stderr,"### Fatal error [%s]: ",main->get("argv0"));  
    
    va_start(ap, fmt);              // ap starts with string 'fmt'
    vfprintf(stderr, fmt, ap);      // print out on stderr

    if (fmt[strlen(fmt)-1] != '\n') // be nice if no newline supplied
        fprintf(stderr,"\n");       // and print it anyhow
    fflush(stderr);                 // flush it NOW 
    va_end(ap);                     // end varargs 


    int e = main->dec_error();
    if (e > 0) return;
   
    if (main->get_debug()>5) abort(); // produce coredump if requested
      exit(-1);                       // close the shop and pass this to parent
}


void Warning(char *fmt ...)
{
    va_list ap;

    init_main();

    fprintf(stderr,"### Warning [%s]: ",main->get("argv0"));  
    
    va_start(ap, fmt);              // ap starts with string 'fmt'
    vfprintf(stderr, fmt, ap);      // print out on stderr
    if (fmt[strlen(fmt)-1] != '\n') // be nice if no newline supplied
        fprintf(stderr,"\n");       // and print it anyhow
    fflush(stderr);                 // flush it NOW 

#if 0
    //  using strings only, not very useful
    for (;;) {
        char *p = va_arg(ap, char *);
        if (p == 0) break;
        cerr << p ;
    }
    cerr << endl;
#endif
    va_end(ap);                     // end varargs 

}

void Dprintf(int level, char *fmt ...)
{
    va_list ap;

    init_main();

    if (level > main->get_debug()) return;

    va_start(ap, fmt);              // ap starts with string 'fmt'
    vfprintf(stderr, fmt, ap);      // print out on stderr
    fflush(stderr);                 // flush it NOW 
    va_end(ap);                     // end varargs 

}




