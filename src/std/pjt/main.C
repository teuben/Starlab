// main.C
//
//  By default a user program defines starlab_main
//  The user could define their own main(), but they need to make sure
//  the user interface is added in a consistent fashion
//

#include "program.h"

Program p;                          // define the Program!

extern void starlab_main(void);	    // to be defined by user

int main(int argc, char **argv)     // entry point for C++
{       
    p.init(argc,argv);              // save the CLI of your program
    p.parse();                      // parse the UI
    if (p.go()) starlab_main();     // if ok to go, call your 'main'
    p.fini();                       // cleanup and report trouble
    return 0;                       // always return ok
}
