//
// "Program" is the command line user interface class.
//
// It currently lives in $STARLAB_PATH/src/std/pjt but will move into
// $STARLAB_PATH/src/std when things are ironed out.
// The code however is already being added to the 'std' library libstd.a
//
//
#include <iostream.h>
#include "stdinc.h"
#include "stdarg.h"

#ifndef _PROGRAM_H
#define _PROGRAM_H


typedef char *string;		// this should really become a class

				// "Keyword" should also become a class

typedef struct keyword {    // holds old all keyword information            
    string keyval;              // pointer to original R/O [not used]
    char option;                // short option
    string key;                 // keyword (same as long option)
    string val;                 // value
    string help;                // help
    int count;                  // update/reference count; 0=original default
    int upd;                    // 0=read 1=unread original 2=unread updated
    int system;                 // system keyword ? (0=no 1=yes)
//  struct keyword *next;	// next one in linked list (not used)
} Keyword;


class Program {			// NOTE: one instance of this class will exist

  private:

    // static strings to be defined by programmer - or 
    // dummy lib inherited, see p_*.C

    static string keywords[];   // array of keywords + defaults + help
    static string usage;        // one line usage
    static string version;      // N.M date author
    static string description;  // longer description
    static string examples;     // example(s) of usage

    //

    string package;             // optional package name (NEMO, STARLAB)
    string progname;            // short program name without path
    int argc;                   // argument count, for completeness
    string *argv;               // 0 terminated array of CL strings

    int nkeys;			// number of keywords + 1
    Keyword *keys;		// 0=not used 1=first keyword etc.

    static Program *self;       // trick to find yourself again in classless

  public:

    Program();                  
    Program(string);            
    ~Program();

    // Initialize and program flow

    void init(int argc, string *argv);
    void append(string skeyvalhelp);
    void parse(void);
    int go(void);
    void fini(void);
    Program *main(void);

    // debugging

    void show(void);            // for debugging

    // getting and setting parameters

    string get(string key);    	// get keyword value
    int count(string key);	// get keyword access count
    int get_debug(void);	// get debug_level
    int dec_error(void);	// decrement and return error_level

    int get_argc(void);
    char **get_argv(void);

  private:

    int debug_level;            //  debug= or $XXX_DEBUG
    int error_level;            //  error= or $XXX_ERROR
    int help_level;             //  help=  or $XXX_HELP

    void set_debug(string);	// set 
    void set_error(string);
    void set_help(string);

    void help(void);
    void addkey(int, string, int);
    int findkey(string);

};

//  getparam.C

extern int    hasvalue(string);
extern int    isdefault(string);

extern string getparam(string);
extern real   getrparam(string);
extern int    getiparam(string);

extern int    getflag(string);

extern int    get_argc(void);
extern char **get_argv(void);

extern void   Error(char * ...);
extern void   Warning(char *, ...);
extern void   Dprintf(int, char *, ...);

#endif
