// Program.C   -   command line user interface class
//
// History:
//  0.0 16-aug-96  Initial release to get us going.                         PJT
//                 Lot's of missing features,
//                 such as the options mode (only key=val mode works) and
//                 not much system keyword support. 
//                 The code has a lot of memory leaks.
//                 Also still missing a powerful String class. 
//  0.1 10-apr-97  doc fixes, and initial integration into Starlab          PJT
//  0.2  5-jun-97  enabled minimum match
//
//

#include "program.h"
#include <string.h>
#include <stdlib.h>

#define PROGRAM_VERSION "Program:: Version 0.1 5-jun-97 PJT"

#define MINMATCH


// local functions, currently defined near the end of this file.
//      they are all old NEMO routines and exist because of a lack
//      of a real string class

static int    xstrlen(void *xspt,int nbyt);
static string parname(string arg);
static string parvalue(string arg);
static string parhelp(string arg);

Program *Program::self = 0;               // catch those who don't call init


Program::Program() {
    package = 0;
    progname = 0;
}

Program::Program(string p) {
    package = p;
    progname = 0;
}

Program::~Program() {
    if (debug_level>9) {
        cerr << "[DEBUG " << debug_level << "] dtor " << progname << endl;
    }
}


void Program::init(int argc, string *argv) {
    char *sp = strrchr(argv[0],'/');
    if (sp == 0)            // note progname points to 
        progname = argv[0];
    else
        progname = sp+1;

    Program::argc = argc;
    Program::argv = argv;
    self = this; 

    debug_level = 0;        // should also read from $STARLAB_DEBUG
    error_level = 0;        // should also read from $STARLAB_ERROR
    help_level = 0;         // should also read from $STARLAB_HELP

    nkeys = xstrlen(keywords, sizeof(string));		// fixed size array!!!!
    keys = new keyword[nkeys];
    addkey(0,"  argv0=\n Program name",1);
    keys[0].val = progname;
    for (int i=1; i<nkeys; i++)
        addkey(i,keywords[i-1],0);
}

void Program::append(string skeyvalhelp) {
    Warning("Program::append - not implemented yet");
}

void Program::parse(void) {
#if 0
    int posflag = FALSE;            // force key=val, don't allow positional
#else
    int posflag = TRUE; 
#endif
    string name;

    for (int i=1; argv[i] ; i++) {
        if (debug_level > 0)
            cout << "Parsing " << argv[i] << endl;
        name = parname(argv[i]);
        posflag = posflag && (name == 0);
        if (posflag) {
            if (i >= nkeys) Error("Too many un-named arguments");
            // should free old one here too
            keys[i].val = strdup(argv[i]);
            keys[i].count++;
        } else {
            if (name == 0) {
                cerr << "Parameter " << name << "must be named" << endl;
                exit(1);
            }
            int j = findkey(name);
            if  (j >= 0) {
                if(keys[j].count)
                    Error("Parameter duplicated");
                // free old value here too
                keys[j].val = strdup(parvalue(argv[i]));
                keys[j].count++;
            } else {                // not listed - check if system key
                if (strcmp(name,"help")==0)
                    cout << "SET_HELP" << endl;
                else if (strcmp(name,"debug")==0)
                    set_debug(parvalue(argv[i]));
                else if (strcmp(name,"error")==0)
                    set_error(parvalue(argv[i]));
                else
                    Error("Parameter %s unknown",name);
            } // j>0
        } // if(posflag)
    } // for (i)
}

void Program::fini(void) {
    int i, n=0;

    for (i=1; i<nkeys; i++)
        n += keys[i].upd ? 1 : 0;
    if (n && debug_level > 0) {
        Warning("The following keywords have never been read:");
        for (i=1; i<nkeys; i++)
            if (keys[i].upd) cerr << " " << keys[i].key;
        cerr << endl;
    }
}

int Program::go(void) {

    if (argc > 1) {
	if (strcmp(argv[1],"--help") == 0) {
            cerr << PROGRAM_VERSION << endl;
            cerr << "Current system options:" << endl;
            cerr << " --help         this help" << endl;
            cerr << " --keywords     show program keys, values and help" << endl;
            cerr << " --version      show program version" << endl;
            cerr << " --usage        show program usage line" << endl;
            cerr << " --description  show program description" << endl;
            cerr << " --examples     show program examples" << endl;
            cerr << " --show         show some debugging info" << endl;
            if (argc==2) return 0;
	}
        if (strcmp(argv[1],"--keywords") == 0) {
            for (int i=1; i<nkeys; i++) {
                cout << keys[i].key << "=" << keys[i].val << endl;
                cout << "\t" << keys[i].help << endl;
            }
            if (argc==2) return 0;
        }
        if (strcmp(argv[1],"--version") == 0) {
            cout << progname << " : " << version << endl;
            if (argc==2) return 0;
            return 0;
        }
        if (strcmp(argv[1],"--usage") == 0) {
            cout <<  usage << endl;
            if (argc==2) return 0;
        }
        if (strcmp(argv[1],"--description") == 0) {
            cout << description << endl;
            if (argc==2) return 0;
        }
        if (strcmp(argv[1],"--examples") == 0) {
            cout << examples << endl;
            if (argc==2) return 0;
        }
        if (strcmp(argv[1],"--show") == 0) {
            show();
            if (argc==2) return 0;
        }
    }

    // check if any  '???' values remain


    int missing = 0;
    for (int i=1; i<nkeys; i++)
        if (strcmp(keys[i].val,"???")==0) {
            missing++;
            break;
        }
    if (missing) {
        cerr << "Insufficient parameters" << endl;
        cerr << "Usage: " << progname;
        int otherargs = FALSE;
        for (int i=1; i<nkeys; i++)
            if (strcmp(keys[i].val,"???") == 0)
                printf(" %s=???", keys[i].key);
            else
                otherargs = TRUE;
        cerr << (otherargs ? " ...\n" : "\n");
        Error("The above required keywords have no value");
    }
    return 1;
}

void Program::show(void) { 
	string *sp;
        int i;
        cout << "Program: " << progname << endl;
        cout << "Usage: " << usage << endl;
        cout << "Description: " << endl << description << endl;
	cout << "Version: " << version << endl;
	cout << "Examples: " << examples << endl;
	for (i=0, sp=keywords; *sp; i++, sp++)
	    cout << "key(" << i+1 << ") = " << *sp << endl;
        for (i=1; i<argc; i++)
	    cout << "arg(" << i << ") = " << argv[i] << endl;
        for (i=1; i<nkeys; i++) {
            cout << "key(" << i << ") = " << keys[i].key << endl;
            cout << "def(" << i << ") = " << keys[i].val << endl;
            cout << "help(" << i << ") = " << keys[i].help << endl;
        }
	cout << "debug_level: " << debug_level << endl;
	cout << "error_level: " << error_level << endl;
            
}

Program *Program::main(void) { 
    return self; 
}

int Program::get_argc(void) {
    return argc;
}

char **Program::get_argv(void) {
    return argv;
}





/* private functions of Program:: */

void Program::help(void) {

}

void Program::addkey(int i, string skeyvalhelp, int system=0)
{
    int j;
    if (i >= nkeys) Error("addkey internal error");
    
    keys[i].keyval = skeyvalhelp;
    keys[i].option = *skeyvalhelp;
    keys[i].key = strdup(parname(skeyvalhelp));
    keys[i].val = strdup(parvalue(skeyvalhelp));
    keys[i].help = strdup(parhelp(skeyvalhelp));
    keys[i].count = 0;
    keys[i].upd = 1;
    keys[i].system = system;

    // test internal consistencies, duplicate keys etc.

    for (j=0; j<i; j++) {
        if (keys[j].option == ' ') continue;
        if (keys[j].option == keys[i].option)
            Warning("Option %c duplicated (%s,%s)",
                keys[i].option, keys[i].key, keys[j].key);
    }
    
}


string Program::get(string key) {
    int i = findkey(key);
    if (i<0) Error("%s: Unknown keyword",key);
    keys[i].upd = 0;        // mark keyword as read
    return keys[i].val;
}

int Program::count(string key) {
    int i = findkey(key);
    if (i<0) Error("Illegal key");
    return keys[i].count;
}

int Program::get_debug(void) {
    return debug_level;
}

void Program::set_debug(string debug_value) {
    debug_level = atoi(debug_value);
    cerr << "DEBUG level set to: " << debug_level << endl;
}

void Program::set_error(string error_value) {
    error_level = atoi(error_value);
    if (error_level < 0) Error("%d: Cannot set error < 0",error_level);
    cerr << "ERROR level set to: " << error_level << endl;
}

int Program::dec_error(void) {
    if (error_level == 0) return 0;
    // tricky: we return error_level, then decrement it (saves a temp)
    return error_level--;
}


// local helper functions because we don't have a fancy string class (yet)
// these come pretty much straight from NEMO

static int xstrlen(void *xspt, int nbyt) 
{
    int nval, i;
    int lpflg;
    char *cp = (char *) xspt;

    nval = 0;                                   /* init count of values */
    do {                                        /* loop over values */
        nval++;                                 /*   count one more */
        lpflg = FALSE;                          /*   init loop flag */
        for (i = 0; i < nbyt; i++)              /*   loop over bytes */
            if (*cp++ != 0)                     /*     a byte of data? */
                lpflg = TRUE;                   /*       set loop flag */
    } while (lpflg);                            /* until a NULL value */
    return (nval);                              /* return total count */
}


static string parname(string arg)
{
    static char namebuf[64];
    char *ap, *np;

    ap = (char *) arg;

    if (arg[1] == ' ' || arg[1] == ':') {  // old (NEMO) vs. new (STARLAB) way
        ap++;           // skip the short option
        ap++;           // and the whitespace or colon
    }

#if 0
	// should not have to do this, but ok for error / syntax checking
    while (*ap == ' ')          /* skip initial whitespace */
        ap++;
#endif

    np = &namebuf[0];               /* point to spot where to copy to */
    while ((*np = *ap) != 0) {       /* copy until '=' or ' ' */
#if 0
        if (*np == '=' || *np == ' ') {
#else
        if (*np == '=') {
#endif
            *np = 0;
            return ((string) namebuf);
        }
        np++;
        ap++;
    }
#if 0
    namebuf[0] = 0;
    return ((string) namebuf);  /* return pointer to null string */
#else
    return NULL;                /* return NULL (nothing found) */
#endif    
}

/*
 * PARVALUE: extract value from name=value string, skipping initial whitespace
 *           if no value, pointer to 0 (last element) is returned!
 *           if no key, pointer to (local copy of ) first element is returned
 *  ???      BUG:  when HELPVEC is not defined, helpvec also returned
 *           Note: returns unsafe pointer into the input string
 *           
 */

static string parvalue(string arg)
{
#if 0
    permanent char *valbuf = NULL;      /* dynamic array, with (re)allocate */
    permanent int size = 0;
#else
    static char valbuf[256];         /* static array: dangerous */
#endif

    char *ap;

    ap = (char *) arg;
    while (*ap) {
        if (*ap++ == '=') {
            while (*ap && *ap == ' ')      /* skip whitespace after '=' */
                ap++;
            strncpy(valbuf,ap,255);
            valbuf[255] = 0;
            ap = valbuf;
            while (*ap) {
                if (*ap == '\n') {
                    *ap = 0;
                    return valbuf;          /* return patched "value\nhelp" */
                }
                ap++;
            }
            return valbuf;                  /* return unmodified "value" */
        }
    }
#if 1        
    return ((string) ap);	    /* return pointer to null string */
#else
    return NULL;                    /* return nothing found */
#endif    
}
/* 
 * PARHELP: extract help from a defv[] "keyword=value\n help" string
 *          If no help part (or \n), returns zero string.
 *          Note: returns pointer into original string, assumed
 *          to be R/O (might be in text space of code)
 */
static string parhelp(string arg)
{
    char *cp = (char *) arg;

    while (*cp  && *cp != '\n')      /* scan until EOS or BOH */
        cp++;
    if(*cp == '\n')
        cp++;
    while (*cp && (*cp == ' ' || *cp == '\t'))   /* skip white space > BOH */
        cp++;
    return cp;      // no need to copy to buffer, since this is last 
}

/*
 * FINDKEY:  scan valid keywords and return index (>=0) for a keyname
 *           return -1 if none found
 *	     Optionally match can be done in minimum match mode
 */

int Program::findkey(string name)
{
    int i, j, l, count, last;

    //    Warning("Hello world");

    if (nkeys<=0) return -1;                /* quick: no keywords at all */
    for (i = 0; i < nkeys; i++)             /* try exact match */
        if (strcmp(keys[i].key,name) == 0) return i;

#if defined(MINMATCH)
    l = strlen(name);                       /* try minimum match */
    count = 0;                              /* count # matches on */
    for (i=1; i<nkeys; i++) {               /* any of the program keys */
        if (strncmp(keys[i].key, name,l)==0) {
            last = i;
            count++;
        }
    }
    if (count==1) {
      Warning("Resolving partially matched keyword %s= into %s=",
               name, keys[last].key);
        return last;
    } else if (count > 1) {
        Dprintf(0,"### Minimum match failed, found: ");
        for (j=0; j<nkeys; j++)
            if (strncmp(keys[j].key,name,l)==0)
                Dprintf(0,"%s ",keys[j].key);
	Dprintf(0,"\n");
	Error("Ambiguous keyword %s=",name);
    }
#endif

    return -1;          /* if all else fails: return -1 */
}

