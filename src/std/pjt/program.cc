// Program.C   -   command line user interface class
//
// History:
//  0.0 16-aug-96  Initial release for starlab - to get us going.           PJT
//                 Lot's of missing features,
//                 such as the options mode (only key=val mode works) and
//                 not much system keyword support. 
//                 The code has a lot of memory leaks.
//                 Also still missing a powerful String class. 
//  0.1 10-apr-97  doc fixes, and initial integration into Starlab          PJT
//  0.2  5-jun-97  enabled minimum match
//  1.0 24-nov-02  for CARMA, and STL                                       PJT
//
//
//  ToDo:  - make arguments and (member) functions 'const' where appropriate
//

#include "program.h"
#include "stdlib.h"

#define PROGRAM_VERSION "Program:: Version 1.0 24-nov-02 PJT"

#define MINMATCH


// local functions, currently defined near the end of this file.
//      they are all old NEMO routines and exist because of a lack
//      of a real string class

static string parname(string arg);
static string parvalue(string arg);
static string parhelp(string arg);

Program *Program::self = 0;               // catch those who don't call init


Program::Program() {
  //    package = NULL;
  //    progname = NULL;
}

Program::Program(string p) {
  //    package = p;
  //    progname = 0;
}

Program::~Program() {
  if (debug_level>9) {
    cerr << "[DEBUG " << debug_level << "] dtor " << progname << endl;
  }
}


//  init:  this  is usually the first place the code arrives at

void Program::init(int argc, char *argv[]) {
  string pn(argv[0]);
  int i, islash;

  islash = pn.find_last_of('/');
  if (islash >= 0)
    progname = pn.substr(islash+1);
  else
    progname = pn;

  Program::argc = argc;
  Program::argv = new string[argc];
  for (i=0; i<argc; i++) Program::argv[i] = argv[i];
  self = this; 

  debug_level = 0;        // should also read from $STARLAB_DEBUG
  error_level = 0;        // should also read from $STARLAB_ERROR
  help_level = 0;         // should also read from $STARLAB_HELP

  for (nkeys=0; keywords[nkeys] != end_of_keywords; nkeys++)
    ;
  nkeys++;                         // one extra for the "argv0"
  keys = new Keyword[nkeys];       
  addkey(0,"argv0=\n Program name",1);
  keys[0].val = progname;
  for (int i=1; i<nkeys; i++)
    addkey(i,keywords[i-1],0);
}

void Program::append(string skeyvalhelp) {
  cerr << "Program::append - not implemented yet" << endl;
}

void Program::parse(void) {
#if 1
  bool posflag = false;            // force key=val, don't allow positional
#else
  bool posflag = true; 
#endif
  string name;

  for (int i=1; i<argc ; i++) {
    if (debug_level > 0)
      cout << "Parsing " << argv[i] << endl;
    name = parname(argv[i]);
    posflag = posflag && (name.size() == 0);
    if (posflag) {
      if (i >= nkeys) {
	cerr << "Too many un-named arguments" << endl;
	exit(1);
      }
      // should free old one here too
      keys[i].val = argv[i];
      keys[i].count++;
    } else {
      if (name.size() == 0) {
	cerr << "Parameter " << name << "must be named" << endl;
	exit(1);
      }
      int j = findkey(name);
      if  (j >= 0) {
	if(keys[j].count) {
	  cerr << "Parameter duplicated" << endl;
	  exit(1);
	}
	// free old value here too
	keys[j].val = parvalue(argv[i]);
	keys[j].count++;
      } else {                // not listed - check if system key
	if (name == "help")
	  cout << "SET_HELP" << endl;
	else if (name == "debug")
	  set_debug(parvalue(argv[i]));
	else if (name == "error")
	  set_error(parvalue(argv[i]));
	else {
	  cerr << "Parameter " << name << "unknown" << endl;
	  exit(1);
	}
      } // j>0
    } // if(posflag)
  } // for (i)
}

void Program::fini(void) {
  int i, n=0;

  for (i=1; i<nkeys; i++)
    n += keys[i].upd ? 1 : 0;
  if (n && debug_level > 0) {
    cerr << "The following keywords have never been read:" ;
    for (i=1; i<nkeys; i++)
      if (keys[i].upd) cerr << " " << keys[i].key;
    cerr << endl;
  }
}

int Program::go(void) {

    if (argc > 1) {
      if (argv[1] == "--help") {
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
        if ( argv[1] == "--keywords") {
            for (int i=1; i<nkeys; i++) {
                cout << keys[i].key << "=" << keys[i].val << endl;
                cout << "\t" << keys[i].help << endl;
            }
            if (argc==2) return 0;
        }
        if (argv[1] == "--version") {
            cout << progname << " : " << version << endl;
            if (argc==2) return 0;
            return 0;
        }
        if (argv[1] == "--usage") {
            cout <<  usage << endl;
            if (argc==2) return 0;
        }
        if (argv[1] == "--description") {
            cout << description << endl;
            if (argc==2) return 0;
        }
        if (argv[1] == "--examples") {
            cout << examples << endl;
            if (argc==2) return 0;
        }
        if (argv[1] == "--show") {
            show();
            if (argc==2) return 0;
        }
    }

    // check if any  '???' values remain


    int missing = 0;
    for (int i=1; i<nkeys; i++)
        if (keys[i].val == "???") {
            missing++;
            break;
        }
    if (missing) {
        cerr << "Insufficient parameters" << endl;
        cerr << "Usage: " << progname;
        bool otherargs = false;
        for (int i=1; i<nkeys; i++)
	    if (keys[i].val == "???")
	      cerr << keys[i].key <<  "=???" ;
            else
	      otherargs = true;
        cerr << (otherargs ? " ...\n" : "\n");
        cerr << "The above required keywords have no value" << endl;
	exit(1);
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
  for (i=0; keywords[i].size(); i++)
    cout << "key(" << i+1 << ") = " << keywords[i] << endl;
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

string *Program::get_argv(void) {
    return argv;
}





/* private functions of Program:: */

void Program::help(void) {

}

void Program::addkey(int i, string skeyvalhelp, int system=0)
{
  int j;
  if (i >= nkeys) {
    cerr << "addkey internal error" << endl;
    exit(1);
  }
  
  keys[i].keyval = skeyvalhelp;
  keys[i].option = skeyvalhelp[0];
  keys[i].key = parname(skeyvalhelp);
  keys[i].val = parvalue(skeyvalhelp);
  keys[i].help = parhelp(skeyvalhelp);
  keys[i].count = 0;
  keys[i].upd = 1;
  keys[i].system = system;

  cerr << "addkey: " << keys[i].key << endl;
  
  // test internal consistencies, duplicate keys etc.
#if 0
  for (j=0; j<i; j++) {
    if (keys[j].option == ' ') continue;
    if (keys[j].option == keys[i].option)
      cerr  << "Option " <<  keys[i].option << " duplicate (" <<
	keys[i].key << "," <<  keys[j].key << endl;
  }
#endif
}


string Program::get(string key) {
    int i = findkey(key);
    if (i<0) {
      cerr << key << " Unknown keyword" << endl;
      exit(1);
    }
    keys[i].upd = 0;        // mark keyword as read
    return keys[i].val;
}

int Program::count(string key) {
    int i = findkey(key);
    if (i<0) {
      cerr <<  key << " Illegal key" << endl;
      exit(1);
    }
    return keys[i].count;
}

int Program::get_debug(void) {
    return debug_level;
}

void Program::set_debug(const string debug_value) {
  const char *dv = debug_value.c_str();
  debug_level = atoi(dv);
  cerr << "DEBUG level set to: " << debug_level << endl;
}

void Program::set_error(const string error_value) {
  const char *ev = error_value.c_str();
  error_level = atoi(ev);
  if (error_level < 0) {
    cerr << error_level << " Cannot set error < 0" << endl;
    exit(1);
  }
  cerr << "ERROR level set to: " << error_level << endl;
}

int Program::dec_error(void) {
    if (error_level == 0) return 0;
    // tricky: we return error_level, then decrement it (saves a temp)
    return error_level--;
}


// local helper functions because we don't have a fancy string class (yet)
// these come pretty much straight from NEMO

static string parname(string arg)
{
  int ieq = arg.find_first_of('=');
  if (ieq < 0)
    cerr << "parname " << arg << " has no equals" << endl;
  return arg.substr(0,ieq);
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
  int ieq = arg.find_first_of('=');
  if (ieq < 0)
    cerr << "parvalue " << arg << " has no equals" << endl;
  return arg.substr(ieq+1);
}

/* 
 * PARHELP: extract help from a defv[] "keyword=value\n help" string
 *          If no help part (or \n), returns zero string.
 *          Note: returns pointer into original string, assumed
 *          to be R/O (might be in text space of code)
 */
static string parhelp(string arg)
{
  int ieq = arg.find_first_of('\n');
  if (ieq < 0)
    cerr << "parhelp " << arg << " has no newline" << endl;
  return arg.substr(ieq+1);
}

/*
 * FINDKEY:  scan valid keywords and return index (>=0) for a keyname
 *           return -1 if none found
 *	     Optionally match can be done in minimum match mode
 */

int Program::findkey(const string name)
{
    int i, j, l, count, last;

    if (nkeys<=0) return -1;                /* quick: no keywords at all */
    for (i = 0; i < nkeys; i++)             /* try exact match */
        if (keys[i].key == name) return i;

#if 0
    // lotsa compile booboos still....
#if defined(MINMATCH)
    l = strlen(name);                       /* try minimum match */
    count = 0;                              /* count # matches on */
    for (i=1; i<nkeys; i++) {               /* any of the program keys */
        if (keys[i].key.compare(name,l)) {
            last = i;
            count++;
        }
    }
    if (count==1) {
      cerr << "Resolving partially matched keyword " <<  name 
	   << " into ", keys[last].key << "=" << endl;
      return last;
    } else if (count > 1) {
      cerr << "### Minimum match failed, found: " ;
        for (j=0; j<nkeys; j++)
            if (strncmp(keys[j].key,name,l)==0)
	      cerr << keys[j].key << " ";
	cerr << endl;
	cerr << "Ambiguous keyword %s=" << name << endl;
	exit(1);
    }
#endif
#endif

    return -1;          /* if all else fails: return -1 */
}

