// test2.C:     test all the CLUI options

#include "program.h"

string Program::keywords[] = {
    "nbody=???\n  number of bodies",
    "nrad=10\n    number of radii",       // NOTE: bad, option is duplicated
    "seed=0\n     a seed integer either in the range [1,2147483647],\n"
                    "or 0 which causes a clock related seed to be used",
    "mass=1\n     Total Mass",
    "radius=1\n   Radius",
    "comment=\n   Optional comment",
    "verbose=f\n  Verbose mode",
    "test=f\n     Test fatal error levels",
    "zero=0\n     Testing optionless keyword",
    end_of_keywords,
};

string Program::version = "1.0 6-jun-97 PJT";

string Program::usage = "testing out various CLUI options";

string Program::description = 
 "This program does not have a final name, but simply tests out all\n"
 "the command line user interface (CLUI) options";

void program_main()
{
  string comment;
  int  n, nr, seed, argc;
  string *argv;
  double  m, r, z;
  bool  v_flag, t_flag;
  char  seedlog[128];

  n = getiparam("nbody");
  nr = getiparam("nrad");
  m = getdparam("mass");
  r = getdparam("radius");
  z = getdparam("zero");
  if (hasvalue("comment")) comment = getparam("comment");
  
  seed = getiparam("seed");
    
  v_flag = getflag("verbose");		// undocumented features
  t_flag = getflag("test");

  if (t_flag) {
    cerr << "t_flag is true" << endl;
  }
  
  if (v_flag) {
    argc = get_argc();
    argv = get_argv();
    cout << "argc=" << argc << endl;
    for (int i=0; i<argc; i++)
      cout << "argv[" << i << "] = " << argv[i] << endl;
  }
  
  //  sprintf(seedlog, "       random number generator seed = %d",seed);
  //  cout << "COMMENT:" << seedlog << endl;
  cout << "nbody = " << n << endl;
  cout << "nrad = " << nr << endl;
  cout << "m = " << m <<endl;
  cout << "r = " << r <<endl;
  if (comment.size()) cout << "comment: " << comment << endl;
}

