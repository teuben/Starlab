// mkplummer.C:     example of a proposed format

#include "program.h"
// #include <cstring>
#include <cstdio>

string Program::keywords[] = {
    "nbody=???\n  number of bodies in the Nbody system",
    "seed=0\n     an integer either in the range [1,2147483647], or \n"
                   "0 which causes a seed to be chosen as the value of\n"
                   "the UNIX clock in seconds; this guarantees that no two\n"
                   "calls will give the same value for the seed if they\n"
                   "are more than 2.0 seconds apart.",
    "mfrac=0.999\n    Mass cutoff",
    "rfrac=22.8042468\n   Radial cutoff",
    "comment=\n   Optional comments added to output stream",
    "iflag=f\n    undocumented flag",
    "oflag=f\n    undocumented flag",
    end_of_keywords,
};

string Program::version = "1.0 13-aug-96 PH";

string Program::usage = "make a Plummer model, with a spatial or mass cut-off";

string Program::description = 
 "mkplummer builds a nbody system according to a Plummer model, in virial\n"
 "units (M=G=-4E=1, with E the total energy), and finite spatial extent\n"
 "which can be regulated by specifying mfrac or rfrac or using their default\n"
 "values.\n"
 "     litt: S.J. Aarseth, M. Henon and R. Wielen (1974),\n"
 "           Astron. and Astrophys. 37, p. 183.";

void program_main()
{

    string comment;
    int  n, seed, argc;
    string *argv;
    double  mfrac, rfrac;
    bool  i_flag, o_flag;
    char  seedlog[128];

    n = getiparam("nbody");
    if (hasvalue("comment")) comment = getparam("comment");
    mfrac = getdparam("mfrac");
    rfrac = getdparam("rfrac");

    //    seed = srandinter(getiparam("seed"));

    i_flag = getflag("iflag");		// undocumented features
    o_flag = getflag("oflag");		// undocumented features

    // NOTE:
    // b->log_history(argc, argv);
    argc = get_argc();
    argv = get_argv();
    cout << "argc=" << argc << endl;
    for (int i=0; i<argc; i++)
        cout << "argv[" << i << "] = " << argv[i] << endl;


    if (o_flag) cerr << "mkplummer: random seed = " << seed << endl;
    sprintf(seedlog, "       random number generator seed = %d",seed);
    cout << "COMMENT:" << seedlog << endl;
    cout << "nbody = " << n << endl;
    cout << "mfrac = " << mfrac <<endl;
}

