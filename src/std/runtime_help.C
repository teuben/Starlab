#include "stdinc.h"

//// runtime_help:  test Starlab run-time help function
////
//// Options:       -help     print this message

#ifndef TOOLBOX

local char* stredit(char* s, char c1, char c2)	// duplicate in kira_init.C
{
    if (!s) return NULL;

    char* s1 = new char[strlen(s)+1];
    if (!s1) return NULL;

    strcpy(s1, s);

    char* ss = s1;
    while (*ss != '\0') {
	if (*ss == c1) *ss = c2;
	ss++;
    }

    return s1;
}

void get_runtime_help(char* source_file, char* compile_date, int level)
{
    // Extract help information from the specified source file.
    // source_file was the location of the file in question at the
    // time the program was compiled, but we should allow for the
    // possibilities that the name has changed or the file has
    // become inaccessible.

    // As of Sep 2001, source_file has the STARLAB_PATH environment
    // removed.

    cerr << endl
	 << "    Starlab version " << STARLAB_VERSION << endl;

    char* s = stredit(compile_date, '_', ' ');
    if (s) {
	cerr << "    program created " << s << endl;
	delete s;
    }

    // Look for the source file.  First try the name verbatim, then
    // try prepending the current STARLAB_PATH environment string.

    char src[1024];
    bool exists = false;

    strcpy(src, source_file);
    ifstream file(src);

    if (file) {

	exists = true;
	file.close();

    } else {

	// See if we can construct a file name from the environment.

	if (getenv("STARLAB_PATH")) {
	    sprintf(src, "%s/%s", getenv("STARLAB_PATH"), source_file);
	    ifstream file(src);
	    if (file) {
		exists = true;
		file.close();
	    }
	}
    }

    cerr << "    source file ";
    if (!exists) {
	cerr << "unavailable" << endl;
	exit(0);
    } else
	cerr << "    " << src << endl << endl;

    // First, find and print out level-1 help lines.

    char cmd[1024];
    strcpy(cmd, "grep '^////' ");
    strcat(cmd, src);
    strcat(cmd, " | sed s%////%\"   \"%");

    // Assume that help lines begin with "//// " (level 1)
    // or "//++ " (level 2).

    system(cmd);
    cerr<< endl;

    // Now check for level-2 help.

    if (level > 1) {
	strcpy(cmd, "grep '^//++' ");
	strcat(cmd, src);
	strcat(cmd, " | sed s%//++%\"    + \"%");

	system(cmd);
	cerr << endl;
    }

    exit(0);
}

void check_runtime_help(int argc, char** argv,
			char* source_file, char* compile_date)
{
    int help_level = 0;

    for (int i = 1; i < argc; i++) {
	if (strstr(argv[i], "-help")) help_level++;
	if (strstr(argv[i], "-HELP")) help_level = 2;
    }

    if (help_level) get_runtime_help(source_file, compile_date, help_level);
}

#else

main(int argc, char** argv)
{
    check_help();
}

#endif
