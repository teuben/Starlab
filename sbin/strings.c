#include <stdio.h>
#include <string.h>

/* Manipulate strings for UNIX (but indices start at ZERO!!).
 *
 *	length(of) string		returns the length of the given string
 *	substr(ing) string i j		returns the substring string[i:j]
 *	strindex string1 string2	returns the location of the first
 *					occurrence of string2 in string1,
 *					or -1 if string2 is not found
 */

main(argc, argv)
int argc;
char **argv;
{
    int i, j, iend;
    int strindex = -1, found;

    if (argc < 2) exit(1);

    if (strncmp(argv[0], "length", 6) == 0)			/* length */
	fprintf(stdout, "%d", strlen(argv[1]));

    else if (strncmp(argv[0], "subst", 5) == 0) {		/* substring */
	if (argc < 3) exit(1);

	if (argc < 4)
	    iend = strlen(argv[1]) - 1;
	else {
	    iend = atoi(argv[3]);
	    if (iend >= (int)strlen(argv[1])) iend = strlen(argv[1]) - 1;
	}

	for (i = atoi(argv[2]); i <= iend; i++)
	    fprintf(stdout, "%c", argv[1][i]);

    } else if (strncmp(argv[0], "strindex", 5) == 0) {		/* strindex */
	if (argc < 3) exit(1);

	for (j = 0; j < (int)strlen(argv[1]); j++) {
	    if (argv[1][j] == argv[2][0]) {
	        found = 1;
	        for (i = j+1; i < j + (int)strlen(argv[2]); i++)
		    if (argv[1][i] != argv[2][i-j]) {
		        found = 0;
		        break;
		    }
	        if (found) {
		    strindex = j;
		    break;
	        }
	    }
	}
	fprintf(stdout, "%d", strindex);
    }
}
