#include <stdio.h>

/* Echo the command line to stderr. */

main(int argc, char *argv[])
{
    int i;
    for (i = 1; i < argc; i++) {
	fprintf(stderr, "%s", argv[i]);
	if (i < argc-1) fprintf(stderr, " ");
	else fprintf(stderr, "\n");
    }
}
