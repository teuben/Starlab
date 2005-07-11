#include <stdio.h>
#include <stdlib.h>

/* Print out a specific range of bytes from the input file/stream.
 * For use with kiraindex.  In order to preserve a simple input
 * format, the first argument must be a file name ("-" ==> use stdin).
 * Subsequent pairs or arguments are start and finish values.
 */

void skipbytes(long long count, FILE *in) {
  int c;
  while(count-- > 0 && (c = getc(in)) != EOF)
    ;
}

void printbytes(long long count, FILE *in) {
  int c;
  while(count-- > 0 && (c = getc(in)) != EOF)
    putc(c, stdout);
}

main(int argc, char *argv[])
{
    int iarg;
    long long start, finish, current = 0;

    if (argc < 4) return 1;
    if (strcmp(argv[1], "-")) {
	if(freopen(argv[1], "r", stdin) == NULL) {
	    fprintf(stderr, "%s: %s: cannot open input: ", argv[0], argv[1]);
	    exit(1);
	}
    }

    iarg = 2;
    while (iarg < argc - 1) {
	start  = atoll(argv[iarg++]);		/* indices from kiraindex */
	finish = atoll(argv[iarg++]);
	/* fprintf(stderr, "start = %lld, finish = %lld\n", start, finish); */
	skipbytes(start-current, stdin);
	printbytes(finish-start, stdin);
	current = finish;
    }
    return 0;
}
