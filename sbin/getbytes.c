#include <stdio.h>

/* Print out a specific range of bytes from the input file/stream.
 * For use with kiraindex.
 */

void skipbytes(int count, FILE *in) {
  int c;
  while(count-- > 0 && (c = getc(in)) != EOF)
    ;
}

void printbytes(int count, FILE *in) {
  int c;
  while(count-- > 0 && (c = getc(in)) != EOF)
    putc(c, stdout);
}

main(int argc, char *argv[])
{
    long long start, finish;

    if (argc == 3) {
	start = atoi(argv[1]);		/* from kiraindex */
	finish = atoi(argv[2]);
    } else if (argc == 4) {
	if (freopen(argv[1], "r", stdin) == NULL) return 1;
	start = atoi(argv[2]);
	finish = atoi(argv[3]);
    } else
	return 1;

    skipbytes(start, stdin);
    printbytes(finish-start, stdin);
}
