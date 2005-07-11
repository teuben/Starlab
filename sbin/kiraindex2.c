#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Reformat kira output: convert (ascii or binary) to
 * (binary or prettily-indented ascii).
 * Written Stuart Levy, July 2001.
 * Copied and modified by Steve McMillan, July 2005.
 */

/*
 * Would be nice to have options to
 *     (1) print the snapshot as it is read, not just the indices
 *     (2) buffer if last snap is desired -- better, use file only
 *         and seek to starting index when needed.
 */

#define USAGE \
"Usage: %s -f file -n snapshots [all] -t time [0]\n\
Scan an ASCII or binary kira file and emit an index listing\n\
time and starting and ending byte offsets for each full snapshot.\n\
Optionally start listing at the specified time and continue for\n\
the specified number of snapshots\n"

struct shortform {
    char nfields;
    char intag[5];
    char outtag[5];
    char counts[4];
    char as32;
    double v[8];
    int fp, vp;
};

struct shortform s_tmpv = { 4, "tmrv", "tmpv", {1,1,3,3}, 0 };
struct shortform s_TL = { 2, "TL", "TL", {1,1}, 1 };

long long at;
int nesting = 0;

void skipbytes(int count, FILE *inf) {
  int c;
  at += count;
  while(count-- > 0 && (c = getc(inf)) != EOF)
    ;
}

void err_exit(char *prog)
{
    fprintf(stderr, USAGE, prog);
    exit(1);
}

main(int argc, char *argv[])
{
    char line[512];
    long long pat, start = -1;
    double systime = -1;
    int postprefix = 0;
    char *prog = argv[0];

    double time = 0;
    int time_set = 0, nsnap = 1, nsnap_set = 0, count = 0;

    if (argc <= 1) {

	if (isatty(0)) err_exit(argv[0]);

    } else {

	int i;
	for (i = 1; i < argc; i++)
	    if (argv[i][0] == '-') 
		switch (argv[i][1]) {

		case 'f':	if (freopen(argv[++i], "r", stdin) == NULL) {
				    fprintf(stderr,
					    "%s: %s: cannot open input: ",
					    prog, argv[i]);
				    perror("");
				    err_exit(argv[0]);
				}
				/* printf("# %s\n", argv[i]); */
				break;

		case '-':
		case 'h':	err_exit(argv[0]);

		case 'n':	nsnap = atoi(argv[++i]);
				nsnap_set = 1;
				break;
		case 't':	time = atof(argv[++i]);
				time_set = 1;
				break;

		} else {

		    // Interpret an unknown argument as a file name.

		    if (freopen(argv[i], "r", stdin) == NULL) {
			fprintf(stderr,
				"%s: %s: cannot open input: ",
				prog, argv[i]);
			perror("");
			err_exit(argv[0]);
		    }
		}
    }

/*
    fprintf(stderr, "time_set = %d, nsnap_set = %d\n", time_set, nsnap_set);
    fprintf(stderr, "nsnap = %d, time = %f\n", nsnap, time);
*/
    
    if (nsnap <= 0) return 1;

    at = 0LL;
    while(fgets(line, sizeof(line), stdin) != NULL) {

	char *s = line;
	int slen, linelen;
	char *val = strchr(s, '=');
	while(*s == ' ') s++;
	for(slen = 1; s[slen] > ' ' && s[slen] != '='; slen++)
	    ;
	linelen = strlen(line);
	at += linelen;
	if(!memcmp(s, ")P", 2)) {
	    nesting--;
	    if(nesting < 0) nesting = 0;
	    if(nesting == 0 && start >= 0) {

	        if (!time_set || systime >= time) {
		    printf("%lg %lld %lld\n", systime, start, at);
		    count++;
		}
		if (nsnap_set && count >= nsnap) return 0;

		start = -1;
		systime = 0;
	    }

	} else if(!memcmp(s, "(P", 2)) {
	    nesting++;
	    postprefix = 1;
	    pat = at - linelen;

	} else if(slen == 4 && !memcmp(s, "tmpv", 4)) {
	    skipbytes( 8 + 8 + 3*8 + 3*8, stdin );
	} else if(slen == 8 && !memcmp(s, "t64mpv32", 8)) {
	    skipbytes( 8 + 4 + 3*4 + 3*4, stdin );
	} else if(slen == 2 && !memcmp(s, "TL", 2)) {
	    skipbytes( 4 + 4, stdin );

	} else if(slen == 11 && !memcmp(s, "system_time", 11)) {
	    if(val) systime = atof(val+1);

	} else if(slen == 4 && !memcmp(s, "name", 4) && val) {
	    do ++val; while(*val == ' ');
	    if(!memcmp(val, "root", 4)) {
		start = pat;
	    }
	}
    }
    return 0;
}
