#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Reformat kira output: convert (ascii or binary) to
 * (binary or prettily-indented ascii).  Stuart Levy, July 2001.
 */

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


void skipbytes( int count, FILE *inf ) {
  int c;
  at += count;
  while(count-- > 0 && (c = getc(inf)) != EOF)
    ;
}

main(int argc, char *argv[])
{
    char line[512];
    long long pat, start = -1;
    double systime = -1;
    int postprefix = 0;
    char *prog = argv[0];

    if((argc <= 1 && isatty(0)) || (argc>1 && argv[1][0] == '-' && argv[1][1] != 'a')) {
	fprintf(stderr, "Usage: %s [file.kira] > indexfile\n\
Scans an ASCII or binary kira file and emits an index listing\n\
time, and starting and ending byte offsets for its synchronizing snapshot.\n",
			argv[0]);
	exit(1);
    }
    if(argc > 1 && 0!=strcmp(argv[1], "-")) {
	if(freopen(argv[1], "r", stdin) == NULL) {
	    fprintf(stderr, "%s: %s: cannot open input: ", prog, argv[1]);
	    perror("");
	    fprintf(stderr, "Run %s with no arguments for help.\n", prog);
	    exit(1);
	}
    }
    if(argc > 1)
	printf("# %s\n", argv[1]);

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
	    	printf("%lg %lld %lld\n", systime, start, at);
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
