#include <stdio.h>
#define INFINITY 1000000000
FILE *fp;

void skip_lines(nskip)
int nskip;
{
    int i;

    while ( nskip > 0 && (i = getc(fp)) != EOF )
	if (i == '\n') nskip--;
}



/*
   Print out a range of lines from a file or stdin.
*/

main(argc, argv)
int argc;
char *argv[];
{
    int len = 0, nums = 0;
    int i, iarg = -1, istart, ifinish, nlines;
    int iloc[3];

    if (argc < 2) exit(1);

    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
		case 'l':
		case 'L':
		    len = 1;
		    break;
		case 'n':
		case 'N':
		    nums = 1;
		    break;
		default:
		    if ( (argv[i][1] >= '0' && argv[i][1] <= '9')
			 || argv[i][2] != '\0')
			iloc[++iarg] = i;
	    }
	} else {
	    iloc[++iarg] = i;
	}
    }

    if (iarg < 0) exit(1);

    if ( (fp = fopen(argv[iloc[0]], "r")) == NULL) {
	printf("Can\'t open %s\n", argv[iloc[0]]);
	exit(1);
    }

    if (len != 0) {
	printf("%d lines\n", count_lines());
	exit(0);
    }

    if (iarg > 0 && argc > iloc[1])
	istart = decode(argv[iloc[1]], 1);
    else {
	istart = 1;
	ifinish = INFINITY;
    }

    if (iarg > 1 && argc > iloc[2]) {
	ifinish = decode(argv[iloc[2]], istart);
	if (ifinish == -INFINITY) ifinish = -ifinish;
    } else
	ifinish = INFINITY;

    if (istart < 0 || ifinish < 0 ) {
	nlines = count_lines();
	if (istart == -INFINITY) {
	    istart = nlines;
	    ifinish = nlines;
	} else {
	    if (istart < 0) istart = nlines + istart;
	    if (ifinish < 0) ifinish = nlines + ifinish;
	}
    }

    skip_lines(istart-1);
    for (nlines = istart; nlines <= ifinish; nlines++)
	print_line(nums*nlines);
}

decode(arg, izero)
char arg[];
int izero;
{
    int ioff = 0, iline;

    if (arg[0] == 'e' || arg[0] == 'E'
	       || arg[0] == 'l' || arg[0] == 'L'
                         || (arg[0] == '-' && arg[1] == '\0') )
	iline = -INFINITY;
    else {
	if (arg[0] == '+' || arg[0] == '#' || arg[0] == '^') ioff = 1;
	iline = atoi(&arg[ioff]) + ioff*(izero-1);
    }
    return iline;
}

char line[500];

count_lines()
{
    char c1 = '\0';
    int i;
    int num = 0;

    while ( (i = getc(fp)) != EOF ) if ( (c1 = i) == '\n') num++;

    if (num > 0 && c1 != '\n') num++;
    rewind(fp);
    return num;
}


print_line(n)
int n;
{
    int i = -1, j;

    while ( (j = getc(fp)) != EOF ) {
	line[++i] = j;
	if (line[i] == '\n') break;
    }
    if (j == EOF && i < 0) exit(0);

    line[i] = '\0';
    if (n > 0)
	printf("%6d: %s\n", n, line);
    else
	printf("%s\n", line);

    return i-1;
}
