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
int fullform = 0;
int nesting = 0;
char prefix[8];
char *equals = "";

void setprefix() {
  prefix[0] = '\0';
  if(fullform > 0 && nesting > 0) {
    switch(nesting) {
    case 0: break;
    case 1: prefix[0] = ' '; prefix[1] = '\0'; break;
    case 2: prefix[0] = ' '; prefix[1] = ' '; prefix[2] = '\0';  break;
    case 3: prefix[0] = '\t'; prefix[1] = '\0';  break;
    case 4: prefix[0] = '\t'; prefix[1] = ' '; prefix[2] = '\0';  break;
    case 5: prefix[0] = '\t'; prefix[1] = ' '; prefix[2] = ' '; prefix[3] = '\0';  break;
    default: prefix[0] = '\t'; prefix[1] = '\t'; prefix[2] = '\0';  break;
    }
  }
}

void vflush( struct shortform *sp ) {
  int i, k;
  if(sp->fp == 0) return;
  for(i = k = 0; i < sp->fp; k += sp->counts[i], i++) {
    if(sp->as32) {
	at += fprintf(stdout,
		sp->counts[i]==3 ? "%s%c%s%.8g %.8g %.8g\n"
				 : "%s%c%s%.8g\n",
	    prefix, sp->intag[i], equals, sp->v[k], sp->v[k+1], sp->v[k+2]);
    } else {
	at += fprintf(stdout,
		sp->counts[i]==3 ? "%s%c%s%.17lg %.17lg %.17lg\n"
				 : "%s%c%s%.17lg\n",
	    prefix, sp->intag[i], equals, sp->v[k], sp->v[k+1], sp->v[k+2]);
    }
  }
  sp->fp = sp->vp = 0;
}

double getfloat64( FILE *f ) {
    static int one = 1;
    double v;
    if(fread(&v, sizeof(v), 1, f) <= 0) {
	fprintf(stderr, "Error reading double!\n");
	exit(1);
    }
    if(*(char *)&one == 1) {
	unsigned long long lv = *(unsigned long long *)&v;
	lv = (lv>>32) | (lv<<32);
	lv = (lv&0x0000FFFF0000FFFFLL)<<16
	   | (lv>>16)&0x0000FFFF0000FFFFLL;
	lv = (lv&0x00FF00FF00FF00FFLL)<<8
	   | (lv>>8)&0x00FF00FF00FF00FFLL;
	return *(double *)&lv;
    } else {
	return v;
    }
}

float getfloat32( FILE *f ) {
    static int one = 1;
    float v;
    if(fread(&v, sizeof(v), 1, f) <= 0) {
	fprintf(stderr, "Error reading float!\n");
	exit(1);
    }
    if(*(char *)&one == 1) {
	unsigned int l = (*(unsigned int *)&v)>>16 | (*(unsigned int *)&v)<<16;
	l = (l&0x00FF00FF)<<8 | (l>>8)&0x00FF00FF;
	return *(float *)&l;
    } else {
	return v;
    }
}

void vputfull( struct shortform *sp ) {
    int i, k;
    for(i = 0; i < sp->nfields; i++) {
	if(sp->as32) {
	    if(sp->counts[i] == 1)
		at += fprintf(stdout, "%s%c%s%.8g\n",
			prefix, sp->intag[i], equals, getfloat32(stdin));
	    else
		at += fprintf(stdout, "%s%c%s%.8g %.8g %.8g\n",
			prefix, sp->intag[i], equals, getfloat32(stdin),
				getfloat32(stdin), getfloat32(stdin));
	} else {
	    if(sp->counts[i] == 1)
		at += fprintf(stdout, "%s%c%s%.17lg\n",
			prefix, sp->intag[i], equals, getfloat64(stdin));
	    else
		at += fprintf(stdout, "%s%c%s%.17lg %.17lg %.17lg\n",
			prefix, sp->intag[i], equals, getfloat64(stdin),
				getfloat64(stdin), getfloat64(stdin));
	}
    }
}

void putfloat64( FILE *f, double v ) {
  static int one = 1;
  if(*(char *)&one == 1) {
    unsigned long long lv = *(unsigned long long *)&v;
    lv = (lv>>32) | (lv<<32);
    lv = (lv&0x0000FFFF0000FFFFLL)<<16
       | (lv>>16)&0x0000FFFF0000FFFFLL;
    lv = (lv&0x00FF00FF00FF00FFLL)<<8
       | (lv>>8)&0x00FF00FF00FF00FFLL;
    fwrite( (char *)&lv, 8, 1, f );
  } else {
    fwrite( (char *)&v, 8, 1, f );
  }
  at += 8;
}

void putfloat32( FILE *f, float fv ) {
  static int one = 1;
  if(*(char *)&one == 1) {
    unsigned int l = (*(unsigned int *)&fv)>>16 | (*(unsigned int *)&fv)<<16;
    l = (l&0x00FF00FF)<<8
      | (l>>8)&0x00FF00FF;
    fwrite( (char *)&l, 4, 1, f );
  } else {
    fwrite( (char *)&fv, 4, 1, f );
  }
  at += 4;
}

void copybytes( int count, FILE *inf, FILE *outf ) {
  int c;
  while(count-- > 0 && (c = getc(inf)) != EOF)
    putc(c, outf);
  at += count;
}

int vaccum( char *tok, int toklen, char *val, struct shortform *sp ) {
  char *which = strchr(sp->intag, tok[0]);
  int cur, k;

  if(fullform)
    return 0;

  if(which == NULL || val == NULL || toklen != 1) {
    vflush( sp );
    return 0;
  }

  cur = which - sp->intag;
  if(sp->fp != cur) {
    vflush( sp );
    at += 3 + strlen(val);
    putc( tok[0], stdout );
    putc( ' ', stdout );
    putc( '=', stdout );
    fputs( val, stdout );
    return 1;
  }

  sp->v[sp->vp++] = strtod( val, &val );
  if(sp->counts[cur] == 3) {
      sp->v[sp->vp++] = strtod(val, &val);
      sp->v[sp->vp++] = strtod(val, &val);
  }
  if(++cur >= sp->nfields) {
    at += printf("%s%s%s\n", prefix, sp->outtag, equals);
    if(sp->as32) {
	for(cur = 0; cur < sp->vp; cur++)
	    putfloat32( stdout, sp->v[cur] );
    } else {
	for(cur = 0; cur < sp->vp; cur++)
	    putfloat64( stdout, sp->v[cur] );
    }
    sp->fp = sp->vp = 0;
  } else {
    sp->fp = cur;
  }
  return 1;
}

void fixkeyword( char *s ) {
    if(s[0] != '(' && s[0] != ')')
	return;

    if(fullform) {
	int trunc = 0;
	switch(s[1]) {
	case 'P': if(!memcmp(s+1, "Particle", 8)) trunc = 1; break;
	case 'S': if(!memcmp(s+1, "Star", 4)) trunc = 1; break;
	case 'D': if(!memcmp(s+1, "Dynamics", 8)) trunc = 1; break;
	case 'H': if(!memcmp(s+1, "Hydro", 5)) trunc = 1; break;
	}
	if(trunc) {
	    s[2] = '\n';
	    s[3] = '\0';
	}
    } else {
	if(s[2] == '\n') {
	    switch(s[1]) {
	    case 'P': strcpy(s+1, "Particle\n"); break;
	    case 'S': strcpy(s+1, "Star\n"); break;
	    case 'D': strcpy(s+1, "Dynamics\n"); break;
	    case 'H': strcpy(s+1, "Hydro\n"); break;
	    }
	}
    }
}


main(int argc, char *argv[]) {
    char line[512];
    long long pat, start = -1;
    double systime = -1;
    int postprefix = 0;
    char *prog = argv[0];

    if((argc <= 1 && isatty(0)) || (argc>1 && argv[1][0] == '-' && argv[1][1] != 'a')) {
	fprintf(stderr, "Usage: %s [-a] [file.kira] > outfile.kira\n\
Converts an ASCII kira (tdyn) output stream to binary (tmpv, TL) form.\n\
With -a option, converts to (indented) ASCII form instead.\n\
In any case, for each root-level synchronizing snapshot, reports (to stderr)
the time, and starting and ending byte offsets in the output stream,
for easy indexing.\n",
			argv[0]);
	exit(1);
    }
    if(argc > 1 && !strcmp(argv[1], "-a")) {
	fullform = 1;
	argc--, argv++;
    }
    if(argc > 1 && 0!=strcmp(argv[1], "-")) {
	if(freopen(argv[1], "r", stdin) == NULL) {
	    fprintf(stderr, "%s: %s: cannot open input: ", prog, argv[1]);
	    perror("");
	    fprintf(stderr, "Run %s with no arguments for help.\n", prog);
	    exit(1);
	}
    }

    equals = fullform ? " = " : "=";

    at = 0LL;
    while(fgets(line, sizeof(line), stdin) != NULL) {
	char *s = line;
	int slen;
	char *val = strchr(s, '=');
	while(*s == ' ') s++;
	for(slen = 1; s[slen] > ' ' && s[slen] != '='; slen++)
	    ;
	if(val) {
	    while(*++val == ' ')
		;
	}
	if(!memcmp(s, ")P", 2)) {
	    nesting--;
	    if(nesting < 0) nesting = 0;
	    setprefix();
	    if(nesting == 0 && start >= 0) {
		fixkeyword( s );
	    	fprintf(stderr, "%lg %lld %lld\n", systime, start, at + strlen(prefix) + strlen(s));
		start = -1;
		systime = 0;
	    }

	} else if(!memcmp(s, "(P", 2)) {
	    nesting++;
	    postprefix = 1;
	    pat = at;

	} else if(vaccum( s, slen, val, &s_tmpv )) {
	    continue;
	} else if(vaccum( s, slen, val, &s_TL )) {
	    continue;
	} else if(!memcmp(s, "system_time", 11)) {
	    if(val) systime = atof(val);

	} else if(!memcmp(s, "name", 4)) {
	    if(val && memcmp(val, "root", 4) == 0) {
		start = pat;
	    }
	} else if(!memcmp(s, "tmpv", 4)) {
	    if(fullform) {
		vputfull( &s_tmpv );
	    } else {
		at += fprintf(stdout, "tmpv=\n");
		copybytes(4*8, stdin, stdout);
	    }
	    continue;

	} else if(!memcmp(s, "TL", 2)) {
	    if(fullform) {
		vputfull( &s_TL );
	    } else {
		at += fprintf(stdout, "TL=\n");
		copybytes(2*4, stdin, stdout);
	    }
	    continue;
	}

	vflush( &s_tmpv );
       	vflush( &s_TL );

	if(fullform) {
	    if(val) {
		at += fprintf(stdout, "%s%.*s%s%s", prefix, slen, s, equals, val);
	    } else {
		fixkeyword( s );
		at += fprintf(stdout, "%s%s", prefix, s);
	    }
	} else {
	    /* go fast */
	    if(val) {
		int i;
		at += slen + 1 + strlen(val);
		for(i = 0; i < slen; i++)
		    putc(s[i], stdout);
		putc('=', stdout);
		fputs(val, stdout);		/* includes trailing newline */
	    } else {
		fixkeyword( s );
		at += strlen(s);
		fputs(s, stdout);
	    }
	}
	if(postprefix) {
	    setprefix();
	    postprefix = 0;
	}
    }
    return 0;
}
