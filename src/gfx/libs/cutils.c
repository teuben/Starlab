
/*
 * Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
 * by Steve McMillan, Drexel University, Philadelphia, PA.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the author named above.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 *	CUTILS: C-files to allow FORTRAN to do I/O, communicate
 *		with the operating system, etc.  Note the annoying
 *		trailing underscores in some names...
 */

#include <stdio.h>

/*
 *	MYPUTC, MYGETC:	I/O interfaces.
 */

void myputc(c, l)
char *c;
long l;
{
    putchar(*c);
    fflush(stdout);
}

void myputc_(c, l)
char *c;
long l;
{
    myputc(c, l);
}

void mygetc(c, l)
char *c;
long l;
{
    *c = getchar();
}

void mygetc_(c, l)
char *c;
long l;
{
    mygetc(c, l);
}

/*----------------------------------------------------------------------*/

/*
 *      The following routines are unnecessary on a Sun, but
 *      are necessary on other systems (e.g. HP, linux).
 *
 *	Define them so long as they don't conflict with a built-in.
 */

/*
 *	CFLUSH:  FORTRAN-accessible call to flush.
 */

void cflush()
{
    fflush(stdout);
}

void cflush_()
{
    cflush();
}

#include <stdlib.h>

/*
 *	CGETENV:  FORTRAN-accessible call to getenv.
 *		  Accepts and returns FORTRAN character strings.
 */

void cgetenv(name, value, ln, lv)
char *name, *value;
long ln, lv;
{
    int i;
    char *c;
    char temp[100];

    /* Translate name from FORTRAN to C */

    for (i = 0; i < ln; i++) temp[i] = name[i];
    temp[ln] = '\0';

    lv = 0;
    value[0] = '\0';

    if ( c = getenv(temp) )
      do {value[lv++] = *c;} while ( *(++c) > '\0' );
}

void cgetenv_(name, value, ln, lv)
char *name, *value;
long ln, lv;
{
    cgetenv(name, value, ln, lv);
}

/*
 *	CSYSTEM:  Temporary patch for the FORTRAN "system" routine.
 *
 */

int csystem(string, length)
char *string;
long int length;
{
    char *str;
    int i;
    
    str = (char*) malloc(length+1);
    for (i = 0; i < length; i++) str[i] = string[i];
    str[length] = '\0';

    return system(str);
}

int csystem_(string, length)
char *string;
long int length;
{
    csystem(string, length);
}

/*----------------------------------------------------------------------*/

/*
 *	UWAIT:	Wait a specified number of microseconds.
 */

uwait(iwait)
int *iwait;
{
#ifdef HAS_USLEEP
    usleep (*iwait);
#endif
}

uwait_(iwait)
int *iwait;
{
    uwait(iwait);
}
