
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

#include <stdio.h>
#include <sys/file.h>

ctox(file, x, n, nmax, lf)
char *file;
float *x;
int *n, *nmax;
long int lf;			/* SUN FORTRAN convention! */
{
    int i;
    FILE *fpin, *fopen();	/* file pointer and function declaration */
    char *name, *malloc();
    unsigned char c;

    *n = 0;

    if ( (name = malloc(lf+1)) == NULL) return;

    for (i = 0; i < lf; i++) name[i] = file[i];
    name[lf] = '\0';

    if ( (fpin = fopen(name, "r")) == NULL) return;

    /* Discard the first 8 characters. */

    if ((int)fread(x, sizeof(unsigned char), 8, fpin) < 8) return;

    /* Read up to nmax characters from the file. */

    while (fread(&c, sizeof(unsigned char), 1, fpin) == 1 && *n < *nmax )
	x[(*n)++] = (float) c;

    fclose(fpin);
}

ctox_(file, x, n, nmax, lf)
char *file;
float *x;
int *n, *nmax;
long int lf;			/* SUN FORTRAN convention! */
{
    ctox(file, x, n, nmax, lf);
}
