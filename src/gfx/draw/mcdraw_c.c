
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
  Call FORTRAN routine mcdraw, but allocate working space first.
*/

#include <stdio.h>

mcdrawc (size, string, which, len)
int *size, *which;
char *string;
long int len;
{
	char *x, *y, *z, *a, *w, *malloc();

	if ((a = malloc(12*(*size))) == NULL) err_quit(1);
	x = y = z = a;
	y += 4*(*size);
	z += 8*(*size);

	if ((w = malloc(4)) == NULL)   err_quit(2);
/*	For now, w really is only one word long. */

	/* FORTRAN call: */

#ifdef FORTRAN_TRAILING_UNDERSCORE
	mcdraw_(x, y, z, a, w, size, string, which, len);
#else
	mcdraw (x, y, z, a, w, size, string, which, len);
#endif
}

mcdrawc_(size, string, which, len)
int *size, *which;
char *string;
long int len;
{
    mcdrawc(size, string, which, len);
}

err_quit(i)
int i;
{
	printf("Error allocating memory block #%d\n", i);
	exit(0);
}
