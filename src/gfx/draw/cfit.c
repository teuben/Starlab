
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
	CFIT: Intermediate storage allocation routine for SVD fitting.
*/

#define NULL (void*)0

void cfit(x, y, z, n, a, m, poly, nw, chisq, status)
float *x, *y, *z, *a, *chisq;
int *n, *m, *nw, *status;
void (*poly)();
{
    float *sig, *u, *v, *w, *b;
    int i, usize;

    /* First dimension of U must actually be MAX(m, n) */

    usize = (*n > *m ? *n : *m);

    if ( (sig = (float *) malloc((*n)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for SIG...\n");
	*status = 1;
	return;
    }

    if ( (u = (float *) malloc(usize*(*m)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for U...\n");
	*status = 1;
	return;
    }

    if ( (v = (float *) malloc((*m)*(*m)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for V...\n");
	*status = 1;
	return;
    }

    if ( (w = (float *) malloc((*m)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for W...\n");
	*status = 1;
	return;
    }

    if ( (b = (float *) malloc(usize*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for B...\n");
	*status = 1;
	return;
    }

    /* Set up weighting for NR routine. */

    if ( *nw <= 1)
	for (i = 0; i < *n; i++) sig[i] = 1.0;
    else
	for (i = 0; i < *n; i++) sig[i] = 1.0/z[i];

    /* Call the FORTRAN routine */

#ifdef FORTRAN_TRAILING_UNDERSCORE
    svdfit_(x, y, sig, n, a, m, u, v, w, b, &usize, m, chisq, poly);
#else
    svdfit(x, y, sig, n, a, m, u, v, w, b, &usize, m, chisq, poly);
#endif

    free(sig);
    free(u);
    free(v);
    free(w);

    *status = 0;
}

void cfit_(x, y, z, n, a, m, poly, nw, chisq, status)
float *x, *y, *z, *a, *chisq;
int *n, *m, *nw, *status;
void (*poly)();
{
  cfit(x, y, z, n, a, m, poly, nw, chisq, status);
}
