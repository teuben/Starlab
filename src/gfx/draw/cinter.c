
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

/*      CINTER: Interface routines to avoid use of FORTRAN "malloc" */

/*
 *	CREBIN: Intermediate storage allocation routine for rebinning.
 */

#define NULL 0

void crebin(string,array,nmax,narr,status,prompt,length)
char *string;
float *array;
int *nmax, *narr, *status, *prompt;
long length;
{
    float *temp;
    int *itemp;

    /* Allocate workspace. */

    if ( (temp = (float *) malloc((*nmax)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for TEMP...\n");
	*status = 1;
	return;
    }

    if ( (itemp = (int *) malloc((*nmax)*sizeof(int))) == NULL) {
	printf("Can't allocate workspace for ITEMP...\n");
	*status = 2;
	return;
    }

    /* Call the FORTRAN routine */

    *status = 0;

#ifdef FORTRAN_TRAILING_UNDERSCORE
    rebin_(string,array,nmax,narr,status,prompt,temp,itemp,length);
#else
    rebin(string,array,nmax,narr,status,prompt,temp,itemp,length);
#endif

    /* Free workspace. */

    free(temp);
    free(itemp);
}

void crebin_(string,array,nmax,narr,status,prompt,length)
char *string;
float *array;
int *nmax, *narr, *status, *prompt;
long length;
{
    crebin(string,array,nmax,narr,status,prompt,length);
}

/*
 *	CAUTOCORREL: Intermediate storage allocation routine for autocorrel.
 */

void cautocorrel(array,nmax)
float *array;
int *nmax;
{
    float *temp;

    /* Allocate workspace. */

    if ( (temp = (float *) malloc((1+(*nmax)/2)*sizeof(float))) == NULL) {
	printf("Can't allocate workspace for TEMP...\n");
	return;
    }

    /* Call the FORTRAN routine */

#ifdef FORTRAN_TRAILING_UNDERSCORE
    autocorrel_(array, nmax, temp);
#else
    autocorrel(array, nmax, temp);
#endif

    /* Free workspace. */

    free(temp);
}

void cautocorrel_(array,nmax)
float *array;
int *nmax;
{
    cautocorrel(array,nmax);
}
