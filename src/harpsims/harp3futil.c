/* 
 *
 * HARP3FUTIL.C : harp3 low-level fortran interface
 *
 DEBUG MEMO

 95/4/9 multi-board on steinberger STILL PROBLEMATIC...
 SYMPTOM: go into infinite loop in wait_until...
 after h3jpdma_indirect (in main force loop)
 -- occures *only* with nbody4
 tried : change wait status in h3calc_lasthalf.
 Wait for nboard*5 times before dma_result and
 also get_result seems to cure the problem up to
 nboards = 4.... Still not okay for nboards = 6

 The above was illusion. at least nboards*5 still
 fails for nboards = 4. Try nboards*10 again

 nboards*10 and nboards = 4 still failes...

 95/4/10 changed host to alexandria (a 3400)

 Seems to work fine...

 95/4/11 change on steinberger : Increase the number of MB calls in
 -- wait_harp3_until_ready
 -- raw_status_of_board

 seems to fix the problem. This must mean very early acces of
 status after command issuing can cause trouble.
 */

#include <stdio.h>
#include <sys/types.h>
#include "harp3.h"
#include "harp3defects.h"
#ifdef INTERNAL_OUT
#   undef INTERNAL_OUT
#   define  INTERNAL_OUT 1
#else
#   define  INTERNAL_OUT 0
#endif

#ifndef NBOARDS
#  define NBOARDS 1
#endif

#define NWAIT 2+4
/*
 Up to 95/11/08, NWAIT used to be
#define NWAIT 2+2
This changed to see the effect on S8 code
*/
/*
   NWAIT = 1+2 caused (?) trouble on nbody4... (95/5/6)
 */
unsigned int * base;

unsigned int h3wait_();
static int nboards;

static int ndefects;

void jpdma_index_initialize();

int kill_dum;
void kill_time(i)
    int i;
{
    int j;
    for(j=0;j<i*3;j++)kill_dum = rand();
}

    static int njold = 0;
    static int njset;

void h3open_()
{
    unsigned int * open_harp3();
    base =open_harp3();
    nboards = get_nboards();
    ndefects = get_number_of_defect_chips();
    jpdma_index_initialize();
    harp3_testrun_();
    set_mode(base,255,15);
    h3wait_();
    /* code added at 95/10/27 */
    njold = 0;
    
}
void h3open_notest_()
{
    unsigned int * open_harp3();
    base =open_harp3();
    nboards = get_nboards();
    ndefects = get_number_of_defect_chips();
    jpdma_index_initialize();
    set_mode(base,255,15);
    h3wait_();
}

void h3close_()
{
    close_harp3();
}

int h3npipe_()
{
    return get_number_of_pipelines(base);
}

void h3setnboards_(nb)
    int * nb;
{
    nboards =  *nb;
}

int h3getnboards_()
{
    return nboards;
}

unsigned int h3wait_()
{
    return wait_harp3_until_idle(base);
}

void h3setti_(ti)
    double * ti;
{
    raw_set_ti(base, ti);
    raw_set_ti(base, ti);
}
void h3randomwrite_(addr, dat)
    int * addr;
    int *dat;
{
    dma_random_write(base,  *addr, *dat);
}

h3setmode_(mode, iboard)
    int *mode;
    int *iboard;
{
    set_mode(base, *mode, *iboard);
}

h3setled_(mode, iboard)
    int *mode;
    int *iboard;
{
    set_led(base, *mode, *iboard);
}

int h3jpmax_()
{
    int jpdma_max_particles();
    return jpdma_max_particles();
}



h3jpdma_(nj,xj,vj,aj,jj,mj,tj,pj)
    int * nj;
    double xj[][3];
    double vj[][3];
    double aj[][3];
    double jj[][3];
    double mj[];
    double tj[];
    int    pj[];
{
#define JPMAX 100    
    int i, k;
    int index[JPMAX];
    register int p;
    if(*nj > JPMAX){
	fprintf(stderr,"(h3jpdma) Too large nj (%d > %d\n",nj, JPMAX);
	exit(-1);
    }
    for (i=0; i< *nj; i++){
#if 0	
	index[i]= pj[i]*3+0x2000000;
#endif
	p = pj[i];
	index[i]=0x2000000*(p%nboards+1) + 3*(p/nboards);
#if INTERNAL_OUT	
	printf("index, i,  mj = 0x%x %d %g\n",index[i],  i, mj[i]);
#endif	
    }
    jpdma_for_array(base, *nj, index,xj,vj,aj,jj,tj,mj);
}

static int jpdma_index_array[2000000];
static int jpdma_index_initialized = 0;
void jpdma_index_initialize()
{
    int i;
    if(jpdma_index_initialized) return;
    jpdma_index_initialized =1;
    for (i = 0;i < 2000000; i++){
        jpdma_index_array[i+1] =
        0x2000000*(i%nboards+1) + 3*(i/nboards);
    }
}
h3jpdma_indirect_(nj,hostindex,xj,vj,aj,jj,mj,tj,mode)
    int * nj;
    int hostindex[];
    double xj[][3];
    double vj[][3];
    double aj[][3];
    double jj[][3];
    double mj[];
    double tj[];
    int *mode;
{
#define JPMAX 100    
    int i, k;
    int index[JPMAX];
    register int p;
    if(*nj > JPMAX){
	fprintf(stderr,"(h3jpdma) Too large nj (%d > %d\n",nj, JPMAX);
	exit(-1);
    }
    for (i=0; i< *nj; i++){
#if 0
	index[i]= hostindex[i]*3+0x2000000-3;
	p = hostindex[i] - 1;
	index[i]=0x2000000*(p%nboards+1) + 3*(p/nboards);
#endif
        index[i]=jpdma_index_array[hostindex[i]];
#if INTERNAL_OUT	
	printf("j, mj = 0x%x %d %g\n", index[i], i, mj[i]);
#endif	
    }
    jpdma_for_array_indirect(base, *nj, hostindex,index,xj,vj,aj,jj,tj,mj,*mode);
}


static int buff_used = 0;
h3mjpdma_indirect_(nj,hostindex,xj,vj,aj,jj,mj,tj,mode,buff_id)
    int * nj;
    int hostindex[];
    double xj[][3];
    double vj[][3];
    double aj[][3];
    double jj[][3];
    double mj[];
    double tj[];
    int *mode;
    int * buff_id;
{
#define JPMAX 100    
    register int i;
    static int index[JPMAX];
    int bufno;
    if(*nj > JPMAX){
	fprintf(stderr,"(h3jpdma) Too large nj (%d > %d\n",nj, JPMAX);
	exit(-1);
    }
    for (i=0; i< *nj; i+=4){
        index[i]=jpdma_index_array[hostindex[i]];
        index[i+1]=jpdma_index_array[hostindex[i+1]];
        index[i+2]=jpdma_index_array[hostindex[i+2]];
        index[i+3]=jpdma_index_array[hostindex[i+3]];
#if INTERNAL_OUT	
	printf("j, mj = 0x%x %d %g\n", index[i], i, mj[i]);
#endif	
    }
    if(*buff_id >= 0){
	bufno = *buff_id;
    }else{
	bufno = buff_used;
	buff_used ++;
    }
    mjpdma_for_array_indirect(base, *nj, hostindex,index,xj,vj,aj,jj,tj,mj,
			      *mode, bufno);
}

h3mjpdma_start_(buff_id)
    int * buff_id;
{
    jpdma_start(base, *buff_id);
}

h3jpdma_flush_()
{
    int i;
    for(i=0;i<buff_used; i++){
	jpdma_start(base, i);
    }
    buff_used = 0;
}

h3jpdma_reset_()
{
    buff_used = 0;
}


h3jpdma_xvtm_(nj,xj,vj,mj,tj,pj)
    int * nj;
    double xj[][3];
    double vj[][3];
    double mj[];
    double tj[];
    int    pj[];
{
#define JPMAX 100    
    int i, k;
    int index[JPMAX];
    register int p;
    if(*nj > JPMAX){
	fprintf(stderr,"(h3jpdma) Too large nj (%d > %d\n",nj, JPMAX);
	exit(-1);
    }
    for (i=0; i< *nj; i++){
/*	index[i]= pj[i]*3+0x2000000;*/
	p = pj[i] - 1;
	index[i]=0x2000000*(p%nboards+1) + 3*(p/nboards);
#if INTERNAL_OUT	
	printf("j, mj = 0x%x %d %g\n", index[i], i, mj[i]);
#endif	
    }
    jpdma_for_array_xvtm(base, *nj, index,xj,vj,tj,mj);
}



void clear_sprious_particles(nj, nb)
    int nj;
    int nb;
{
    static double x[3];
    static double v[3];
    static double a[3];
    static double j[3];
    static double m;
    static double t;
    int j1, jw, k;
    for(k=0;k<3;k++){
	x[k] =1e50;
	v[k] =0.0;
	a[k] =0.0;
	j[k] =0.0;
    }
    m = 0.0;
    t = 0.0;
    if( nj % nb){
	j1 = nb-(nj%nb); /* number of unused locations */
	for(jw = 0; jw < j1; jw ++){
	    int jtmp = 0x2000000*(nb-jw) + 3*(nj/nb);
	    /* Next h3wait is for test of 3700 */
	    h3wait_();
	    jpdma_for_array(base, 1, &jtmp,x, v, a, j, &t, &m);
	}
    }
}

	    
    
h3calc_(nj,ni,xi,vi,eps2,h2,acc,jerk,pot)
    int * nj;
    int * ni;
    double xi[][3];
    double vi[][3];
    double eps2[];
    double h2[];
    double acc[][3];
    double jerk[][3];
    double pot[];
{
#define IPMAX 100
    static ipdma_packet i_particles[IPMAX];
    static result_packet real_result[IPMAX];
    ipdma_packet * ip;
    int i, k;
    if(*nj != njold){
	clear_sprious_particles(*nj, nboards);
	njold = *nj;
	njset = ((*nj+nboards - 1)/nboards)*3;
    }
	    
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }

    ipdma_for_array_with_real_ni(base, *ni, *ni,xi,vi,eps2,h2);
    wait_harp3_until_idle(base);
    dma_result_initialize(base);
#ifndef ONEBUS
#ifndef TWOBUS    
#ifndef THREEBUS    
#ifndef FOURBUS
#ifndef SIXBUS    
    Present code works only for ONEBUS case....
    aho-----------------
#endif    
#endif    
#endif    
#endif    
#endif
#ifdef ONEBUS    
    calc(base, njset, nboards, *ni - 1);
#endif
#ifdef TWOBUS    
    calc(base, njset, nboards, (*ni)/2 - 1);
#endif
#ifdef THREEBUS    
    calc(base, njset, nboards, (*ni + ndefects)/3 - 1);
#endif
#ifdef FOURBUS    
    calc(base, njset, nboards, (*ni + ndefects)/3 - 1);
#endif
    {
	int i;
	int k;
	kill_time(1);
	
    }
    wait_harp3_until_idle(base);
    wait_harp3_until_idle(base);
    for(i = 0; i<nboards*NWAIT; i++){
	wait_harp3_until_idle(base);
	wait_harp3_until_idle(base);
    }
    {
	int i;
	i = 1;
	kill_time(i);
    }
    dma_result(base);
#if 0    
    get_result(base, *ni, real_result);
    for(i=0;i< *ni;i++){
	for(k=0;k<3;k++)acc[i][k] = real_result[i].acc[k];
	for(k=0;k<3;k++)jerk[i][k] = real_result[i].jerk[k];
	pot[i] = real_result[i].phi;
    }
#else
    get_result_for_array(base, *ni, acc,jerk,pot);
#endif    

}

h3calc_firsthalf_(nj,ni,xi,vi,eps2,h2)
    int * nj;
    int * ni;
    double xi[][3];
    double vi[][3];
    double eps2[];
    double h2[];
{
#define IPMAX 100    
    static ipdma_packet i_particles[IPMAX];
    ipdma_packet * ip;
    int i, k;
    if(*nj != njold){
	clear_sprious_particles(*nj, nboards);
	njold = *nj;
	njset = ((*nj+nboards - 1)/nboards)*3;
    }
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }
    ipdma_for_array(base, *ni, xi,vi,eps2,h2);
    wait_harp3_until_idle(base);
    dma_result_initialize(base);

#ifndef ONEBUS
#ifndef TWOBUS    
#ifndef THREEBUS    
#ifndef FOURBUS    
    Present code works only for ONEBUS case....
    aho-----------------
#endif    
#endif    
#endif    
#endif
#ifdef ONEBUS    
    calc(base, njset, nboards, *ni - 1);
#endif
#ifdef TWOBUS    
    calc(base, njset, nboards, (*ni)/2 - 1);
#endif
#ifdef THREEBUS    
    calc(base, njset, nboards, (*ni + ndefects)/3 - 1);
#endif
#ifdef FOURBUS    
    calc(base, njset, nboards, (*ni + ndefects)/3 - 1);
#endif
}

h3calc_lasthalf_(ni,acc,jerk,pot)
    int * ni;
    double acc[][3];
    double jerk[][3];
    double pot[];
{
    static result_packet real_result[IPMAX];
    ipdma_packet * ip;
    int i, k;
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }
    for(i = 0; i<nboards*NWAIT; i++){
	wait_harp3_until_idle(base);
	wait_harp3_until_idle(base);
    }
    dma_result(base);
    for(i = 0; i<NWAIT; i++){
	wait_harp3_until_idle(base);
	wait_harp3_until_idle(base);
    }
    get_result_for_array(base, *ni, acc,jerk,pot);
}

#define REAL double
void accel_by_harp3_separate_(ni,xi,nj,xj,m, a, pot,eps2)
    int *ni;
    REAL xi[][3];
    int *nj;
    REAL xj[][3];
    REAL m[];
    REAL a[][3];
    REAL pot[];
    REAL *eps2;
{
    h3open_();
    calculate_accel_by_harp3_separate_trial_noopen(*ni,xi,*nj,xj,m,
						   a, pot,*eps2,0);
    h3close_();
}


void h3nbread_(nboards)
    int * nboards;
{
    int get_nb_retval,  get_nb_mode;
    int iwait;
    iwait = 1;
#if INTERNAL_OUT
    printf("(h3nbread) nboards = %d\n", *nboards);
#endif    
    get_nb_retval = get_nb_mode = 1;

    while(get_nb_retval){
	nb_read_low(base, 250, *nboards);

	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();

	h3wait_();
	h3wait_();
	h3wait_();

	get_nb_retval = get_nb_low(base,get_nb_mode);

	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	h3wait_();
	get_nb_mode =0;
    }
}

int compare(i,j)
    int *i;
    int *j;
{
    return *i - *j;
}

int h3nblist_(board, chip, nblist)
    int * board;
    int * chip;
    int * nblist;
{
    int construct_nblist(), nnb, compare(), i;
    nnb = construct_nblist_low(*board, *chip,1, nblist);
#if INTERNAL_OUT    
    printf("h3nblist: nnb = %d\n", nnb);
    for(i=0;i<nnb;i++){
	printf("h3nblist: nb[%d]= %d\n", i, nblist[i]);
    }
    printf("sort...\n");
#endif    
    qsort(nblist, nnb, sizeof(int), compare);
#if INTERNAL_OUT    
    printf("h3nblist: nnb = %d\n", nnb);
    for(i=0;i<nnb;i++){
	printf("h3nblist: nb[%d]= %d\n", i, nblist[i]);
    }
#endif    
    while ((nnb > 0) && (nblist[nnb-1] >= njold)) nnb--;
    for(i=0;i<nnb-1;i++){
        /* This part is not complete yet : possible THREE occation */
	while((nblist[i] == nblist[i+1]) && (i<nnb-1)){
	    int j;
#if INTERNAL_OUT	    
            printf("correcting nblist, %d %d %d %d\n",
                   i, nnb, nblist[i], nblist[i+1]);
#endif	    
	    for(j=i+1; j<nnb-1;j++){
		nblist[j] = nblist[j+1];
	    }
	    nnb--;
	}
    }
#if INTERNAL_OUT    
    printf("h3nblist: adjusted nnb = %d\n", nnb);
#endif    
    return nnb;
}


accel_by_harp3_separate_noopen_(ni,xi,nj,xj,m, a, pot,eps2)
    int *ni;
    double xi[][3];
    int *nj;
    double xj[][3];
    double m[];
    double a[][3];
    double pot[];
    double *eps2;
{
    calculate_accel_by_harp3_separate_trial_noopen(*ni,xi,*nj,xj,m, a, pot,*eps2,0); 
}

int dummy_wait_on_host(idummy)
    int * idummy;
{
    (*idummy) ++;
    return *idummy;
}

int dummy_wait_on_host2(idummy)
    int * idummy;
{
    (*idummy) ++;
    return *idummy;
}

