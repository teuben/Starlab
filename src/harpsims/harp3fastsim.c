/* 
 *
 * HARP3FASTSIM.C : harp3 simulater which performs the operations all using
 * IEEE-754 arithmetics
 *
 */

#include <stdio.h>
#include <sys/types.h>
#ifdef INTERNAL_OUT
#   undef INTERNAL_OUT
#   define  INTERNAL_OUT 1
#else
#   define  INTERNAL_OUT 0
#endif

unsigned int h3wait_();


#define harpreal double

#define NMAX 1000
#define IPMAX 3
#define JPMAX 26
#define BUFMAX  50
#define NDIM 3

#define NNBMAX 1024

int npipes = IPMAX;
double xjmem[NMAX][NDIM];
harpreal vjmem[NMAX][NDIM];
double xjpred[NMAX][NDIM];
harpreal vjpred[NMAX][NDIM];
harpreal ajmem[NMAX][NDIM];
harpreal jjmem[NMAX][NDIM];
harpreal mjmem[NMAX];
double tjmem[NMAX];

double xjbuf[BUFMAX][JPMAX][NDIM];
harpreal vjbuf[BUFMAX][JPMAX][NDIM];
harpreal ajbuf[BUFMAX][JPMAX][NDIM];
harpreal jjbuf[BUFMAX][JPMAX][NDIM];
harpreal mjbuf[BUFMAX][JPMAX];
double tjbuf[BUFMAX][JPMAX];
int indexbuf[BUFMAX][JPMAX];
int njbuf[BUFMAX];
double timem;

double ximem[IPMAX][NDIM];
harpreal vimem[IPMAX][NDIM];
double hmem[IPMAX];
double epsmem[IPMAX];

int nnb[IPMAX];
int nb[IPMAX][NNBMAX];
static int debug_level = 0;
void h3open_()
{
    
}

void h3close_()
{
}

int h3npipe_()
{
    return npipes;
}

void h3setnboards_(nb)
    int * nb;
{
    
}

int h3getnboards_()
{
    return 1;
    /* simulator always return 1 for now... */
}

unsigned int h3wait_()
{
    return 0;
}

void h3setti_(ti)
    double * ti;
{
    timem= *ti;
}

int h3jpmax_()
{

    return JPMAX;
}


#ifdef G77
#define h3jpdma_indirect_ h3jpdma_indirect__
#define h3mjpdma_indirect_ h3mjpdma_indirect__
#define h3mjpdma_start_ h3mjpdma_start__
#define h3mjpdma_flush_ h3mjpdma_flush__
#endif

void h3setmode_()
{}


void h3jpdma_indirect_(int *,int*,double[][3],double[][3],double[][3],
		  double[][3],double*,double* ,int*);

void h3mjpdma_indirect_(int *,int*,double[][3],double[][3],double[][3],
		  double[][3],double*,double* ,int*,int*);
void h3mjpdma_start_(int *);
void h3mjpdma_flush_();

void h3jpdma_indirect_(nj,hostindex,xj,vj,aj,jj,mj,tj,mode)
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
    int zero = 0;
    h3mjpdma_indirect_(nj,hostindex,xj,vj,aj,jj,mj,tj,mode,&zero);
    h3mjpdma_start_(&zero);
}


static int buff_used = 0;
void h3mjpdma_indirect_(nj,hostindex,xj,vj,aj,jj,mj,tj,mode,buff_id)
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
    register int i;
    int bufno;
    if (debug_level) fprintf(stderr,"(h3mjpdma) called %d %d\n",*nj, *mode);
    if(*nj > JPMAX){
	fprintf(stderr,"(h3jpdma) Too large nj (%d > %d\n",nj, JPMAX);
	exit(-1);
    }
    if(*buff_id >= 0){
	bufno = *buff_id;
    }else{
	bufno = buff_used;
	buff_used ++;
    }
    njbuf[bufno] = *nj;
    for(i=0;i<*nj;i++){
	int k;
	int hostid;
	hostid = hostindex[i]-1;
	for(k=0;k<NDIM;k++){
	    xjbuf[bufno][i][k] = xj[hostid][k];
	    vjbuf[bufno][i][k] = vj[hostid][k];
	    ajbuf[bufno][i][k] = aj[hostid][k];
	    jjbuf[bufno][i][k] = jj[hostid][k];
	}
	tjbuf[bufno][i] = tj[hostid];
	mjbuf[bufno][i] = mj[hostid];
	indexbuf[bufno][i] = hostid;
    }
}

void h3mjpdma_start_(buff_id)
    int * buff_id;
{
    register int i,k, j, buff;
    if (debug_level> 0) fprintf(stderr,"(h3mjpdma_start) called %d \n",*buff_id);
    buff = *buff_id;
    for(i=0;i<njbuf[buff];i++){
	int k;
	int j;
	j= indexbuf[buff][i];
	if (debug_level > 2){
            printf("jpdma_start %d %d %d %f\n",
                   i, j, buff, xjbuf[buff][i][0]);
        }
	for(k=0;k<NDIM;k++){
	    xjmem[j][k] = xjbuf[buff][i][k];
	    vjmem[j][k] = vjbuf[buff][i][k];
	    ajmem[j][k] = ajbuf[buff][i][k];
	    jjmem[j][k] = jjbuf[buff][i][k];
	}
	tjmem[j]= tjbuf[buff][i];
	mjmem[j]= mjbuf[buff][i];
    }

}

void h3jpdma_flush_()
{
    int i;
    for(i=0;i<buff_used; i++){
	h3mjpdma_start_(&i);
    }
    buff_used = 0;
}

void h3jpdma_reset_()
{
    buff_used = 0;
}

static int njset;
h3calc_firsthalf_(nj,ni,xi,vi,eps2,h2)
    int * nj;
    int * ni;
    double xi[][3];
    double vi[][3];
    double eps2[];
    double h2[];
{
    int i,j,k;
    njset = *nj;
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }
    for (i=0;i < *ni; i++){
	for(k=0;k<3;k++)ximem[i][k] = xi[i][k];
	for(k=0;k<3;k++)vimem[i][k] = vi[i][k];
	hmem[i] = 2*h2[i];
	epsmem[i] = 2*eps2[i];
	/* factor of two just to mimic HARP chips...*/
    }
    if (debug_level > 0){
	for (i=0; i< *ni; i+= 2){
	    if ((hmem[i] != hmem[i+1]) ||(hmem[i] != hmem[i+1])){
		fprintf(stderr,"h3calc_firsthalf, eps/h2 error for %d\n", i);
	    }
	}
    }
    for (j=0;j< *nj; j++){
	double s, s1, s2;
	s = timem - tjmem[j];
	s1 = 1.5*s;
	s2 = 2.0*s;
	for (k=0;k<3;k++){
            xjpred[j][k]= ((jjmem[j][k]*s + ajmem[j][k])*s + vjmem[j][k])*s
	    + xjmem[j][k];
            vjpred[j][k]= (jjmem[j][k]*s1 + ajmem[j][k])*s2 + vjmem[j][k];
	    if (debug_level > 2){
		printf("predict %d %d %22.14g %22.14g\n",
		       j,k,xjpred[j][k], vjpred[j][k]);
	    }
	}
    }
}

calculate_force_on_one_particle(i,acc,jerk,pot)
    int i;
    double acc[][3];
    double jerk[][3];
    double pot[];
{
    register int j, k;
    double sqrt();
    double xi,yi,zi,vxi,vyi,vzi;
    register double dx,dy,dz,dvx,dvy,dvz,  dr2, drdv, dr2i, dr3i;
    register double ax,ay,az,jx,jy,jz, epsi;
    ax=ay=az=jx=jy=jz=0.0;
    nnb[i] = 0;
    pot[i] = 0.0;
    xi=ximem[i][0];
    yi=ximem[i][1];
    zi=ximem[i][2];
    vxi=vimem[i][0];
    vyi=vimem[i][1];
    vzi=vimem[i][2];
    epsi = epsmem[i];
    for(j=0;j<njset;j++){
	dr2 = epsi;
	drdv = 0.0;
	dx = xjpred[j][0] - xi;
	dy = xjpred[j][1] - yi;
	dz = xjpred[j][2] - zi;
        dvx = vjpred[j][0] - vxi;
        dvy = vjpred[j][1] - vyi;
        dvz = vjpred[j][2] - vzi;
	dr2 = epsi + dx*dx +  dy*dy +  dz*dz ;
	drdv =  dx*dvx +  dy*dvy +  dz*dvz ;
	if (dr2 < hmem[i]){
	    nb[i][nnb[i]] = j;
	    nnb[i]++;
	}
	if (dr2 > 0.0){
	    dr2i = 1.0/dr2;
	    dr3i = mjmem[j]*dr2i*sqrt(dr2i);
	    drdv = 3.*drdv*dr2i;
	    ax +=  dx*dr3i;
	    ay +=  dy*dr3i;
	    az +=  dz*dr3i;
	    jx += (dvx - dx*drdv)*dr3i;
	    jy += (dvy - dy*drdv)*dr3i;
	    jz += (dvz - dz*drdv)*dr3i;
	}
    }	
    acc[i][0] = ax;
    acc[i][1] = ay;
    acc[i][2] = az;
    jerk[i][0] = jx;
    jerk[i][1] = jy;
    jerk[i][2] = jz;
}
h3calc_lasthalf_(ni,acc,jerk,pot)
    int * ni;
    double acc[][3];
    double jerk[][3];
    double pot[];
{
    register int i, j, k;
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }
    if (debug_level) fprintf(stderr,"lasthalf ni=%d\n",*ni);
    for (i=0; i< *ni; i++)calculate_force_on_one_particle(i,acc,jerk,pot);
}
h3calc_lasthalf_old(ni,acc,jerk,pot)
    int * ni;
    double acc[][3];
    double jerk[][3];
    double pot[];
{
    register int i, j, k;
    if(*ni > IPMAX){
	fprintf(stderr,"(h3ipdma) Too large ni (%d > %d\n",ni, IPMAX);
	exit(-1);
    }
    if (debug_level) fprintf(stderr,"lasthalf ni=%d\n",*ni);
    for (i=0; i< *ni; i++){
	double sqrt();
	register double dx[3], dv[3], dr2, drdv, dr2i, dr3i;
	for(k=0;k<3 ; k++){
	    acc[i][k] = 0.0;
	    jerk[i][k] = 0.0;
	}
	nnb[i] = 0;
	pot[i] = 0.0;
	for(j=0;j<njset;j++){
	    dr2 = epsmem[i];
	    drdv = 0.0;
	    for(k=0;k<3 ; k++){
		dx[k] = xjpred[j][k] - ximem[i][k];
		dv[k] = vjpred[j][k] - vimem[i][k];
		dr2 +=  dx[k]*dx[k];
		drdv += dx[k]*dv[k];
		if (debug_level > 2){
		    fprintf(stderr,"rijk %d %d %d %g %21.14g %21.14g\n",i, j, k, dx[k],
		       xjpred[j][k],ximem[i][k]);
		}
	    }
	    if (dr2 < hmem[i]){
		nb[i][nnb[i]] = j;
		nnb[i]++;
	    }
	    if (debug_level > 1){
		fprintf(stderr, "i, j, nnb[i], dr2, h = %d %d %d %g %g\n",
		       i, j, nnb[i], dr2, hmem[i]);
	    }
	    if (dr2 > 0.0){
		dr2i = 1.0/dr2;
		dr3i = mjmem[j]*dr2i*sqrt(dr2i);
		drdv = 3.*drdv*dr2i;
		for(k=0;k<3 ; k++){
		    acc[i][k] +=  dx[k]*dr3i;
		    jerk[i][k] += (dv[k] - dx[k]*drdv)*dr3i;
		    if (debug_level > 1){
			fprintf(stderr,"force i,j,k %d %d %d %20.10g \n",i, j, k, acc[i][k]);
		    }
		}
	    }else{
		if (debug_level > 1){
		    fprintf(stderr,"force skip i,j %d %d  \n",i, j);
		}
	    }
	}
    }	
}

h3calc_(nj,ni,xi,vi,eps2,h2, acc, jerk, pot)
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
    h3calc_firsthalf_(nj,ni,xi,vi,eps2,h2);
    h3calc_lasthalf_(ni,acc,jerk,pot);
}

void h3nbread_(nboards)
    int * nboards;
{

}
#ifdef G77
# define set_debug_level_ set_debug_level__
#endif
void set_debug_level_(int * level)
{
    debug_level = *level;
}

void h3setdebuglevel_(int * level)
{
    set_debug_level_(level);
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

    int i;
    for(i=0; i< nnb[*chip]; i++){
	nblist[i] = nb[*chip][i];
    }
    return nnb[*chip];
}


