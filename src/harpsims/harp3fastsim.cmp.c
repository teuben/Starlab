typedef unsigned long size_t;
typedef long fpos_t;
#pragma pack( 8 )
struct _Kee1 { 
    int _cnt;
    unsigned char  *_ptr;
    unsigned char  *_base;
    int _bufsiz;
    short _flag;
    short _file;
    char  *__newbase;
    void  *_lock;
    unsigned char  *_bufendp;
    } ;
typedef struct _Kee1 FILE;
extern struct _Kee1 _iob[];
extern int fread( );
extern int fwrite( );
extern int _flsbuf( );
extern int _filbuf( );
extern int ferror( );
extern int feof( );
extern void clearerr( );
extern int putchar( );
extern int getchar( );
extern int putc( );
extern int getc( );
extern int remove( );
extern int rename( );
extern struct _Kee1  *tmpfile( );
extern char  *tmpnam( );
extern int fclose( );
extern int fflush( );
extern struct _Kee1  *fopen( );
extern struct _Kee1  *freopen( );
extern void setbuf( );
extern int setvbuf( );
extern int fprintf( );
extern int fscanf( );
extern int printf( );
extern int scanf( );
extern int sprintf( );
extern int sscanf( );
struct _Kee2 { 
    char  * *_a0;
    int _offset;
    } ;
typedef struct _Kee2 va_list;
extern int vfprintf( );
extern int vprintf( );
extern int vsprintf( );
extern int fgetc( );
extern char  *fgets( );
extern int fputc( );
extern int fputs( );
extern char  *gets( );
extern int puts( );
extern int ungetc( );
extern int fgetpos( );
extern int fseek( );
extern int fsetpos( );
extern long ftell( );
extern void rewind( );
extern void perror( );
typedef signed long ptrdiff_t;
typedef unsigned int wchar_t;
typedef unsigned int wctype_t;
typedef int time_t;
typedef int clock_t;
typedef long ssize_t;
typedef unsigned char uchar_t;
typedef unsigned short ushort_t;
typedef unsigned int uint_t;
typedef unsigned long ulong_t;
typedef volatile unsigned char vuchar_t;
typedef volatile unsigned short vushort_t;
typedef volatile unsigned int vuint_t;
typedef volatile unsigned long vulong_t;
struct _Kee3 { 
    long r[1];
    } ;
typedef struct _Kee3  *physadr_t;
struct label_t { 
    long val[10];
    } ;
typedef struct label_t label_t;
typedef int level_t;
typedef int daddr_t;
typedef char  *caddr_t;
typedef long  *qaddr_t;
typedef char  *addr_t;
typedef unsigned int ino_t;
typedef short cnt_t;
typedef int dev_t;
typedef int chan_t;
typedef long off_t;
typedef unsigned long rlim_t;
typedef int paddr_t;
typedef unsigned short nlink_t;
typedef int key_t;
typedef unsigned int mode_t;
typedef unsigned int uid_t;
typedef unsigned int gid_t;
typedef void  *mid_t;
typedef int pid_t;
typedef char slab_t[12];
typedef unsigned long shmatt_t;
typedef unsigned long msgqnum_t;
typedef unsigned long msglen_t;
typedef unsigned int wint_t;
typedef unsigned long sigset_t;
typedef long timer_t;
typedef void  (*sig_t)( );
typedef int id_t;
typedef unsigned int major_t;
typedef unsigned int minor_t;
typedef unsigned int devs_t;
typedef unsigned int unit_t;
typedef unsigned long vm_offset_t;
typedef unsigned long vm_size_t;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef struct _Kee3  *physadr;
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
typedef volatile unsigned char vu_char;
typedef volatile unsigned short vu_short;
typedef volatile unsigned int vu_int;
typedef volatile unsigned long vu_long;
struct _quad { 
    int val[2];
    } ;
typedef struct _quad quad;
typedef long swblk_t;
typedef unsigned long fixpt_t;
typedef int fd_mask;
struct fd_set { 
    int fds_bits[128];
    } ;
typedef struct fd_set fd_set;
extern void bzero( );
#pragma pack( 0 )
struct timeval ;
int select( );
extern int fileno( );
extern struct _Kee1  *fdopen( );
extern char  *cuserid( );
extern int getopt( );
extern char  *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern char  *ctermid( );
extern int getw( );
extern int pclose( );
extern int putw( );
extern struct _Kee1  *popen( );
extern char  *tempnam( );
extern void setbuffer( );
extern void setlinebuf( );
unsigned int h3wait_( );
int npipes = 3;
double xjmem[1000][3];
double vjmem[1000][3];
double xjpred[1000][3];
double vjpred[1000][3];
double ajmem[1000][3];
double jjmem[1000][3];
double mjmem[1000];
double tjmem[1000];
double xjbuf[50][26][3];
double vjbuf[50][26][3];
double ajbuf[50][26][3];
double jjbuf[50][26][3];
double mjbuf[50][26];
double tjbuf[50][26];
int indexbuf[50][26];
int njbuf[50];
double timem;
double ximem[3][3];
double vimem[3][3];
double hmem[3];
double epsmem[3];
int nnb[3];
int nb[3][1024];
static int debug_level = 0;
void h3open_(  )
    
{
    
}
void h3close_(  )
    
{
    
}
int h3npipe_(  )
    
{
    
    return npipes; 
}
void h3setnboards_( nb )
    int  *nb;
    
{
    
}
int h3getnboards_(  )
    
{
    
    return 1; 
}
unsigned int h3wait_(  )
    
{
    
    return 0; 
}
void h3setti_( ti )
    double  *ti;
    
{
    
    timem = ti[0];
}
int h3jpmax_(  )
    
{
    
    return 26; 
}
void h3setmode_(  )
    
{
    
}
void h3jpdma_indirect_( int  *, int  *, double [][3], double [][3], double [][3], double [][3], double  *, double  *, int  * );
void h3mjpdma_indirect_( int  *, int  *, double [][3], double [][3], double [][3], double [][3], double  *, double  *, int  *, 
    int  * );
void h3mjpdma_start_( int  * );
void h3mjpdma_flush_( );
void h3jpdma_indirect_( nj, hostindex, xj, vj, aj, jj, mj, tj, mode )
    int  *nj;
    int  *hostindex;
    double  (*xj)[3];
    double  (*vj)[3];
    double  (*aj)[3];
    double  (*jj)[3];
    double  *mj;
    double  *tj;
    int  *mode;
    
{
    int zero;
    
    zero = 0;
    h3mjpdma_indirect_( nj, hostindex, xj, vj, aj, jj, mj, tj, mode, &zero );
    h3mjpdma_start_( &zero );
}
static int buff_used = 0;
void h3mjpdma_indirect_( nj, hostindex, xj, vj, aj, jj, mj, tj, mode, buff_id )
    int  *nj;
    int  *hostindex;
    double  (*xj)[3];
    double  (*vj)[3];
    double  (*aj)[3];
    double  (*jj)[3];
    double  *mj;
    double  *tj;
    int  *mode;
    int  *buff_id;
    
{
    register int i;
    int bufno;
    int k;
    int hostid;
    long _Kii1;
    
    if (debug_level) {
        fprintf( (_iob + 2), "(h3mjpdma) called %d %d\n", (nj[0]), (mode[0]) );
    } 
    if (nj[0] > 26) {
        fprintf( (_iob + 2), "(h3jpdma) Too large nj (%d > %d\n", nj, 26 );
        exit( ((-1)) );
    } 
    if (buff_id[0] >= 0) {
        bufno = buff_id[0];
    } else {
        bufno = buff_used;
        buff_used++;
    } 
    njbuf[bufno] = nj[0];
    i = 0;
    _Kii1 = ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(vjbuf[bufno][i]) || (void  *)(xjbuf[bufno][i]) > (void  *)(vjbuf[bufno][i]
         + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(ajbuf[bufno][i]) || (void  *)(xjbuf[bufno][i]) > (void  *)(ajbuf[
        bufno][i] + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(jjbuf[bufno][i]) || (void  *)(xjbuf[bufno][i]) > (void  *)(
        jjbuf[bufno][i] + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(xj[hostid]) || (void  *)(xjbuf[bufno][i]) > (void  *)
        (xj[hostid] + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(vj[hostid]) || (void  *)(xjbuf[bufno][i]) > (void  *)(vj[
        hostid] + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(aj[hostid]) || (void  *)(xjbuf[bufno][i]) > (void  *)(aj[
        hostid] + 2)) && ((void  *)(xjbuf[bufno][i] + 2) < (void  *)(jj[hostid]) || (void  *)(xjbuf[bufno][i]) > (void  *)(jj[
        hostid] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(ajbuf[bufno][i]) || (void  *)(vjbuf[bufno][i]) > (void  *)(
        ajbuf[bufno][i] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(jjbuf[bufno][i]) || (void  *)(vjbuf[bufno][i]) > (
        void  *)(jjbuf[bufno][i] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(xj[hostid]) || (void  *)(vjbuf[bufno][i]) > 
        (void  *)(xj[hostid] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(vj[hostid]) || (void  *)(vjbuf[bufno][i]) > (
        void  *)(vj[hostid] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(aj[hostid]) || (void  *)(vjbuf[bufno][i]) > (
        void  *)(aj[hostid] + 2)) && ((void  *)(vjbuf[bufno][i] + 2) < (void  *)(jj[hostid]) || (void  *)(vjbuf[bufno][i]) > (
        void  *)(jj[hostid] + 2)) && ((void  *)(ajbuf[bufno][i] + 2) < (void  *)(jjbuf[bufno][i]) || (void  *)(ajbuf[bufno][i]) > 
        (void  *)(jjbuf[bufno][i] + 2)) && ((void  *)(ajbuf[bufno][i] + 2) < (void  *)(xj[hostid]) || (void  *)(ajbuf[bufno][i])
         > (void  *)(xj[hostid] + 2)) && ((void  *)(ajbuf[bufno][i] + 2) < (void  *)(vj[hostid]) || (void  *)(ajbuf[bufno][i]) > (
        void  *)(vj[hostid] + 2)) && ((void  *)(ajbuf[bufno][i] + 2) < (void  *)(aj[hostid]) || (void  *)(ajbuf[bufno][i]) > (
        void  *)(aj[hostid] + 2)) && ((void  *)(ajbuf[bufno][i] + 2) < (void  *)(jj[hostid]) || (void  *)(ajbuf[bufno][i]) > (
        void  *)(jj[hostid] + 2)) && ((void  *)(jjbuf[bufno][i] + 2) < (void  *)(xj[hostid]) || (void  *)(jjbuf[bufno][i]) > (
        void  *)(xj[hostid] + 2)) && ((void  *)(jjbuf[bufno][i] + 2) < (void  *)(vj[hostid]) || (void  *)(jjbuf[bufno][i]) > (
        void  *)(vj[hostid] + 2)) && ((void  *)(jjbuf[bufno][i] + 2) < (void  *)(aj[hostid]) || (void  *)(jjbuf[bufno][i]) > (
        void  *)(aj[hostid] + 2)) && ((void  *)(jjbuf[bufno][i] + 2) < (void  *)(jj[hostid]) || (void  *)(jjbuf[bufno][i]) > (
        void  *)(jj[hostid] + 2));
    if (!_Kii1) {
        while ( i < nj[0] ) {
            hostid = hostindex[i] - 1;
            xjbuf[bufno][i][0] = xj[hostid][0];
            vjbuf[bufno][i][0] = vj[hostid][0];
            ajbuf[bufno][i][0] = aj[hostid][0];
            jjbuf[bufno][i][0] = jj[hostid][0];
            xjbuf[bufno][i][1] = xj[hostid][1];
            vjbuf[bufno][i][1] = vj[hostid][1];
            ajbuf[bufno][i][1] = aj[hostid][1];
            jjbuf[bufno][i][1] = jj[hostid][1];
            xjbuf[bufno][i][2] = xj[hostid][2];
            vjbuf[bufno][i][2] = vj[hostid][2];
            ajbuf[bufno][i][2] = aj[hostid][2];
            jjbuf[bufno][i][2] = jj[hostid][2];
            tjbuf[bufno][i] = tj[hostid];
            mjbuf[bufno][i] = mj[hostid];
            indexbuf[bufno][i] = hostid;
            i++;
            }
    } else {
        while ( i < nj[0] ) {
            hostid = hostindex[i] - 1;
            xjbuf[bufno][i][0] = xj[hostid][0];
            xjbuf[bufno][i][1] = xj[hostid][1];
            xjbuf[bufno][i][2] = xj[hostid][2];
            vjbuf[bufno][i][0] = vj[hostid][0];
            vjbuf[bufno][i][1] = vj[hostid][1];
            vjbuf[bufno][i][2] = vj[hostid][2];
            ajbuf[bufno][i][0] = aj[hostid][0];
            ajbuf[bufno][i][1] = aj[hostid][1];
            ajbuf[bufno][i][2] = aj[hostid][2];
            jjbuf[bufno][i][0] = jj[hostid][0];
            jjbuf[bufno][i][1] = jj[hostid][1];
            jjbuf[bufno][i][2] = jj[hostid][2];
            tjbuf[bufno][i] = tj[hostid];
            mjbuf[bufno][i] = mj[hostid];
            indexbuf[bufno][i] = hostid;
            i++;
            }
    } 
}
void h3mjpdma_start_( buff_id )
    int  *buff_id;
    
{
    register int i;
    register int k;
    register int j;
    register int buff;
    int _Kii2;
    int _Kii1;
    
    if (debug_level > 0) {
        fprintf( (_iob + 2), "(h3mjpdma_start) called %d \n", (buff_id[0]) );
    } 
    buff = buff_id[0];
    for ( i = 0; i<njbuf[buff]; i++ ) {
        _Kii1 = indexbuf[buff][i];
        if (debug_level > 2) {
            printf( "jpdma_start %d %d %d %f\n", i, _Kii1, buff, xjbuf[buff][i][0] );
        } 
        xjmem[_Kii1][0] = xjbuf[buff][i][0];
        xjmem[_Kii1][1] = xjbuf[buff][i][1];
        xjmem[_Kii1][2] = xjbuf[buff][i][2];
        vjmem[_Kii1][0] = vjbuf[buff][i][0];
        vjmem[_Kii1][1] = vjbuf[buff][i][1];
        vjmem[_Kii1][2] = vjbuf[buff][i][2];
        ajmem[_Kii1][0] = ajbuf[buff][i][0];
        ajmem[_Kii1][1] = ajbuf[buff][i][1];
        ajmem[_Kii1][2] = ajbuf[buff][i][2];
        jjmem[_Kii1][0] = jjbuf[buff][i][0];
        jjmem[_Kii1][1] = jjbuf[buff][i][1];
        jjmem[_Kii1][2] = jjbuf[buff][i][2];
        tjmem[_Kii1] = tjbuf[buff][i];
        mjmem[_Kii1] = mjbuf[buff][i];
    }
}
void h3jpdma_flush_(  )
    
{
    int i;
    
    i = 0;
    while ( i < buff_used ) {
        h3mjpdma_start_( &i );
        i++;
    }
    buff_used = 0;
}
void h3jpdma_reset_(  )
    
{
    
    buff_used = 0;
}
static int njset;
int h3calc_firsthalf_( nj, ni, xi, vi, eps2, h2 )
    int  *nj;
    int  *ni;
    double  (*xi)[3];
    double  (*vi)[3];
    double  *eps2;
    double  *h2;
    
{
    int i;
    int j;
    int k;
    double s;
    double s1;
    double s2;
    long _Kii1;
    long _Kii2;
    
    njset = nj[0];
    if (ni[0] > 3) {
        fprintf( (_iob + 2), "(h3ipdma) Too large ni (%d > %d\n", ni, 3 );
        exit( ((-1)) );
    } 
    i = 0;
    _Kii1 = (void  *)(ximem[i] + 2) < (void  *)(xi[i]) || (void  *)(ximem[i]) > (void  *)(xi[i] + 2);
    _Kii2 = (void  *)(vimem[i] + 2) < (void  *)(vi[i]) || (void  *)(vimem[i]) > (void  *)(vi[i] + 2);
    if (!_Kii1) {
        if (!_Kii2) {
            while ( i < ni[0] ) {
                ximem[i][0] = xi[i][0];
                ximem[i][1] = xi[i][1];
                ximem[i][2] = xi[i][2];
                vimem[i][0] = vi[i][0];
                vimem[i][1] = vi[i][1];
                vimem[i][2] = vi[i][2];
                hmem[i] = h2[i] * 2;
                epsmem[i] = eps2[i] * 2;
                i++;
                }
        } else {
            while ( i < ni[0] ) {
                ximem[i][0] = xi[i][0];
                ximem[i][1] = xi[i][1];
                ximem[i][2] = xi[i][2];
                vimem[i][0] = vi[i][0];
                vimem[i][1] = vi[i][1];
                vimem[i][2] = vi[i][2];
                hmem[i] = h2[i] * 2;
                epsmem[i] = eps2[i] * 2;
                i++;
                }
        } 
    } else {
        if (!_Kii2) {
            while ( i < ni[0] ) {
                ximem[i][0] = xi[i][0];
                ximem[i][1] = xi[i][1];
                ximem[i][2] = xi[i][2];
                vimem[i][0] = vi[i][0];
                vimem[i][1] = vi[i][1];
                vimem[i][2] = vi[i][2];
                hmem[i] = h2[i] * 2;
                epsmem[i] = eps2[i] * 2;
                i++;
                }
        } else {
            while ( i < ni[0] ) {
                ximem[i][0] = xi[i][0];
                ximem[i][1] = xi[i][1];
                ximem[i][2] = xi[i][2];
                vimem[i][0] = vi[i][0];
                vimem[i][1] = vi[i][1];
                vimem[i][2] = vi[i][2];
                hmem[i] = h2[i] * 2;
                epsmem[i] = eps2[i] * 2;
                i++;
                }
        } 
    } 
    if (debug_level > 0) {
        i = 0;
        while ( i < ni[0] ) {
            if (hmem[i] != hmem[i+1] != 0) {
                fprintf( (_iob + 2), "h3calc_firsthalf, eps/h2 error for %d\n", i );
            } 
            i +=  2;
            }
    } 
    j = 0;
    while ( j < nj[0] ) {
        s = timem - tjmem[j];
        s1 = 1.5 * s;
        s2 = s * 2.e0;
        for ( k = 0; k<=2; k++ ) {
            xjpred[j][k] = xjmem[j][k] + (vjmem[j][k] + (ajmem[j][k] + jjmem[j][k] * s) * s) * s;
            vjpred[j][k] = vjmem[j][k] + (ajmem[j][k] + jjmem[j][k] * s1) * s2;
            if (debug_level > 2) {
                printf( "predict %d %d %22.14g %22.14g\n", j, k, xjpred[j][k], vjpred[j][k] );
            } 
        }
        j++;
        }
}
int calculate_force_on_one_particle( i, acc, jerk, pot )
    int i;
    double  (*acc)[3];
    double  (*jerk)[3];
    double  *pot;
    
{
    register int j;
    register int k;
    double sqrt( );
    register double xi;
    register double yi;
    register double zi;
    register double vxi;
    register double vyi;
    register double vzi;
    register double dx;
    register double dy;
    register double dz;
    register double dvx;
    register double dvy;
    register double dvz;
    register double dr2;
    register double drdv;
    register double dr2i;
    register double dr3i;
    register double ax;
    register double ay;
    register double az;
    register double jx;
    register double jy;
    register double jz;
    register double epsi;
    
    jz = 0.e0;
    jy = 0.e0;
    jx = 0.e0;
    az = 0.e0;
    ay = 0.e0;
    ax = 0.e0;
    nnb[i] = 0;
    pot[i] = 0.e0;
    xi = ximem[i][0];
    yi = ximem[i][1];
    zi = ximem[i][2];
    vxi = vimem[i][0];
    vyi = vimem[i][1];
    vzi = vimem[i][2];
    epsi = epsmem[i];
    j = 0;
    while ( j < njset ) {
        dx = xjpred[j][0] - xi;
        dy = xjpred[j][1] - yi;
        dz = xjpred[j][2] - zi;
        dvx = vjpred[j][0] - vxi;
        dvy = vjpred[j][1] - vyi;
        dvz = vjpred[j][2] - vzi;
        dr2 = epsi + dx * dx + dy * dy + dz * dz;
        drdv = dx * dvx + dy * dvy + dz * dvz;
        if (dr2 < hmem[i]) {
            nb[i][nnb[i]] = j + 1;
            nnb[i]++;
        } 
        if (dr2 > 0.e0) {
            dr2i = 1.e0 / dr2;
            dr3i = mjmem[j] * dr2i * sqrt( dr2i );
            drdv *=  3.e0 * dr2i;
            ax +=  dx * dr3i;
            ay +=  dy * dr3i;
            az +=  dz * dr3i;
            jx +=  (dvx - dx * drdv) * dr3i;
            jy +=  (dvy - dy * drdv) * dr3i;
            jz +=  (dvz - dz * drdv) * dr3i;
        } 
        j++;
    }
    acc[i][0] = ax;
    acc[i][1] = ay;
    acc[i][2] = az;
    jerk[i][0] = jx;
    jerk[i][1] = jy;
    jerk[i][2] = jz;
}
int h3calc_lasthalf_( ni, acc, jerk, pot )
    int  *ni;
    double  (*acc)[3];
    double  (*jerk)[3];
    double  *pot;
    
{
    register int i;
    register int j;
    register int k;
    int _Kii1;
    
    if (ni[0] > 3) {
        fprintf( (_iob + 2), "(h3ipdma) Too large ni (%d > %d\n", ni, 3 );
        exit( ((-1)) );
    } 
    if (debug_level) {
        fprintf( (_iob + 2), "lasthalf ni=%d\n", (ni[0]) );
    } 
    for ( _Kii1 = 0; _Kii1<ni[0]; _Kii1++ ) {
        calculate_force_on_one_particle( _Kii1, acc, jerk, pot );
    }
}
int h3calc_lasthalf_old( ni, acc, jerk, pot )
    int  *ni;
    double  (*acc)[3];
    double  (*jerk)[3];
    double  *pot;
    
{
    register int i;
    register int j;
    register int k;
    double sqrt( );
    register double dx[3];
    register double dv[3];
    register double dr2;
    register double drdv;
    register double dr2i;
    register double dr3i;
    long _Kii1;
    int _Kii2;
    
    if (ni[0] > 3) {
        fprintf( (_iob + 2), "(h3ipdma) Too large ni (%d > %d\n", ni, 3 );
        exit( ((-1)) );
    } 
    if (debug_level) {
        fprintf( (_iob + 2), "lasthalf ni=%d\n", (ni[0]) );
    } 
    i = 0;
    _Kii1 = (void  *)(acc[i] + 2) < (void  *)(jerk[i]) || (void  *)(acc[i]) > (void  *)(jerk[i] + 2);
    _Kii2 = !_Kii1;
    while ( i < ni[0] ) {
        if (_Kii2) {
            acc[i][0] = 0.e0;
            acc[i][1] = 0.e0;
            acc[i][2] = 0.e0;
            jerk[i][0] = 0.e0;
            jerk[i][1] = 0.e0;
            jerk[i][2] = 0.e0;
        } else {
            acc[i][0] = 0.e0;
            acc[i][1] = 0.e0;
            acc[i][2] = 0.e0;
            jerk[i][0] = 0.e0;
            jerk[i][1] = 0.e0;
            jerk[i][2] = 0.e0;
        } 
        nnb[i] = 0;
        pot[i] = 0.e0;
        j = 0;
        while ( j < njset ) {
            dr2 = epsmem[i];
            drdv = 0.e0;
            for ( k = 0; k<=2; k++ ) {
                dx[k] = xjpred[j][k] - ximem[i][k];
                dv[k] = vjpred[j][k] - vimem[i][k];
                dr2 +=  dx[k] * dx[k];
                drdv +=  dx[k] * dv[k];
                if (debug_level > 2) {
                    fprintf( (_iob + 2), "rijk %d %d %d %g %21.14g %21.14g\n", i, j, k, dx[k], xjpred[j][k], ximem[i][k] );
                } 
            }
            if (dr2 < hmem[i]) {
                nb[i][nnb[i]] = j + 1;
                nnb[i]++;
            } 
            if (debug_level > 1) {
                fprintf( (_iob + 2), "i, j, nnb[i], dr2, h = %d %d %d %g %g\n", i, j, nnb[i], dr2, hmem[i] );
            } 
            if (dr2 > 0.e0) {
                dr2i = 1.e0 / dr2;
                dr3i = mjmem[j] * dr2i * sqrt( dr2i );
                drdv *=  3.e0 * dr2i;
                for ( k = 0; k<=2; k++ ) {
                    acc[i][k] +=  dx[k] * dr3i;
                    jerk[i][k] +=  (dv[k] - dx[k] * drdv) * dr3i;
                    if (debug_level > 1) {
                        fprintf( (_iob + 2), "force i,j,k %d %d %d %20.10g \n", i, j, k, acc[i][k] );
                    } 
                }
            } else {
                if (debug_level > 1) {
                    fprintf( (_iob + 2), "force skip i,j %d %d  \n", i, j );
                } 
            } 
            j++;
        }
        i++;
        }
}
int h3calc_( nj, ni, xi, vi, eps2, h2, acc, jerk, pot )
    int  *nj;
    int  *ni;
    double  (*xi)[3];
    double  (*vi)[3];
    double  *eps2;
    double  *h2;
    double  (*acc)[3];
    double  (*jerk)[3];
    double  *pot;
    
{
    
    h3calc_firsthalf_( nj, ni, xi, vi, eps2, h2 );
    h3calc_lasthalf_( ni, acc, jerk, pot );
}
void h3nbread_( nboards )
    int  *nboards;
    
{
    
}
void set_debug_level_( int  *level )
{
    
    debug_level = level[0];
}
void h3setdebuglevel_( int  *level )
{
    
    set_debug_level_( level );
}
int compare( i, j )
    int  *i;
    int  *j;
    
{
    
    return i[0] - j[0]; 
}
int h3nblist_( board, chip, nblist )
    int  *board;
    int  *chip;
    int  *nblist;
    
{
    int i;
    
    i = 0;
    while ( i < nnb[chip[0]] ) {
        nblist[i] = nb[chip[0]][i];
        i++;
        }
    return nnb[chip[0]]; 
}
