
/* mytime.c:  Time a series of (mainly floating-point) vector loops. */

#include <stdio.h>
#include <math.h>

typedef double real;
typedef void (*action)(int);

#define NTEST 10
#define NMAX  30000

/* Work vectors: */

static real v1[NMAX], v2[NMAX], v3[NMAX], v4[NMAX];
static int  iv1[NMAX], iv2[NMAX], iv3[NMAX];
static int  ind1[NMAX], ind2[NMAX];

static real s1, s2;
static real x = 1.0, y = -2.5, z = 3.14;

/*------------------------------------------------------------------------*/

#include <sys/times.h>
#include <unistd.h>

struct tms  buffer;
static long ticks_per_sec = 0;
static real initial_cpu = 0.0;

void cpu_init()
{
  times(&buffer);

  /* Use both system and user time because of ambiguities
     with Linux multiprocessing... */

    initial_cpu = (real) (buffer.tms_utime + buffer.tms_stime);

    ticks_per_sec = sysconf(_SC_CLK_TCK);	/* Clock ticks per second */
}

real cpu_time()
{
    if (!ticks_per_sec) cpu_init();

    times(&buffer);
    return ((real) (buffer.tms_utime + buffer.tms_stime - initial_cpu))
      			/ ticks_per_sec;
}

/*------------------------------------------------------------------------*/

void reset_v1(int l)
{
    int i;
    for (i = 0; i < l; i++)
	v1[i] = i / (i + 1.0);
}

void initialize()
{
    int i;
    real jreal = 58427.0;	/* Random numbers! */

    reset_v1(NMAX);

    for (i = 0; i < NMAX; i++) {

      v2[i] = 1.0 + 0.00001/(i+1); /* Make sure v2 > 0 for sqrt and / tests,
				      and close to 1 for repeated multiplies. */
      v3[i] = v2[i] - v1[i];

      iv1[i] = i + 42;
      iv2[i] = -2*i*i;
      iv2[i] = i*i*i - 1000;

      jreal *= 48828125.0;
      while (jreal >= 2147483648.0) jreal /= 2147483648.0;

      ind1[i] = i;
      ind2[i]  = (int)(NMAX*jreal/2147483648.0);
      if (ind2[i] < 0 || ind2[i] >= NMAX) {
	  printf("error: ind2[%d] = %d\n", i, ind2[i]);
	  exit(1);
      }

      v4[i] = jreal/2147483648.0 - 0.5;

  }

  s1 = v1[NMAX/2];
  s2 = v2[NMAX/2];

  cpu_init();
  printf("sizeof(real) = %d\n", sizeof(real));
}

int dummy(int n)
{
  return n;
}

void time_action(action a, int factor, char *name, int ntmax)
{
  real dtime0, dtime[NTEST], t_start, t_end;
  int i, j, maxtst;

  int vlen[NTEST]  = {0, 3, 10, 30, 100, 300, 1000, 3000, 30000, 300000};
  int count[NTEST] = {10000000, 3000000, 1000000, 300000, 100000,
		      30000, 10000, 3000, 300, 30};

  printf("\n%s\n", name);

  for (i = 0; i < NTEST; i++)
    if (vlen[i] <= ntmax) maxtst = i;

  /* Function call timing. */

  i = 0;
  t_start = cpu_time();
  for (j = 0; j < count[i]; j++) dummy(vlen[i]);
  t_end = cpu_time();
  dtime0 = (t_end - t_start) / count[i];

  for (i = 0; i <= maxtst; i++) {
    real sum;

    t_start = cpu_time();
    for (j = 0; j < count[i]; j++) {
      s1 = 1./s1;
      a(vlen[i]);
    }
    t_end = cpu_time();
    dtime[i] = (t_end - t_start) / count[i] - dtime0;

    printf("vlen = %6d, time = %.3e, time/loop = %.3e, Mflops = %.4f\n",
	   vlen[i], t_end - t_start, dtime[i],
	   vlen[i]*1.0e-6 / dtime[i] * factor);
  }

  printf("function overhead = %.3e/call\n", dtime0);

  fflush(stdout);

}
/*------------------------------------------------------------------------*/

/* Routines to time... */

#define for_all(i)	register int i; for (i = 0; i < n; i++)

void itor(int n)	{for_all(i) v1[i]  = iv2[i];}
void rtoi(int n)	{for_all(i) iv1[i] = v2[i];}
void iadd(int n)	{for_all(i) iv1[i] = iv2[i] + iv3[i];}

void vmove(int n)	{for_all(i) v1[i] = v2[i];}

void ssum1(int n)	{for_all(i) s1 += v1[i];}
void ssum2(int n)	{for_all(i) s1 += v1[i] + v2[i];}
void ssum3(int n)	{for_all(i) s1 += v1[i] * v2[i];}

void vsadd1(int n)	{for_all(i) v1[i] += s1;}
void vsmul1(int n)	{for_all(i) v1[i] *= s1;}
void vsdiv1(int n)	{for_all(i) v1[i] /= s2;}

void vsmul1a(int n)
{
    register int i;
    for (i = 0; i < n; i += 4) {
        v1[i] *= s1;
        v1[i+1] *= s1;
        v1[i+2] *= s1;
        v1[i+3] *= s1;
    }
}

void vsadd2(int n)	{for_all(i) v1[i] = v2[i] + s1;}
void vsmul2(int n)	{for_all(i) v1[i] = v2[i] * s1;}
void vsdiv2(int n)	{for_all(i) v1[i] = v2[i] / s2;}

void vsum1(int n)	{for_all(i) v1[i] += v2[i];}
void vsum2(int n)	{for_all(i) v1[i] = v2[i] + v3[i];}
void vmul1(int n)	{for_all(i) v1[i] *= v2[i];}
void vmul2(int n)	{for_all(i) v1[i] = v2[i] * v3[i];}
void vdiv1(int n)	{for_all(i) v1[i] /= v2[i];}
void vdiv2(int n)	{for_all(i) v1[i] = v3[i] / v2[i];}

void saxpy1(int n)	{for_all(i) v1[i] = s1*v2[i] + s2;}
void saxpy2(int n)	{for_all(i) v1[i] = s1*v2[i] + v3[i];}
void saxpy3(int n)	{for_all(i) v1[i] = v2[i]*v3[i] + v4[i];}

void vsqrt(int n)	{for_all(i) v1[i] = sqrt(v2[i]);}
void vabs(int n)	{for_all(i) v1[i] = abs(v3[i]);}
void vsin(int n)	{for_all(i) v1[i] = sin(v3[i]);}
void vexp(int n)	{for_all(i) v1[i] = exp(v2[i]);}
void vpow(int n)	{for_all(i) v1[i] = pow(v2[i], s2);}

void scatter1(int n)	{for_all(i) v1[ind1[i]] = v2[i];}
void scatter2(int n)	{for_all(i) v1[ind2[i]] = v2[i];}
void gather1(int n)	{for_all(i) v1[i] = v2[ind1[i]];}
void gather2(int n)	{for_all(i) v1[i] = v2[ind2[i]];}

void vif(int n)		{for_all(i) if (v4[i] > 0.0) v1[i] = v4[i];}

void vforce(int n)
{
    real d0, d1, d2, rij2, fac;

    for_all(i) {
	d0 = v2[i] - x;
	d1 = v3[i] - y;
	d2 = v4[i] - z;
	rij2 = d0*d0 + d1*d1 + d2*d2 + 0.1;
	fac = v1[i] / (rij2*sqrt(rij2));
	v1[0] += d0*fac;
	v1[1] += d1*fac;
	v1[2] += d2*fac;
    }
}

/*------------------------------------------------------------------------*/

main(int argc, char *argv[])
{
    int n = NMAX;
    if (argc > 1) n = atoi(argv[1]);

    initialize();

    time_action(itor, 1, "vr = vi", n);
    time_action(rtoi, 1, "vi = vr", n);
    time_action(iadd, 1, "vi1 = vi2 + vi3", n);

    time_action(vmove, 1, "v1 = v2", n);

    s1 = 0.0, time_action(ssum1, 1, "s += v", n);
    s1 = 0.0, time_action(ssum2, 2, "s += v + v", n);
    s1 = 0.0, time_action(ssum3, 2, "s += v * v", n);

    s1 = 1.00000123;
    reset_v1(n);
    time_action(vsadd1, 1, "v += s", n);
    reset_v1(n);
    time_action(vsmul1, 1, "v *= s", n);

    reset_v1(n);
    time_action(vsmul1a, 1, "v *= s (partly unrolled)", n);
    time_action(vsadd2, 1, "v1 = v2 + s", n);
    time_action(vsmul2, 1, "v1 = v2 * s", n);
    time_action(vsdiv2, 1, "v1 = v2 / s", n);

    reset_v1(n);
    time_action(vsum1, 1, "v1 += v2", n);
    time_action(vsum2, 1, "v1 = v2 + v3", n);

    reset_v1(n);
    time_action(vmul1, 1, "v1 *= v2", n);
    time_action(vmul2, 1, "v1 = v2 * v3", n);
    reset_v1(n);
    time_action(vdiv1, 1, "v1 /= v2", n);
    time_action(vdiv2, 1, "v1 = v3 / v2", n);

/*     {int i; for (i = 0; i < n; i+=1000) printf("%d %f\n", i, v1[i]);} */

    time_action(saxpy1, 2, "v1 = s*v2 + s", n);
    time_action(saxpy2, 2, "v1 = s*v2 + v3", n);
    time_action(saxpy3, 2, "v1 = v2*v3 + v4", n);
    time_action(vforce, 18, "v1 = pseudo_force(v2, v3, v4)", n);

    time_action(vsqrt, 1, "v1 = sqrt(v2)", n);
    time_action(vabs, 1, "v1 = abs(v2)", n);
    time_action(vsin, 1, "v1 = sin(v2)", n);
    time_action(vexp, 1, "v1 = exp(v2)", n);
    time_action(vpow, 1, "v1 = pow(v2, s)", n);

    time_action(scatter1, 1, "v1[iv] = v2", n);
    time_action(scatter2, 1, "v1[iv] = v2", n);
    time_action(gather1, 1, "v1 = v2[iv]", n);
    time_action(gather2, 1, "v1 = v2[iv]", n);

    time_action(vif, 1, "if (v2 > 0) v1 = v2", n);
}
