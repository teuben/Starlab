
/* Version of hydro2 with ximage calls built in... */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double real;
#define local static

/*------------------------------------------------------------------------*/

/* Interface routine between the hydro and the X display.
 * The job of this function is to take the rho data and scale
 * appropriately into a character array, then send the array
 * the X-display routines.
 */

static unsigned char *data = NULL;
static nrep = 0;
static unsigned char red[256], green[256], blue[256];

local void output_grid(int m, int n, real *rho, real rhoi, real amp)
{
  /* Scale and write the rho grid to the X display.
     Note that rho is really a 2-D array of dimensions (m+2) x (n+2). */

  unsigned char c;
  int i, j, k = 0;

  if (!data) {
    nrep = get_nrep();
    get_cmap(red, green, blue);
    data = (unsigned char*)malloc(nrep*m*n);
  }

  for (j = 1; j <= n; j++)
    for (i = 1; i <= m; i++) {

      c = 128 * (1 + (*(rho + i + j*(m+2)) - rhoi) / amp);
      if (c > 255) c = 255;

      if (nrep == 1) {

	*(data + k++) = c;

      } else {

	/* Remap the data:  Turn c into a 3-byte (B, G, R) sequence
	   followed by a null byte. */

	*(data + k++) = blue[c];
	*(data + k++) = green[c];
	*(data + k++) = red[c];
	*(data + k++) = 0;

      }

    }

  ximage(data);

}

/*------------------------------------------------------------------------*/

/* Grid parameters */

#define M	200
#define N	200

#define STANDING_WAVE	1
#define OVERDENSE	1

/* Timestep parameters (smaller courant parameter in 2D) */

#define EPS	0.5
#define TMAX	5.0

/* Equation of state */

#define A	1.0
#define GAMMA	1.666667

/* May be most convenient to work with static data in this case... */

real xmin, xmax, dx, ymin, ymax, dy;
real rhoi = 1.0, amp = 0.1;

/* NOTE: Indices run from 0 to M(N)+1, with 0 and M+1 used to implement
 *	 the B.C.s.  The first "physical" zone is #1, which (in x) runs from
 * 	 xmin to xmin+dx (and similarly for y).  Density and presssure are
 *	 defined at the centers of zones: rho[1] is at x+dx/2, and so on.
 */

real u[M+2][N+2], v[M+2][N+2], rho[M+2][N+2], p[M+2][N+2];
real u0[M+2][N+2], v0[M+2][N+2], rho0[M+2][N+2];

real csound, dtcourant;
real dt, dtdx, dtdy;

inline local real pressure(real rho)
{
  return A*pow(rho, GAMMA);
}

local void initialize_grid()
{
  /* Initial conditions. */

  int i, j;
  real nosc = 2;
  real delx, dely;

  xmin = 0.0;
  xmax = 1.0;
  delx = xmax - xmin;
  dx = delx / M;

  ymin = 0.0;
  ymax = 1.0;
  dely = ymax - ymin;
  dy = dely / M;

  for (i = 1; i <= M; i++) {
    real x = xmin + (i-0.5)*dx;
    for (j = 1; j <= N; j++) {
      real y = ymin + (j-0.5)*dy;

#if STANDING_WAVE

      /* Standing wave. */

      u[i][j] = 0.0;
      v[i][j] = 0.0;
      rho[i][j] = rhoi * (1.0 + amp * sin(2*M_PI*nosc*x/delx)
			  	    * sin(2*M_PI*nosc*y/dely));

#else if OVERDENSE

      /* Overdense region. */

      u[i][j] = 0.0;
      v[i][j] = 0.0;
      if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.01)
	rho[i][j] = rhoi * (1.0 + amp);
      else
	rho[i][j] = rhoi;

#endif

    }
  }

  csound = sqrt(GAMMA*pressure(rhoi)/rhoi);
  dtcourant = dx/csound;

  /* Work at a fraction EPS of the Courant limit. */

  dt = EPS*dtcourant;
  dtdx = 0.5*dt/dx;
  dtdy = 0.5*dt/dy;

}

local void copy_grid()
{
  /* Make copies of u, v and rho for differencing below, and
   * set up the pressure using the equation of state. */

  int i, j;

  for (i = 1; i <= M; i++)
    for (j = 1; j <= N; j++) {
      u0[i][j] = u[i][j];
      v0[i][j] = v[i][j];
      rho0[i][j] = rho[i][j];
      p[i][j] = pressure(rho[i][j]);
    }
}

local void apply_bc()
{
  /* Apply boundary conditions (periodic, here). */

  int i, j;

  /* Top and bottom. */

  for (i = 1; i <= M; i++) {
    u0[i][0] = u[i][N];
    v0[i][0] = v[i][N];
    rho0[i][0] = rho[i][N];
    p[i][0] = p[i][N];

    u0[i][N+1] = u[i][1];
    v0[i][N+1] = v[i][1];
    rho0[i][N+1] = rho[i][1];
    p[i][N+1] = p[i][1];
  }

  /* Left and right edges. */

  for (j = 0; j <= N; j++) {
    u0[0][j] = u[M][j];
    v0[0][j] = v[M][j];
    rho0[0][j] = rho[M][j];
    p[0][j] = p[M][j];

    u0[M+1][j] = u[1][j];
    v0[M+1][j] = v[1][j];
    rho0[M+1][j] = rho[1][j];
    p[M+1][j] = p[1][j];
  }

}

local void single_step()
{
  /* Lax differencing scheme. */

  int i, j;

  for (i = 1; i <= M; i++)
    for (j = 1; j <= N; j++) {

      real dux = u0[i+1][j] - u0[i-1][j];
      real duy = u0[i][j+1] - u0[i][j-1];
      real dvx = v0[i+1][j] - v0[i-1][j];
      real dvy = v0[i][j+1] - v0[i][j-1];
      real dpx = p[i+1][j] - p[i-1][j];
      real dpy = p[i][j+1] - p[i][j-1];

      real drux = rho0[i+1][j]*u0[i+1][j] - rho0[i-1][j]*u0[i-1][j];
      real drvy = rho0[i][j+1]*v0[i][j+1] - rho0[i][j-1]*v0[i][j-1];

      u[i][j] = 0.25*(u0[i][j+1] + u0[i][j-1] + u0[i+1][j] + u0[i-1][j])
		- dtdx*(dpx/rho0[i][j] + u0[i][j]*dux)
		- dtdy*v0[i][j]*duy;

      v[i][j] = 0.25*(v0[i][j+1] + v0[i][j-1] + v0[i+1][j] + v0[i-1][j])
		- dtdx*u0[i][j]*dvx
		- dtdy*(dpy/rho0[i][j] + v0[i][j]*dvy);

      rho[i][j] = 0.25*(rho0[i][j+1] + rho0[i][j-1]
			 + rho0[i+1][j] + rho0[i-1][j])
		  - dtdx*drux - dtdy*drvy;
    }
}


local void step()
{
  /* Take a time step. */

  copy_grid();
  apply_bc();
  single_step();
}

void main(int argc, char *argv[])
{
  int i, nsteps = 200, nout = 2;

  initialize_grid();
  ximage_init(M, N);
  output_grid(M, N, (real*)rho, rhoi, amp);

  for (i = 0; i < nsteps; i++) {
    step();
    fprintf(stderr, "end of step %d\n", i);
    if ((i+1)%nout == 0) output_grid(M, N, (real*)rho, rhoi, amp);
  }

  ximage_quit();

}
