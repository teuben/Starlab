
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Steve's replacement for the pow function (in limited circumstances).

//#include "hdyn.h"
#include "stdinc.h"

#define POW -0.125		// approximations work best for POW close to 0

#ifndef TOOLBOX

#define X0	1.0
#define X1	1.5
#define X2	2.0

static real x[3] = {X0, X1, X2};
static real p[3] = {pow(X0, POW), pow(X1, POW), pow(X2, POW)};
static real c[3];
static bool init_c = false;

local inline real g(real xx)	// 3-point Lagrange approximation
				// to xx^POW on [1,2]
{
    if (!init_c) {
	c[0] = p[0] / ((x[0]-x[1]) * (x[0]-x[2]));
	c[1] = p[1] / ((x[1]-x[0]) * (x[1]-x[2]));
	c[2] = p[2] / ((x[2]-x[0]) * (x[2]-x[1]));
	init_c = true;
    }

    return   c[0] * (xx-x[1]) * (xx-x[2])
	   + c[1] * (xx-x[0]) * (xx-x[2])
	   + c[2] * (xx-x[0]) * (xx-x[1]);
}

#define IMAX	24
static real pow_int[IMAX];
static bool init_p = false;

real pow_approx(real x)
{
    if (x < 1) return 1/pow_approx(1/x);	// !!!

    if (!init_p) {
	int j = 1;
	for (int i = 0; i < IMAX; i++) {
	    pow_int[i] = pow(j, POW);
	    j *= 2;
	}
	init_p = true;
    }

    int i = 0;
    while (x > 2) {		// 2*log2(x) operations
	x *= 0.5;
	i++;
	if (i >= IMAX-1) return pow_int[IMAX-1];
    }
    return pow_int[i] * g(x);	// 12 operations
}

#else

local real pow_true(real x)
{
    return pow(x, POW);
}

main()
{
    real x, dx = 1.e4;
    for (x = 1; x < 1.e8+.5*dx; x += dx)
	printf("%f %f %f\n",
	       log10(x), log10(pow_true(x)), log10(pow_approx(x)));
}

#endif
