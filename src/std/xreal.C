
#include "stdinc.h"

#ifndef TOOLBOX
#ifdef USE_XREAL

inline int xadd(unsigned long long &a, unsigned long long b)
{
    a += b;
    return (a < b ? 1 : 0);	// a < b ==> overflow
}

inline int xsub(unsigned long long &a, unsigned long long b)
{
    unsigned long long aa = a;
    a -= b;
    return (a > aa ? 1 : 0);	// a > aa ==> underflow
}

// Constructors:

xreal::xreal() {i = f = 0;}					// default

//xreal::xreal(xreal x) {					// x = xreal
//    i = x.i;
//    f = x.f;
//}

xreal::xreal(long long ii, unsigned long long ff) {
    i = ii;
    f = ff;
}

xreal::xreal(int ii) {						// x = int
    i = ii;
    f = 0;
}

xreal::xreal(real x)						// x = real
{
    real ii = floor(x);			// split integer/fraction
    real ff = x - ii;

    // Check for over/underflow.

    if (ii <= -TWO63) {			// comparison with TWO63N
					// here fails!

	// i = TWO63N;			// does not work on Linux:
					// RHS is interpreted as int

	// * Ugly! *

	i = 1;
	i = i<<63;			// should work (overflow/wrap)...

	f = 0;				// smallest possible xreal

    } else if (ii >= TWO63) {

	// i = TWO63M;			// does not work on linux:
					// RHS is interpreted as int

	// * Ugly! *

	i = 1;
	i = i <<63;
	i--;				// should work (overflow/wrap)...

	f = 0; f--;			// largest possible xreal

    } else {
	i = (long long)ii;
	f = (unsigned long long)(TWO64*ff);
    }
}

real xreal::to_real()
{
    // Explicitly avoid obvious roundoff problems near zero.

    if (i == -1) {

	return -((0-f)*TWO64I);	// TWO64 <--> 0!
		
    } else
	return i + f*TWO64I;	// should be OK for other real values
}

// Unary -, binary +, -, +=, -= (don't check for over/underflow):

xreal xreal::operator - () {return xreal(-i-1, 0-f);}

xreal xreal::operator + (const xreal y) {
    unsigned long long sumf = f;
    long long sumi = i + y.i + xadd(sumf, y.f);
    return xreal(sumi, sumf);
}

xreal xreal::operator - (const xreal y) {
    unsigned long long sumf = f;
    long long sumi = i - y.i - xsub(sumf, y.f);
    return xreal(sumi, sumf);
}

xreal& xreal::operator += (const xreal y) {
    i += y.i + xadd(f, y.f);
    return *this;
}

xreal& xreal::operator -= (const xreal y) {
    i -= y.i + xsub(f, y.f);
    return *this;
}
xreal xreal::operator + (const real y) {
    return *this + (xreal)y;
}

xreal xreal::operator - (const real y) {
    return *this - (xreal)y;
}

xreal& xreal::operator += (const real y) {
    *this += (xreal)y;
    return *this;
}

xreal& xreal::operator -= (const real y) {
    *this -= (xreal)y;
    return *this;
}

// Logical xreal operators ==, !=, <, <=, >, >=:

bool xreal::operator == (const xreal y) {
    return (i == y.i && f == y.f);
}

bool xreal::operator != (const xreal y) {
    return (i != y.i || f != y.f);
}

bool xreal::operator < (const xreal y) {
    return (i < y.i || (i == y.i && f < y.f));
}

bool xreal::operator <= (const xreal y) {
    return (i < y.i || (i == y.i && f <= y.f));
}

bool xreal::operator > (const xreal y) {
    return (i > y.i || (i == y.i && f > y.f));
}

bool xreal::operator >= (const xreal y) {
    return (i > y.i || (i == y.i && f >= y.f));
}

real fmod2(xreal x, real y)	// limited function:  assumes y is a power
				//		      of 2 less than 1
{
    xfrac fx = x.get_f();
    xreal xy = (xreal)y;
    xfrac fy = xy.get_f();
    xfrac r = fx%fy;
    return r*TWO64I;
}

#else

real fmod2(xreal x, real y)
{
    return fmod(x, y);
}

#endif

void xprint(xreal x,
	    ostream & s,	// default = cerr
	    bool newline)	// default = true
{
#ifdef USE_XREAL
    s << x.get_i() << "+" << x.get_f();
#else
    s << x;
#endif
    if (newline) s << endl;
}

#else

main()
{
    cerr.precision(HIGH_PRECISION);

    xreal x = 2000+M_PI, y, z;
    PRC(x); x.print(); cerr << endl;

    z = (y = x + 1);
    PRC(y); y.print(); cerr << endl;
    PRC(z); z.print(); cerr << endl;

    real dx = 1.e-16;
    y = x + dx;
    PRC(y); y.print(); cerr << endl;

    real dy = y - x;
    PRL(dy);
}

#endif
