
#include "stdinc.h"

#ifndef TOOLBOX
#ifdef USE_XREAL

inline int xadd(xfrac_t &a, xfrac_t b)
{
    a += b;
    return (a < b ? 1 : 0);	// a < b ==> overflow
}

inline int xsub(xfrac_t &a, xfrac_t b)
{
    xfrac_t aa = a;
    a -= b;
    return (a > aa ? 1 : 0);	// a > aa ==> underflow
}

// Constructors:

xreal::xreal() {i = f = 0;}					// default

//xreal::xreal(xreal x) {					// x = xreal
//    i = x.i;
//    f = x.f;
//}

xreal::xreal(xint_t ii, xfrac_t ff) {
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
	i = (xint_t)ii;
	f = (xfrac_t)(TWO64*ff);
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
    xfrac_t sumf = f;
    xint_t sumi = i + y.i + xadd(sumf, y.f);
    return xreal(sumi, sumf);
}

xreal xreal::operator - (const xreal y) {
    xfrac_t sumf = f;
    xint_t sumi = i - y.i - xsub(sumf, y.f);
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
    xfrac_t fx = x.get_f();
    xreal xy = (xreal)y;
    xfrac_t fy = xy.get_f();
    xfrac_t r = fx%fy;
    return r*TWO64I;
}

#else

real fmod2(xreal x, real y)
{
    return fmod(x, y);
}

#endif


//----------------------------------------------------------------------
//
// All the xreal I/O functions should appear here, for consistency.

// Input functions are derived from a single get_xreal function.

#ifdef USE_XREAL

// Read an xreal from a string.
// C++ input apparently doesn't understand 0xa123bcde format for hex...
// Operator not widely used; use scanf until this is resolved.

static bool print_xreal = true;

local inline xreal read_xreal(const char *str)
{
    // Extract an xreal from a string, with tests for various formats.
    // Should be able to understand
    //
    //		int int
    //		int hex
    //		real

    char *sp, *ep;
    // PRL(str);
    long long i = STRTOL(str, &sp, 10);		  // signed integer part
    // PRC(i); PRL(sp);
    unsigned long long f = STRTOUL(sp, &ep, 0);   // unsigned fractional part
						  // "0" here means that we
						  // can read hex or integer
    // PRC(f); PRL(ep);

    if (sp == ep) {				  // if we didn't get both
						  // of above,

	// Hmmm... most likely we have real input data.  Try just reading
	// a real number.  (Steve, 6/00)

	if (print_xreal) {
	    cerr << "read_xreal: error reading xreal input "
		 << "from string" << endl
		 << "    " << str << endl
		 << "Assuming real data." << endl << endl;
	    print_xreal = false;
	}

	return (xreal)strtod(str, NULL);
    }

    return xreal(i, f);
}

xreal get_xreal(char *str)
{
    return read_xreal(str);
}

// Note that the >> operator is not symmetric with <<...

istream & operator >> (istream & s, xreal & x)
{
// C++ input apparently doesn't understand "0xa123bcde" format for hex...

// xint_t i;
// xfrac_t f;
// s >> i >> f;			// s >> x.i >> x.f fails; don't know why...
// x = xreal(i, f);

    // Operator is not widely used; just use read_xreal().

    // Start by reading in two "words".  Easiest to use STL strings for this.

    string i, f;
    s >> i >> f;
    string str = i + " " + f;
    x = read_xreal(str.c_str());

    return s;
}

#endif

xreal get_xreal_from_input_line(char * input_line)
{
    char *val = strchr(input_line, '=');
    if(val == NULL) return (xreal)0;
    val++;

#if defined USE_XREAL

    // "True" xreal:

    return read_xreal(val);

#else

    // xreal is really just real:

    return (xreal)strtod(val, NULL);

#endif
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// All xreal "print" functions are really versions of a single xprint.

void xprint(xreal x,
	    ostream & s,	// default = cerr
	    bool newline)	// default = true
{
#ifdef USE_XREAL

    // Simplest version:
    //
    // s << label << x.get_i() << " " << x.get_f() << endl;
    //
    // Better: Use hex (leading 0x) for the fractional part...

    xfrac_t f = x.get_f();
    char tmp[128];
    if (f == 0)
	sprintf(tmp, "0");		// handy
    else
	sprintf(tmp, "%#16.16llx", f);

    s << x.get_i() << " " << tmp;

#else

    s << x;

#endif

    if (newline) s << endl;
}

#ifdef USE_XREAL

// Member function xprint is rarely used, but convenient to retain it.
// Simply define it in terms of the non-member xprint function.

void xreal::print(ostream& s)
{
    xprint(*this, s, false);
}

// Xreal version of put_real_number is just xprint with a label.

void put_real_number(ostream & s, char * label, xreal x)
{
    s << label;
    xprint(x, s);
}
#endif

//----------------------------------------------------------------------


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
