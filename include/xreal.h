#ifndef  STARLAB_XREAL_H
#  define  STARLAB_XREAL_H

// xreal:  Extended-precision "real" class (SLWM, 5/00).
//
//	   Currently implemented in fixed precision as
//
//		x  =  i  +  2^{-64} f,
//
//	   where the "integer" i is a long long (range 2^63 - 1 to -2^63), and
//	   the "fraction" f is an unsigned long long (range 0 to 2^64 - 1),
//	   giving roughly 20 digits before and after the decimal point.
//
//	   Permitted/defined operations:
//
//		initialized to 0 by default
//		xreal = xreal
//		xreal = int
//		xreal = real
//		int = xreal; Int(xreal), (int) xreal
//		real = xreal; Real(xreal), (real) xreal
//		-xreal
//		xreal +/-/+=/-= xreal
//		xreal +/-/+=/-= int
//		xreal +/-/+=/-= real
//		logical ==, !=, <, <=, >, >=
//		<<, >>
//
// **** Note that long long may not be universally acceptable...  ****
// ****  -- could rewrite to use standard long ints if necessary. ****

//#define TWO64	(pow(2.0, 64))	// too inefficient!

#define TWO64	(18446744073709551616.0)
#define TWO63	( 9223372036854775808.0)
#define TWO64I	(5.42101086242752217e-20)
#define TWO63M	( 9223372036854775807)
#define TWO63N	(-9223372036854775808)

typedef long long		xint_t;
typedef unsigned long long	xfrac_t;

class xreal {

    private:
	xint_t i;					// x = i + f/TWO64
	xfrac_t f;

    public:

	// Constructors:

	xreal();					// default
//	xreal(xreal y);					// x = y
	xreal(xint_t ii, xfrac_t ff);
	xreal(int ii);					// x = int
	xreal(real x);					// x = real

	// Accessors:

	xint_t get_i() {return i;}
	xfrac_t get_f() {return f;}
	real get_frac() {return TWO64I*f;}

	// Conversions:

	int to_int() {return i;}			// always round down!
	real to_real();

	void print(ostream& s = cerr,
		   bool newline = false);		// xprint(x) is better

	xreal::operator int() {return to_int();}
	xreal::operator real() {return to_real();}

	// Unary -, binary +, -, +=, -= (don't check for over/underflow):

	xreal  operator -  ();
	xreal  operator +  (const xreal y);
	xreal  operator -  (const xreal y);
	xreal& operator += (const xreal y);
	xreal& operator -= (const xreal y);
	xreal  operator +  (const real y);
	xreal  operator -  (const real y);
	xreal& operator += (const real y);
	xreal& operator -= (const real y);

	// Logical operators ==, !=, <, <=, >, >=:
	// (see also mixed operators below)

	bool operator == (const xreal y);
	bool operator != (const xreal y);
	bool operator <  (const xreal y);
	bool operator <= (const xreal y);
	bool operator >  (const xreal y);
	bool operator >= (const xreal y);

	friend istream & operator >> (ostream & , xreal &);
};

inline ostream & operator << (ostream & s, xreal x)
{
    return s << (real)x;	// probably what we want most of the time
				// -- use x.xprint() for complete ouptut
}

istream & operator >> (istream & s, xreal & x);

// Various real/xreal comparisons (rule: always promote to xreal).

inline bool operator == (xreal x, real y) {return (x == (xreal)y);}
inline bool operator == (real x, xreal y) {return ((xreal)x == y);}
inline bool operator != (xreal x, real y) {return (x != (xreal)y);}
inline bool operator != (real x, xreal y) {return ((xreal)x != y);}
inline bool operator <  (xreal x, real y) {return (x < (xreal)y);}
inline bool operator <  (real x, xreal y) {return ((xreal)x < y);}
inline bool operator >  (xreal x, real y) {return (x > (xreal)y);}
inline bool operator >  (real x, xreal y) {return ((xreal)x > y);}
inline bool operator <= (xreal x, real y) {return (x <= (xreal)y);}
inline bool operator <= (real x, xreal y) {return ((xreal)x <= y);}
inline bool operator >= (xreal x, real y) {return (x >= (xreal)y);}
inline bool operator >= (real x, xreal y) {return ((xreal)x >= y);}

// Exclude the following definitions in favor of explicit casts if needed:
//
// inline real operator + (xreal x, real y) {return (real)(x + (xreal)y);}
// inline real operator - (xreal x, real y) {return (real)(x - (xreal)y);}

xreal get_xreal(char *str);
void put_real_number(ostream & s, char * label, xreal x);

void identify_xreal(ostream& s = cerr);

#endif
