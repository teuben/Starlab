
#include "stdinc.h"

bool twiddles(real a, real b,
	      real eps)		// default = 1.e-12
{
    if (a == b || 2*abs(a-b) <= eps*(abs(a)+abs(b)))
	return true;
    else
	return false;
}

real adjust_number_to_power(real newstep, real max_step_size)
{
    real tmp = 1.0;
    real limit = max_step_size;
    if(newstep < limit) limit = newstep;
    while(tmp < limit) tmp *= 4.0;
    while(tmp > limit) tmp *= 0.5;
    return tmp;
}

int sign(real x)
{
    if (x > 0)
      return 1;
    else if (x < 0)
      return -1;
    else
      return 0;
}

#if !defined(linux)	// already defined un Linux/math.h

real acosh(real x)  // inverse of cosh()
{
    return log( x + sqrt( x * x - 1 ) );
}

real asinh(real x)  // inverse of sinh()
{
    return sign(x) * log( abs(x) + sqrt( x * x + 1 ) );
}

#endif
