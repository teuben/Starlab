// sint.C

// no negatives allowed

// representation:
//  int capacity  -- how big is dynamic array "digits"
//  int numDigits -- how much of the array has valid data (see invariants)
//  int *digits    -- dynamic array to store digits (see invariants)

// representation invariants: 
//        numDigits <= capacity
//        Data only valid in digits[0..numDigits-1]
//        Any valid digit[i] contains only a single digit value
//              i.e. a number in the range 0..9
//        Least significant digit (lsd) is in position 0.
//        Most sig. digit (msd) is in position numDigits-1.
//                                             0   1   2   3  <-- array index
//        Eg: for the number 8697: digits is | 7 | 9 | 6 | 8 |
//                                 and numDigits is 4
//                                 and capacity >= 4

#include <assert.h>

#include "sint.h"


// note: this is how we initialize a static const member

const int sint::INTDIGITS = 10;  // number of decimal digits in largest
                                   // 32 bit int


sint::sint() : capacity(INTDIGITS),
                   numDigits(1)
{
  digits = new int[capacity];
  digits[0] = 0;  
}


sint::sint(const sint & src): capacity(src.capacity),
                                    numDigits(src.numDigits)
{

  digits = new int[capacity];

  for (int i = 0; i < numDigits; i++) {
    digits[i] = src.digits[i];
  }
}


sint::~sint()
{
  delete [] digits;
}


sint operator +(const sint &a, const sint &b)
{
  sint sum = a;      // implemented in terms of +=
  sum += b;            // can't call until += is implemented
  return sum;
}


// private
// pre: size >= 0 and this is a valid sint (see invariant)
// post: (1) if (capacity < size) before, then (capacity'== size)
//           if (capacity >= size) before, then sint unchanged
//       (2) this is a valid sint with the same "value" as before the call
//       (i.e., all values in 0..numDigits-1 part of the digits array 
//        are unchanged and numDigits is unchanged)
void sint::grow(int size)
{
  assert (size >=0);

  if (capacity < size) {
    int *tmp = new int[size];
    for (int i = 0; i < numDigits; i++) {
      tmp[i] = digits[i];
    }
    delete [] digits;
    digits = tmp;
    capacity = size;

  }
}


ostream & operator <<(ostream &os, const sint &b)
{
  // need to print msd to lsd (i.e., print digits[0] last)

  for (int i = b.numDigits-1; i >= 0; i--) {
    os << b.digits[i];
  } 

  return os;
}


// NOTES on *=: 
//      1. you must *= reimplement for final version
//      2. to call this version *=, need to have implemented other member
//         functions: +=, !=, int->big constr., copy constr, =

void sint::operator *=(const sint & big)
// precondition: sint = a
// postcondition: sint = a*big     
{
  // alg: use repeated addition

    sint k(0);             // used as counter
    sint limit(big);       // in case of aliasing (a*=a;), use copy
    sint copy(*this);      // will keep adding in orig value

    *this = 0;               // set self to 0 for accumulation
    while (k != limit)
    {
        *this += copy;
        k += 1;
    }
}
