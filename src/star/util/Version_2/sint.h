#ifndef BIGINT_H
#define BIGINT_H

#include <iostream.h>
#include <string.h>

// author: Owen Astrachan
// modified by Claire Bono
//
// implements an arbitrary precision integer class
//
// supports non-negative numbers only
//
// constructors:
//
// sint()            -- default constructor, value of integer is 0
// sint(int n)       -- initialize to value of n (C++ int)
//       precondition: n >= 0
//
// sint(const string & s) -- initialize to value specified by s
//        e.g.    sint big("121212121212123123123123123");
//                (1 is most sig. digit of big; 3 is least sig. digit of big)
//        precondition: s has only digit chars and is not empty
//
// sint(const sint & b)  -- copy constructor
//
//
// *****  arithmetic operators:
//   
// void operator += (const sint & b)
//          modifies object by adding a sint to it
//
// sint operator + (const sint & a, const sint & b)
//          adds two sints together, returns result
//
// void operator *= (const sint & b)
//          modifies object by multiplying by sint
//
// sint operator * (const sint & a, const sint & b)
//          multiplies two sints together, returns result
//
//  ***** logical operators:
//
// int operator == (const sint & a, const sint & b)
//          returns 1 if a == b, 0 otherwise
//
// int operator != (const sint & a, const BitInt & b)
//          returns 1 if a != b, 0 otherwise
//
// int operator < (const sint & a, const sint & b)
//          returns 1 if a < b, 0 otherwise
//
// int operator <= (const sint & a, const sint & b)
//          returns 1 if a <= b, 0 otherwise
//
// int operator > (const sint & a, const sint & b)
//          returns 1 if a > b, 0 otherwise
//
// int operator >= (const sint & a, const sint & b)
//          returns 1 if a >= b, 0 otherwise
//
//
//  ***** I/O operators:
//
//  ostream & operator << (ostream & os, const sint & b)
//        stream operator to print value
//

class sint
{
  public:
    
    sint();                  // default constructor, value = 0
    sint(int);               // assign an integer value
    //sint(const string &);    // assign a string 
    sint(const sint &);    // copy constructor
    ~sint();                 // destructor

    // operators: arithmetic, relational

    void operator += (const sint &);
    void operator *= (const sint &);
    sint & operator = (const sint &);

    friend ostream & operator <<(ostream &, const sint &);
    friend int operator == (const sint &, const sint &);
    friend int operator < (const sint &, const sint &);
  private:
    static const int INTDIGITS;
    void grow(int size);
    int capacity;
    int numDigits;
    int *digits;
};


sint operator +(const sint &, const sint &);
sint operator *(const sint &, const sint &);
int operator != (const sint & a, const sint & b);
int operator > (const sint &, const sint &);
int operator >= (const sint &, const sint &);
int operator <= (const sint &, const sint &);


#endif   // BIGINT_H not defined
