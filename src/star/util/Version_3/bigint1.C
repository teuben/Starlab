#include <iostream.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h> 
#include "bigint1.h"
#include "apvector.h"

const int BASE = 10;

// author: Owen Astrachan
//
//****************************************************************
//  bigint1.cpp
//
//  Modified January 8, 1998,  Chris Nevison
//
//  Modification removes operators -= *= (int) * (int, BigInt)
//  for use with lab for creating stubs.
//****************************************************************
//
// BigInts are implemented using a Vector<char> to store
// the digits of a BigInt
// Currently a number like 5,879 is stored as the vector {9,7,8,5}
// i.e., the least significant digit is the first digit in the vector;
// for example, GetDigit(0) returns 9 and getDigit(3) returns 5.
// All operations on digits should be done using private
// helper functions:
//
// int  GetDigit(k)        -- return k-th digit
// void ChangeDigit(k,val) -- set k-th digit to val
// void AddSigDigit(val)   -- add new most significant digit val
//
// by performing all ops in terms of these private functions we
// make implementation changes simpler
// 
// I/O operations are facilitated by the ToString() member function
// which converts a BigInt to its string (ASCII) representation


BigInt::BigInt()
// postcondition: bigint initialized to 0
   : mySign(positive),
     myDigits(1,'0'),
     myNumDigits(1)
{
	// all fields initialized in initializer list
}

BigInt::BigInt(int num)
// postcondition: bigint initialized to num
   : mySign(positive),
     myDigits(1,'0'),
     myNumDigits(0)
{ 
    // check if num is negative, change state and num accordingly
    
    if (num < 0)
    {
        mySign = negative;
        num = -1 * num;
    }

    // handle least-significant digit of num (handles num == 0)
    
    AddSigDigit(num % BASE);
    num = num / BASE;
    
    // handle remaining digits of num
    
    while (num != 0)
    {
        AddSigDigit(num % BASE);
        num = num / BASE;
    }
}


BigInt::BigInt(const apstring & s)
// precondition: s consists of digits only, optionally preceded by + or - 
// postcondition: bigint initialized to integer represented by s
//                if s is not a well-formed BigInt (e.g., contains non-digit
//                characters) then an error message is printed and abort called
  : mySign(positive),
    myDigits(s.length(),'0'),
    myNumDigits(0)
{
    int k;
    int limit = 0;

    if (s.length() == 0)
    {
        myDigits.resize(1);
        AddSigDigit(0);
        return;
    }
    if (s[0] == '-')
    {
        mySign = negative;
        limit = 1;
    }
    if (s[0] == '+')
    {
        limit = 1;
    }
    // start at least significant digit
    
    for(k=s.length() - 1; k >= limit; k--)
    {
        if (! isdigit(s[k]))
        {
            cerr << "badly formed BigInt value = " << s << endl;
            abort();
        }
        AddSigDigit(s[k]-'0');
    }
    Normalize();
}


const BigInt & BigInt::operator +=(const BigInt & rhs)
// postcondition: returns value of bigint + rhs after addition
{

    int sum;
    int carry = 0;

    int k;
    int len = NumDigits();         // length of larger addend

    if (this == &rhs)              // to add self, multiply by 2
    {
    	*this *= 2;
    	return *this;
    }

    if (rhs == 0){
      return *this;
    }

    if (IsPositive() != rhs.IsPositive())    // signs not the same, subtract
    {
        *this -= (-1 * rhs);
        return *this;
    }

    // process both numbers until one is exhausted

    if (len < rhs.NumDigits())
    {
        len = rhs.NumDigits();
    }
    for(k=0; k < len; k++)
    {
        sum = GetDigit(k) + rhs.GetDigit(k) + carry;
        carry = sum / BASE;
        sum = sum % BASE;

        if (k < myNumDigits)
        {
            ChangeDigit(k,sum);
        }
        else
        {
            AddSigDigit(sum);
        }
    }
    if (carry != 0)
    {
        AddSigDigit(carry);
    }
    return *this;    
}

BigInt operator +(const BigInt & lhs, const BigInt & rhs)
// postcondition: returns a bigint whose value is lhs + rhs
{
    BigInt result(lhs);
    result += rhs;   
    return result;
}


//****************************************************************
//
//  Make these stubs

const BigInt & BigInt::operator -=(const BigInt & rhs)
// postcondition: returns value of bigint - rhs after subtraction
{
}

const BigInt & BigInt::operator *=(int num)
// postcondition: returns num * value of BigInt after multiplication
{
}

BigInt operator *(int num, const BigInt & a)
// postcondition: returns num * a     
{
}

//****************************************************************

void BigInt::Print(ostream & os) const
// postcondition: BigInt inserted onto stream os
{
    os << ToString();
}

apstring BigInt::ToString() const
// postcondition: returns string equivalent of BigInt
{
    int k;
    int len = NumDigits();    
    apstring s = "";

    if (IsNegative())
    {
		  s = '-';
    }
    for(k=len-1; k >= 0; k--)
    {
        s += char('0'+GetDigit(k));
    }
    return s;
}

int BigInt::ToInt() const
// precondition: INT_MIN <= self <= INT_MAX
// postcondition: returns int equivalent of self
{
    int result = 0;
    int k;
    if (INT_MAX < *this) return INT_MAX;
    if (*this < INT_MIN) return INT_MIN;
    
    for(k=NumDigits()-1; k >= 0; k--)
    {
        result = result * 10 + GetDigit(k);
    }
    if (IsNegative())
    {
        result *= -1;
    }
    return result;
}

double BigInt::ToDouble() const
// precondition: DBL_MIN <= self <= DLB_MAX
// postcondition: returns double equivalent of self
{
    double result = 0.0;
    int k;
    for(k=NumDigits()-1; k >= 0; k--)
    {
        result = result * 10 + GetDigit(k);
    }
    if (IsNegative())
    {
        result *= -1;
    }
    return result;
}

ostream & operator <<(ostream & out, const BigInt & big)
// postcondition: big inserted onto stream out
{
    big.Print(out);
    return out;
}

istream & operator >> (istream & in, BigInt & big)
// postcondition: big extracted from in, must be whitespace delimited     
{
    apstring s;
    in >> s;
    big = BigInt(s);
    return in;
}

bool operator == (const BigInt & lhs, const BigInt & rhs)
// postcondition: returns true if lhs == rhs, else returns false
{
    return lhs.Equal(rhs);
}

bool BigInt::Equal(const BigInt & rhs) const
// postcondition: returns true if self == rhs, else returns false
{

    if (NumDigits() != rhs.NumDigits() || IsNegative() != rhs.IsNegative())
    {
        return false;
    }
    // assert: same sign, same number of digits;

    int k;
    int len = NumDigits();
    for(k=0; k < len; k++)
    {
        if (GetDigit(k) != rhs.GetDigit(k)) return false;
    }

    return true;
}

bool operator != (const BigInt & lhs, const BigInt & rhs)
// postcondition: returns true if lhs != rhs, else returns false     
{
    return ! (lhs == rhs);
}

bool operator < (const BigInt & lhs, const BigInt & rhs)
// postcondition: return true if lhs < rhs, else returns false     
{
    return lhs.LessThan(rhs);
}

bool BigInt::LessThan(const BigInt & rhs) const
// postcondition: return true if self < rhs, else returns false     
{
    // if signs aren't equal, self < rhs only if self is negative
    
    if (IsNegative() != rhs.IsNegative())
    {
        return IsNegative();         
    }

    // if # digits aren't the same must check # digits and sign
    
    if (NumDigits() != rhs.NumDigits())
    {
        return (NumDigits() < rhs.NumDigits() && IsPositive()) ||
               (NumDigits() > rhs.NumDigits() && IsNegative());
    }

    // assert: # digits same, signs the same

    int k;
    int len = NumDigits();

    for(k=len-1; k >= 0; k--)
    {
        if (GetDigit(k) < rhs.GetDigit(k)) return IsPositive();
        if (GetDigit(k) > rhs.GetDigit(k)) return IsNegative();
    }
    return false;      // self == rhs
}

bool operator > (const BigInt & lhs, const BigInt & rhs)
// postcondition: return true if lhs > rhs, else returns false
{
    return rhs < lhs;    
}

bool operator <= (const BigInt & lhs, const BigInt & rhs)
// postcondition: return true if lhs <= rhs, else returns false
{
    return lhs == rhs || lhs < rhs;
}
bool operator >= (const BigInt & lhs, const BigInt & rhs)
// postcondition: return true if lhs >= rhs, else returns false
{
    return lhs == rhs || lhs > rhs;
}

void BigInt::Normalize()
// postcondition: all leading zeros removed     
{
    int k;
    int len = NumDigits();
    for(k=len-1; k >= 0; k--)        // find a non-zero digit
    {
        if (GetDigit(k) != 0) break;
        myNumDigits--;               // "chop" off zeros
    }
    if (k < 0)    // all zeros
    {
        myNumDigits = 1;
        mySign = positive;
    }
}



int BigInt::NumDigits() const
// postcondition: returns # digits in BigInt
{
    return myNumDigits;
}

int BigInt::GetDigit(int k) const
// precondition: 0 <= k < NumDigits()
// postcondition: returns k-th digit 
//                (0 if precondition is false)
//                Note: 0th digit is least significant digit
{
    if (0 <= k && k < NumDigits())
    {
        return myDigits[k] - '0';
    }
    return 0;
}

void BigInt::ChangeDigit(int k, int value)
// precondition: 0 <= k < NumDigits()
// postcondition: k-th digit changed to value
//                Note: 0th digit is least significant digit
{
    if (0 <= k && k < NumDigits())
    {
        myDigits[k] = char('0' + value);
    }
    else
    {
        cerr << "error changeDigit " << k << " " << myNumDigits << endl;
    }
}

void BigInt::AddSigDigit(int value)
// postcondition: value added to BigInt as most significant digit
//                Note: 0th digit is least significant digit     
{
    if (myNumDigits >= myDigits.length())
    {
        myDigits.resize(myDigits.length() * 2);
    }
    myDigits[myNumDigits] = char('0' + value);
    myNumDigits++;
}

bool BigInt::IsNegative() const
// postcondition: returns true iff BigInt is negative
{
    return mySign == negative;
}

bool BigInt::IsPositive() const
// postcondition: returns true iff BigInt is positive
{
	return mySign == positive;
}
