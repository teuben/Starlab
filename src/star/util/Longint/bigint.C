#include <stdio.h>
#include <string.h>
#include "BigInt.H"

// The BigInt(const char*) constructor.
BigInt::BigInt(const char* digitString) {
  int n = strlen(digitString);
  if (n != 0) {
    digits = new char[ndigits=n];
    char* p = digits;
    const char* q = &digitString[n];
    while (n--) *p++ = *--q - '0';
  }
  else {
    digits = new char[ndigits=1];
    digits[0] = 0;
  }
}

// The BigInt(int) constructor.
BigInt::BigInt(int n) {
  char d[3*sizeof(int)+1];
  char* dp = d;
  ndigits = 0;
  do {
    *dp++ = n%10;
    n /= 10;
    ndigits++;
  } while (n > 0);
  digits = new char[ndigits];
  register int i;
  for (i=0; i<ndigits; i++) digits[i] = d[i];
}

// The BinInt(const&) constructor.
BigInt::BigInt(const BigInt& n) {
  int i = n.ndigits;
  digits = new char[ndigits=i];
  char* p = digits;
  char* q = n.digits;
  while (i--) *p++ = *q++;
}

// The BigInt assignment operator.
void BigInt::operator=(const BigInt& n) {
  if (this==&n) return;
  unsigned i = n.ndigits;
  digits = new char[ndigits=i];
  char* p = digits;
  char* q = n.digits;
  while (i--) *p++ = *q++;
}



// The BigInt addition Operator.
BigInt BigInt::operator+(const BigInt& n) {
  int maxDigits = (ndigits>n.ndigits ? ndigits : n.ndigits) + 1;
  char* sumPtr = new char[maxDigits];
  BigInt sum(sumPtr, maxDigits);
  DigitStream a(*this);
  DigitStream b(n);
  int i = maxDigits;
  int carry = 0;
  while (i--) {
    *sumPtr = (a++) + (b++) + carry;
    if (*sumPtr>9) {
      carry = 1;
      *sumPtr -= 10;
    }
    else carry = 0;
    sumPtr++;
  }
  return sum;
}

// The BigInt::print member function.
void BigInt::print() {
  int i;
  for (i=ndigits-1; i>=0; i--) printf( "%d", digits[i]);
}

