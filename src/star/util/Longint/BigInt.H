// BigInt header file.
// Contains defineiton os BigInt variable class.

#include <stdio.h>
#include <string.h>

// Define BigInt calls variable.
class BigInt {
  char* digits;
  int ndigits;
  BigInt(char* d, int n) {
    digits = d;
    ndigits = n;
  }
  friend DigitStream;
  public:
    BigInt(const char*);
    BigInt(int); 
    BigInt(const BigInt&); 
    void operator = (const BigInt&);
    BigInt operator + (const BigInt&);
    void print();
    ~BigInt() {delete digits;}
};

// Define DigitStream handler.
class DigitStream {
  char* dp;
  int nd;
  public: 
    DigitStream(const BigInt& n) {
      dp = n.digits;
      nd = n.ndigits;
    }
    int operator++() {
      if (nd==0) return 0;
      else {
        nd--;
        return *dp++;
      }
    }
};
