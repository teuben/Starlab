// sint header file.
// Contains definition of sint variable class.

#include <stdinc.h>
//#include <stdio.h>
#include <string.h>

#define NUMERICAL_PRECISION 16

class digit_stream;

class sint {
  
  protected:

  enum {negative, positive} sign;
  
    char* digits;
    int  ndigits;
  
    sint(char* d, int n) {
      sign = positive; 
      digits = d;
      ndigits = n;
    }
    friend digit_stream;

    bool is_negative()  const;    // return true iff number is negative
    bool is_positive()  const;    // return true iff number is positive

    bool equal(const sint & rhs)        const;
    int NumDigits()     const;    // return # digits in number
    int GetDigit(int k) const;
    int set_digit(int k) const;

  public:
    void chop();

    sint(const char*);
    sint(int); 
    sint(real); 
    sint(const sint&); 
    ~sint() {delete digits;}
    
    void operator = (const sint&);
    sint operator + (const sint&);

    bool LessThan(const sint & rhs)     const;
    
    void print(ostream &s = cerr)       const;
    friend ostream & operator <<(ostream &, const sint &);
};

bool operator <  (const sint & lhs, const sint & rhs);
bool operator >  (const sint & lhs, const sint & rhs);

class digit_stream {
  
  protected:
  
  char* dp;
  int nd;
  
  public:
  
    digit_stream(const sint& n) {
      dp = n.digits;
      nd = n.ndigits;
    }
    int operator ++ () {
      if (nd==0) return 0;
      else {
        nd--;
        return *dp++;
      }
    }
};

