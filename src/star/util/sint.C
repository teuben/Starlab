
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//
// sint.C: Super Simon integer.
//         Class definition for an integer long enough to count the
//         number of orbits of a extremeley close binary (orbital period
//         of hours to minutes [maybe even milliseconds]).
//         Ment to be used in hdyn_unpert.C
//
//      S. Portegies Zwart et al.        12/1992
//      S. Portegies Zeart                2/1997
//
//======================================================================

#include "./sint.h"

sint::sint(const char* digit_stream) {
  
  int n = strlen(digit_stream);

  sign = positive;
  
  if (n != 0) {
    
    digits = new char[ndigits=n];
    char* p = digits;
    const char* q = &digit_stream[n];
    while (n--) *p++ = *--q - '0';
    
  }
  else {
    
    digits = new char[ndigits=1];
    digits[0] = 0;
    
  }
}

sint::sint(int n) {
  
  char d[3*sizeof(int)+1];
  char* dp = d;
  ndigits = 0;

  sign = positive;

  do {
    
    *dp++ = n%10;
    n /= 10;
    ndigits++;
    
  } while (n > 0);
  
  digits = new char[ndigits];
  register int i;
  for (i=0; i<ndigits; i++) digits[i] = d[i];
  
}

sint::sint(real r) {
  
  ndigits = 0;
  sign = positive;

  real r_tmp = r;

  ndigits = int(floor(log10(r)+1));
  int precision = min(ndigits, NUMERICAL_PRECISION);
  
  digits = new char[ndigits];

  for(int i=0; i<ndigits; i++) {
    if(i<ndigits-precision) {
        digits[i] = 0;
    }
    else {
      r_tmp = r/pow(10., i+1);
      digits[i] = (int)floor(10*(r_tmp-floor(r_tmp)));
      r_tmp   -= floor(r_tmp);
    }
  }
}

sint::sint(const sint& n) {
  
  int i = n.ndigits;
  digits = new char[ndigits=i];
  char* p = digits;
  char* q = n.digits;
  while (i--) *p++ = *q++;
  
}

void sint::operator=(const sint& n) {
  
  if (this==&n) return;
  
  unsigned i = n.ndigits;
  digits = new char[ndigits=i];
  char* p = digits;
  char* q = n.digits;
  
  while (i--) *p++ = *q++;
}

sint sint::operator+(const sint& n) {
  
  int max_digits = (ndigits > n.ndigits ? ndigits : n.ndigits) + 1;
  
  char* psum = new char[max_digits];
  sint sum(psum, max_digits);
  digit_stream a(*this);
  digit_stream b(n);
  int i = max_digits;
  int carry = 0;
  
  while (i--) {
    
    *psum = (a++) + (b++) + carry;
    
    if (*psum>9) {
      
      carry = 1;
      *psum -= 10;
      
    }
    else carry = 0;
    
    psum++;
    
  }

  return sum;
}

void sint::print(ostream & s) const {
  
  int i;

  if (sign==negative)
     cout << "-";
  
  for (i=ndigits-1; i>=0; i--) printf("%d", digits[i]);
    //    s << digits[i];
  //  printf(s, "%d", digits[i]);
  //s << ToString();
}

void sint::chop()
{
      int k;
      int len = NumDigits();
      for(k=len-1; k >= 0; k--)        // find a non-zero digit
	{
	  cout << "digit "<<k<<" ="<< GetDigit(k)<<endl;
	  if (GetDigit(k) != 0) break;
	  ndigits--;               // "chop" off zeros
	}
      if (k < 0)    // all zeros
	{
	  ndigits = 1;
	  sign = positive;
	}
}

int sint::NumDigits() const
// postcondition: returns # digits in sint
{
    return ndigits;
}

int sint::GetDigit(int k) const
// precondition: 0 <= k < NumDigits()
// postcondition: returns k-th digit 
//                (0 if precondition is false)
//                Note: 0th digit is least significant digit
{
    if (0 <= k && k < NumDigits())
    {
      cout << digits[k]-'0' << endl;
        return digits[k] - '0';
    }
    return 0;
}

// postcondition: returns true iff sint is negative
bool sint::is_negative() const
{
    return sign == negative;
}

bool sint::is_positive() const
// postcondition: returns true iff sint is positive
{
	return sign == positive;
}
bool sint::LessThan(const sint & rhs) const
{
    // if signs aren't equal, self < rhs only if self is negative
    
    if (is_negative() != rhs.is_negative())
    {
        return is_negative();         
    }

    // if # digits aren't the same must check # digits and sign
    
    if (NumDigits() != rhs.NumDigits())
    {
        return (NumDigits() < rhs.NumDigits() && is_positive()) ||
               (NumDigits() > rhs.NumDigits() && is_negative());
    }

    // assert: # digits same, signs the same

    int k;
    int len = NumDigits();

    for(k=len-1; k >= 0; k--)
    {
        if (GetDigit(k) < rhs.GetDigit(k)) return is_positive();
        if (GetDigit(k) > rhs.GetDigit(k)) return is_negative();
    }
    cerr<<"Is it this false?"<<endl;
    return false;      // self == rhs
}

// Overloaded operators.

// postcondition: return true if lhs < rhs, else returns false     
bool operator < (const sint & lhs, const sint & rhs)
{
    return lhs.LessThan(rhs);
}

bool operator > (const sint & lhs, const sint & rhs)
// postcondition: return true if lhs > rhs, else returns false
{
    return (rhs < lhs);    
}

ostream & operator <<(ostream & os, const sint & s)
{
  for (int i=s.NumDigits()-1; i>=0; i--) {
    os << &s.digits[i]; 
  }
  
  return os;
}


