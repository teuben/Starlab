
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
//         of hours to minutes [maybe even miliseconds]).
//         Ment to be used in hdyn_unpert.C
//
//      S. Portegies Zwart et al.        12/1992
//      S. Portegies Zeart                2/1997
//
//======================================================================

//#include <stdio.h>
//#include <string.h>
#include "./sint.h"

sint::sint(const char* digit_stream) {
  
  int n = strlen(digit_stream);
  
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
  
  do {
    
    *dp++ = n%10;
    n /= 10;
    ndigits++;
    
  } while (n > 0);
  
  digits = new char[ndigits];
  register int i;
  for (i=0; i<ndigits; i++) digits[i] = d[i];
  
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
  
  int max_digits = (ndigits>n.ndigits ? ndigits : n.ndigits) + 1;
  
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

void sint::print() {
  
  int i;
  
  for (i=ndigits-1; i>=0; i--) printf( "%d", digits[i]);
  
}

