/*
 *  sint_test.C: test the super simon integer string stream
 *.............................................................................
 *    version 1:  Feb 1997   Simon F. Portegies Zwart 
 *.............................................................................
 */
#include <iostream.h>
#include "sint.h"

#ifdef TOOLBOX

main() {

  char a[80], b[80];

  cout << "\na= "; cin >> a;
  cout << "b= "; cin >> b;

  sint e = a;
  sint f = b;
  sint d = e + f;

  d.print();
  printf("\n");
}

#endif
