#include <iostream.h>
#include "BigInt.H"

main() {

  char a[80], b[80];

  cout << "\na= "; cin >> a;
  cout << "b= "; cin >> b;

  BigInt e = a;
  BigInt f = b;
  BigInt d = e + f;

  d.print();
  printf("\n");
}

