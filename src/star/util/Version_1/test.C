#include <iostream.h>
#include "BigInt.H"

main() {
 
  void test();

  test();
}

void test() {
  BigInt b = 1;
  for (int i=1; i<1000; ++i)
    b = b + i;
}

