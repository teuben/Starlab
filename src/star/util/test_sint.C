/*
 *  sint_test.C: test the super simon integer string stream
 *.............................................................................
 *    version 1:  Feb 1997   Simon F. Portegies Zwart 
 *.............................................................................
 */
#include <iostream.h>
#include "sint.h"

main() {

  char ca[80], cb[80];
  real cc =0;

  cout << "\na= "; cin >> ca;
  cout << "b= "; cin >> cb;
  cout << "c= "; cin >> cc;

  sint a = ca;
  sint b = cb;
  sint c = cc;
  
  a.print();
  cout<<endl;
  b.print();
  cout<<endl;
  c.print();
  cout<<endl;
  cout <<" a= "<< a << " b=" << b << " c=" << c<<endl;

  a = a + b;
  b = a + c;

  a.chop();
  b.chop();
  
  a.print() ;
  cout << endl;
  b.print() ;
  cout << endl;
  
  if (a<b) {
    cout << "a= "; a.print();
    cout << " < b= ";
    b.print(); cout << endl;
  }
  else {
    cout << "a= ";
    a.print();
    cout << "> b= ";
    b.print();
    cout << endl;
  }

  if (a>b) {
    cout << "a= ";
    a.print();
    cout << "> b= ";
    b.print();
    cout << endl;
  }
  else {
    cout << "a= ";
    a.print();
    cout << "< b= ";
    b.print();
    cout << endl;
  }

}
