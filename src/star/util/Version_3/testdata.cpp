//  testdata.cpp
//
//  January 9, 1998
//  Chris Nevison
//
//  This program is designed to test the arithmetic operators + and -
//  and the comparison operator < for the BigInt class.

#include <iostream.h>
#include <fstream.h>
#include <apstring.h>
#include "bigint.h"

int main()
{
  BigInt x, y;
  apstring fname;
  ifstream in;
  ofstream out;

  cout << "Enter input file name: ";
  cin  >> fname;
  in.open(fname.c_str());

  cout << "Enter output file name: ";
  cin  >> fname;
  out.open(fname.c_str());

  while(in >> x >> y){
    out << x  << "   " << y << endl;
    out << "     " << x + y;
    out << "     " << x - y << "    ";
    if(x < y)
      out << " true" << endl;
    else
      out << "false" << endl;
  }

  return 0;
}

