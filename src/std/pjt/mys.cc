#include <iostream>
#include <string>
#include <vector>

using namespace std;

static string aa[] = { "a", "b", "c" , ""};

vector<string> aaa(10);

int main (int argc, char *argv[])
{
  string a,b,c(argv[0]),d;
  int i, n;

  a = argv[0];

  cout << "String testing" << endl;
  cout << "a=" << a << endl;
  cout << "a[1]= " << a[1] << endl;
  cout << "b=" << b << endl;
  cout << "c=" << c << endl;

  i = c.find_last_of('q');
  cout << "i=" << i << endl;
  d = c.substr(i+1);
  cout << "d=" << d << endl;

  for (n=0; argv[n]; n++) {
    cout << "argv[" << n << "] = " << argv[n] << endl;
  }

  cout << sizeof(aa)/sizeof(string) << endl;
  for (i=0; aa[i].size() ; i++)
    cout << "aa: " << i << " " << aa[i] << endl;

}


