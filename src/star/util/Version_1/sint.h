// sint header file.
// Contains defineiton os sint variable class.

#include <stdio.h>
#include <string.h>

class digit_stream;


class sint {
  
  protected:
  
    char* digits;
    int  ndigits;
  
    sint(char* d, int n) {
      digits = d;
      ndigits = n;
    }
    friend digit_stream;
    
  public:
    
    sint(const char*);
    sint(int); 
    sint(const sint&); 
    void operator = (const sint&);
    sint operator + (const sint&);
    void print();
    ~sint() {delete digits;}
};

class digit_stream {
  
  protected:
  
  char* dp;
  int nd;
  
  public:
  
    digit_stream(const sint& n) {
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

