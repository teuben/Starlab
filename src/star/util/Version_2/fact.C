// Factorial using sints
// by Owen Astrachan

#include "sint.h"
#include <iostream.h>
#include <stdlib.h>        // for atoi

// illustrates use of sint class
    
// precondition: n >= 0
// postcondition: returns n!     
sint Fact(int n) {
    sint prod;
    prod = 1;
    int k;
    for(k=1; k <= n; k++)
    {
        prod *= k;
    }
    return prod;
}

#ifdef TOOLBOX

main(int argc, char * argv[])
{
    int k;
    int limit;
    sint val;

    if (argc > 1){              // command line args?
        limit = atoi(argv[1]);
    }
    else{
        cout << "enter limit ";
        cin >> limit;
    }

    cout << "number \t factorial" << endl;
    cout << "------ \t ---------" << endl << endl;
    for(k=limit; k >= 1; k--)
    {
        val = Fact(k);
        cout << k << "\t" << val << "\t" << endl;
    }
    
    return 0;
}

#endif
