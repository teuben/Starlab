
#include "stdinc.h"

local void print(bool a, bool b)
{
    cerr << "a = " << a << "  b = " << b << endl;
}

main()	// Mixing bool, TRUE, true, arithmetic and logical operations!
{
    bool a, b;

    a = true;
    b = false;
    print(a, b);

    a = FALSE;
    b = !b;
    print(a, b);

    a = 1 - a;
    b = 1 - b;
    print(a, b);

    a = !a || !b;
    b = TRUE;
    print(a, b);
}
