
// hdyn_tt.C: hdyn-specific node-handling functions.

#include "hdyn.h"


main()
{
    hdyn *b;
    cout << "Hi \n";
    b = get_hdyn(cin);

    cout << "Hi 2\n";
    pp3(b);
    cout << "Hi 3\n";
    put_node(cout,*b);
}

