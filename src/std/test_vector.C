#include  "starlab_vector.h"

//// test_vector:  test Starlab vector class operators
////
//// No options

#define POST_GCC_30

main(int argc, char** argv)
{
    check_help();

    vector x(1.0);
#ifdef POST_GCC_30
    vector y(1.0+x);
#else
    vector y(x+1.0);
#endif // POST_GCC_30
    vector z(x+y);

    cout << "x = "; x.print();
    cout << "y = "; y.print();
    cout << "z = "; z.print();

    cout << "-x = "; (-x).print();

    cout << "y.z = " << y*z << "\n";

    cout << "x+y = "; (x+y).print();
    cout << "x-y = "; (x-y).print();
#ifdef POST_GCC_30
    cout << "x+1 skipped" << endl;
#else
    cout << "x+1 = "; (x+1).print();
#endif //  POST_GCC_30
    cout << "1+x = "; (1+x).print();
    cout << "x+y+z = "; (x+y+z).print();

    cout << "2*x = "; (2*x).print();
#ifdef POST_GCC_30
    cout << "x*2 skipped" << endl;
#else
    cout << "x*2 = "; (x*2).print();
#endif // POST_GCC_30

#ifdef POST_GCC_30
    cout << "x/2 skipped" << endl;
    cout << "2/x skiped" << endl;
    cout << "1/2 * x = "; (1./2 * x).print();
#else
    cout << "2/x = "; (2/x).print();
    cout << "x/2 = "; (x/2).print();
#endif // POST_GCC_30

    z = x;
    cout << "z = "; z.print();
    cout << "|z| = " << abs(z) << "\n";

    cout << "|2x-4(y+z)| = " << abs(2*x-4*(y+z)) << "\n";

    cout << "test of  x = vector(1001, 1002, 1003);  : ";
    x = vector(1001, 1002, 1003);
    x.print();
    cout << "test of  x = 8;  : ";
    x = 8;
    x.print();

    vector a(1,2,3);

    cout << "test of a = (1,2,3) :    a = " << a << endl;
    cout << "test of a[2]        : a[2] = " << a[2] << endl;
    cout << "test of a[2] = 5    :    a = " ;
    a[2] = 5;
    cout << a << endl;
}
