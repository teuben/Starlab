#include <stdiostream.h>

main()
{
    double x;
    cin.read((unsigned char *)&x+4, 4);		// second 4 bytes
    cout << x << endl;
}
