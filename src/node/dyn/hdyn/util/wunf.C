#include <stdiostream.h>
#include <math.h>

main()
{
    double x = M_PI;
    cout.write((unsigned char *)&x+4, 4);	// second 4 bytes
}
