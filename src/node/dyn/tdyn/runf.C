#include <stdiostream.h>

main()
{
    double t, m, x[3], v[3];

    cin.read((unsigned char *)&t, 8);
    cin.read((unsigned char *)&m, 8);
    cin.read((unsigned char *)x, 24);
    cin.read((unsigned char *)v, 24);

    cout.precision(18);
    cout << t << endl;
    cout << m << endl;
    cout << x[0] << " " << x[1] << " " << x[2] << endl;
    cout << v[0] << " " << v[1] << " " << v[2] << endl;
}
