
#include "stdinc.h"

main() {

    int n = 0;
    real x;
    real sum = 0;

    while (cin >> x) {
	n++;
	sum += x;
    }

    cout << "n, mean = " << n << " " << sum/n << endl;

}
