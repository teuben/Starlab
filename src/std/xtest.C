
#include "stdinc.h"

main()
{
#if 0
    long long i = 1;
    int n = 0;
    while (1) {
	PRC(n), PRL(i);
	if (i < 0) break;
	n++;
	i *= 2;
    }

    i = 1;
    i = i<<63;
    PRL(i);
    PRL(sizeof(i));
#else
    xreal x;
    cin >> x;
    cerr.precision(20);
    PRL((real)x);
#endif
}
