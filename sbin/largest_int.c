#include <stdio.h>

main()
{
    unsigned char c = 1;
    unsigned short s = 1;
    unsigned int  i = 1;
    unsigned long l = 1;
    unsigned long long ll = 1;

#if 0

    /* Display values: */

    while (c) printf("%u\n", (c <<= 1) - 1);
    while (i) printf("%u\n", (i <<= 1) - 1);
    while (l) printf("%lu\n", (l <<= 1) - 1);
    while (ll) printf("%qu\n", (ll <<= 1) - 1);

#endif

    /* Just the largest: */

    c = -1;
    printf("largest unsigned char      = %u\n", c);

    s = -1;
    printf("largest unsigned short     = %u\n", s);

    i = -1;
    printf("largest unsigned int       = %u\n", i);

    l = -1;
    printf("largest unsigned long      = %lu\n", l);

    ll = -1;
#if 1
    printf("largest unsigned long long = %lu\n", ll);	/* Alpha/Dec UNIX */
#else
    printf("largest unsigned long long = %qu\n", ll);	/* PC/Linux */
#endif
}
