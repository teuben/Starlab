#include "stdinc.h"

#define	RAS_MAGIC	0x59a66a95

main()
{
    int i = RAS_MAGIC;
    unsigned char *c;

    c = (unsigned char*) &i;

    for (int j = 0; j < 4; j++)
	printf("%d %x\n", j, *(c+j));
}
