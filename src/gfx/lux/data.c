
#include <stdio.h>
#include <math.h>

main()
{
    int i;
    for (i = 0; i < 10000; i++)
	printf("%d %f %f\n", i, sin((M_PI*i)/1000.), cos((M_PI*i)/1000.));
}
