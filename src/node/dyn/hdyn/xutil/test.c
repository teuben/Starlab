#include <stdio.h>
#include <math.h>
#include <signal.h>

void mysig(int i)
{
    fprintf(stderr, "i = %d\n", i);
    exit(42);
}

void silly(float x, float* y)
{
    fprintf(stderr, "about to do something silly...\n");
    *y = 1/x;
}

main()
{
    float y;
    /* signal(SIGFPE, mysig); */
    silly(0.0, &y);
    fprintf(stderr, "y = %f\n", y);
}
