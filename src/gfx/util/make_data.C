
// Produce data for input to make_image.

#include "stdinc.h"

#define N 256

main(int argc, char *argv[])
{
    int n = N;
    if (argc > 1) n = atoi(argv[1]);

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	  cout << i+j << " ";
}
