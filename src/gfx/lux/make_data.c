
/* Produce a lot of data for input to ximage. */

#include <stdio.h>

main(int argc, char* argv[])
{
  int nx = 128, ny = 128, nk = 128;
  int i, j, k;

  if (argc > 1) nk = atoi(argv[1]);

  for (k = 0; k < nk; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
	putchar((char)(i+j+k));
}
