
// Produce data for input to make_image.

#include "stdinc.h"

main()
{
  int nx = 128, ny = 128;

    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
	  cout << i+j << " ";
}
