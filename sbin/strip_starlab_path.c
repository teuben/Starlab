#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main(int argc, char *argv[])
{
    /* Remove the STARLAB_PATH stem from the input file name, if possible,
       along with the leading "/" in the remainder. */

    if (argc > 1) {
	char *file = argv[1];
	char *slab = getenv("STARLAB_PATH");
	char *s = strstr(file, slab);
	if (s)
	    printf("%s\n", s+strlen(slab)+1);
	else
	    printf("%s\n", file);
    }
}
