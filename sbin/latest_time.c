
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

/*  Print out the (time_t) mtime of the most recently modified file
 *  specified on the command line.
 */

main(int argc, char *argv[])
{
    int i;
    time_t latest_time = -1;

    for (i = 1; i < argc; i++) {
        struct stat buffer;
	if (stat(argv[i], &buffer) == 0) {
	  if (buffer.st_mtime > latest_time)
	    latest_time = buffer.st_mtime;
	}
    }
    printf("%d\n", latest_time);
}
