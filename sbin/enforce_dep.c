
/* Enforce_dep:  Enforce dependencies not easily handled by make,
 *   		 by "touching" any source file older than any one
 *		 of the specified dependencies.
 *
 *		 Command-line arguments are divided into "sources"
 *		 and "dependencies" by use of a separator character.
 *
 * Conventions:	 "Source" files without extensions force trials with
 *		 ".C", ".c", ".f" and ".F" added.
 *
 *		 The separator may have a directory appended to it.
 *		 This directory is prepended to all dependencies.
 */

#define DEBUG		0
#define SEPARATOR	'/'
 
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <utime.h>
#include <string.h>

#define NEXT		5

main(int argc, char* argv[])
{
  char **src = argv + 1, **dep = (char**)NULL;
  char dep_dir[128] = ".";
  char* ext[NEXT] = {".C", ".c", ".f", ".F", ""};

  int i, j, nsrc = 0, ndep = 0;
  time_t dep_date = 0;

  struct stat buf;

  /* Note: We need at least one source, at least one dependency,
           as well as the separator. */

  if (argc > 3) {

    /* Locate the separator, and determine the locations and
       lengths of the two lists. */

    for (i = 1; i < argc; i++) {
      if (argv[i][0] == SEPARATOR) {
	if (argv[i][1] != '\0') {
	  j = 1;
	  do {
	    dep_dir[j-1] = argv[i][j];
	  } while (argv[i][j++] != '\0');
	}
	nsrc = i - 1;
	dep = argv + i + 1;
	break;
      }
    }

    ndep = argc - nsrc - 2;
    if (nsrc == 0 || ndep == 0 || dep == NULL) exit(0);

    if (DEBUG) {
      printf("sources:\n");
      for (i = 0; i < nsrc; i++) printf("    %s\n", *(src+i));
    }
      printf("dependencies:\n");
      for (i = 0; i < ndep; i++) printf("    %s/%s\n", dep_dir, *(dep+i));

    /* Get the date of the newest dependency. */

    for (i = 0; i < ndep; i++) {
      char depend[128];
      sprintf(depend, "%s/%s", dep_dir, *(dep+i));
      if (stat(depend, &buf) == 0)
	if (dep_date < buf.st_ctime) dep_date = buf.st_ctime;
    }

    /* Test the source files against this date. */

    for (i = 0; i < nsrc; i++) {

      int next = NEXT - 1;

      if (index(*(src+i), '.') == NULL) {

	/* No extension.  Go through the defaults. */

	next = 0;

      }

      for (j = next; j < NEXT; j++) {

	char source[128];
	sprintf(source, "%s%s", *(src+i), ext[j]);
	if (DEBUG) printf("source = %s\n", source);

	if (stat(source, &buf) == 0) {

	  if (buf.st_ctime < dep_date) {

	    /* Out of date.  Update the file's timestamp. */

	    /*	    if (DEBUG)*/
	    printf("%s is out of date\n", source);
	    /* 	    else */
	    utime(source, (struct utimbuf *)NULL);

	  } else {

	    if (DEBUG) printf("%s is not out of date\n", source);

	  }

	}

      }

    }
  }
}
