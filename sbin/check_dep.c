
/* Check_dep:  Check whether the first file listed on the command line
 *   	       is older (return 1 as exit status) or younger (return 0)
 *	       than any of the other files listed in the second file.
 *
 *	       NOTE: Return 0 (i.e. do nothing) if the first file does not
 *		     exist, or if any other unrecoverable error occurs..
 */

#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

void quit(int i)
{
  printf("%d\n", i);
  exit(0);
}

main(int argc, char* argv[])
{
  struct stat buf;
  unsigned long source_date;    /* Avoid use of time_t, as some UNICES
                                   don't understand it. */
  FILE *f;
  char filename[1024];

  if (argc < 3) quit(0);
  if (stat(argv[1], &buf) != 0) quit(0);
  source_date = (unsigned int) buf.st_ctime;

  f = fopen(argv[2], "r");
  if (!f) {
      fprintf(stderr, "check_dep: can't open file %s\n", argv[2]);
      quit(0);
  }

  while (fgets(filename, 1024, f)) {
    char *tail = strchr(filename, '\n');
    if(tail) *tail = '\0';
    if (stat(filename, &buf) == 0)
      if (buf.st_ctime > source_date) quit(1);
  }
  quit(0);

}
