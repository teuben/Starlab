
// mycat.C - Write stdin to stdout (to test the g++ piping bug)

#  include  <stdio.h>
#  include  <stdlib.h>
#  include  <iostream.h>

int get_line(istream & str, char * line)
{
    str.get(line, 256, '\n');
    if (str.eof()) return 0;

    char c;
    str.get(c);
    return 1;
}

main()
{
    char line[256];
    while (get_line(cin, line)) {

	cout << line << endl << flush; 		// Doesn't work!

//      fprintf(stdout, "%s\n", line);		// Works!
    }
    fprintf(stdout, "%s\n", line);
    fflush(stdout);
}
