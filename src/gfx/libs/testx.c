
/*
 * testx:	See if we can open an X-window.
 * Author:	Steve McMillan, Drexel University, August 1995.
 */

#include "../lux/win.h"
#include <stdlib.h>

main()
{
    Display *display;
    char *display_name = NULL;
    char *disp;

    disp = getenv("DISPLAY");

    if (disp == '\0') {
	printf("DISPLAY variable not set.\n");
	exit(1);
    } else
	printf("Checking permission to open X-display %s...", disp);

    if ( (display = XOpenDisplay(display_name)) == NULL ){
	printf("denied.\n");
	exit(1);
    } else {
	printf("OK.\n");
	exit(0);
    }
}
