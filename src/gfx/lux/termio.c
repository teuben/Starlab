
/* From old X interface for McDraw by Steve McMillan*/

#include "win.h"
#include <string.h>

#include <termio.h>		/* Going to use getchar with timeout */

/*
 * (Use termio rather than termios here, because the HP can't
 *  get/set termios structures.)
 *
 * Note that some systems use sys/termio.h...
 */

#include <fcntl.h>

#define	MYCMASK 0377
#define BACKSPACE 8
#define TAB 9
#define LINE_FEED 10
#define CARRIAGE_RETURN 13
#define REFRESH 18
#define KILL_LINE 21
#define KILL_WORD 23
#define DEL 127

#define IDLESTRNGCNT_MAX 255

/* struct termio ios0, ios1; */

				/* Hmmm... */

struct termios ios0, ios1;      /* NEED termios here on DEC UNIX... */
                                /* NEED termio  here on HP UNIX... */

char idlestrng[256];
int  fd, idlestrngcnt;

get_timeout()

/*
   Read a character from a terminal stream, with timeout.
   Accumulate string in global variable idlestrng
 */

{
    char c;

    if ( (c = getchar_ub()) > 0) return(process_char(c));

    return(0);

    /* c > 0 ==> a character was read before timeout */
}

getchar_ub()	/* unbuffered input */
{
    char c;

    return ((read(0, &c, 1) > 0) ? c & MYCMASK : EOF);
}

process_char(c)  /* Add char to storage string, deal with display. */
char c;
{
    int i,idletog = 0;

    if (c == '\b' || c == DEL) {

	/* Backspace/delete ==> erase character */

	if (idlestrngcnt > 0)
	    erase_char();
	else {
/*	    putchar('\a');
	    fflush(stdout);	*/	/* Sun CC doesn't understand \a!! */
	}
    } 
    else if (c == KILL_WORD) {

	/* Delete word. */

	if (idlestrngcnt > 0)
	    erase_word();
	else {
/*	    putchar('\a');
	    fflush(stdout);	*/
	}
    }
    else if (c == REFRESH) {

	/* Erase and rewrite the line */

	putchar('^');putchar('R');putchar('\n');
	for (i = 0; i < idlestrngcnt; i++) putchar(idlestrng[i]);
	fflush(stdout);
    } 
    else if (c == KILL_LINE) {

	/* Erase the line */

	if (idlestrngcnt == 0) {
/*	    putchar('\a');
	    fflush(stdout);	*/
	}
	while(idlestrngcnt) erase_char();
	idlestrng[0] = '\0';
    } 
    else if (c == '\n' || (c >= ' '&& c <= '~')) {

	idlestrng[idlestrngcnt] = c;
	idlestrng[++idlestrngcnt] = '\0';
	
	if (c == '\n') idletog++;
	
	if (idlestrngcnt >= IDLESTRNGCNT_MAX) {
	    idlestrng[idlestrngcnt] = '\n';
	    idletog++;
	}
	putchar(c);
	fflush(stdout);
    }
    return(idletog);
}

erase_char()  /* Erase a character from the display. */
{
    idlestrngcnt--;
    putchar('\b');
    putchar(' ');
    putchar('\b');
    fflush(stdout);
    idlestrng[idlestrngcnt] = '\0';
}

erase_word()  /* Erase a work from the display. */
{
    int i,is = 0,len,lens;

    len = idlestrngcnt;

    i = idlestrngcnt-1;
    for(i=idlestrngcnt-1;i>=0;i--) if (idlestrng[i] != ' ') break;
    lens = i;

    for(i=0;i<lens;i++) 
      if (idlestrng[i] == ' ' || idlestrng[i] == ';') is = i;

    if (is!=0) is = is+1;
    for(i=is;i<len;i++) erase_char();
}


set_timeout()
{

    /* Identify the terminal stream and get its current configuration. */

    if ( (fd = open("/dev/tty", O_RDONLY)) < 0 ) exit(1);
    ioctl(fd, TCGETA, &ios0);
    ioctl(fd, TCGETA, &ios1);

    /* Turn off ICANON mode and inhibit ECHO (set bits 1 and 3 = 0). */

    ios1.c_lflag &= 0xfff5;

    /* Set MIN = 0, TIME = 1 * (1/10 secs). */

    ios1.c_cc[VMIN] = 0;
    ios1.c_cc[VTIME] = 1;

    ioctl(fd, TCSETA, &ios1);

    idlestrngcnt = 0;
    idlestrng[0] = '\0';
}

reset_term(strng)
char *strng;
{

    /* Restore original terminal settings. */

    ioctl(fd, TCSETA, &ios0);

    close(fd);

    strcat(strng,idlestrng);
}

clear_buffer()
{
    idlestrng[0] = '\0';
    idlestrngcnt = 0;
}
