
c
c	Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c       by Steve McMillan, Drexel University, Philadelphia, PA.
c
c       All rights reserved.
c
c	Redistribution and use in source and binary forms are permitted
c	provided that the above copyright notice and this paragraph are
c	duplicated in all such forms and that any documentation,
c	advertising materials, and other materials related to such
c	distribution and use acknowledge that the software was developed
c	by the author named above.
c
c	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c

	subroutine nextcmd(string,istart,ns,prompt)
        save
c
c	Simple parser. Return in string(1:ns) the next input command.
c	Commands may be issued several to a line, delimited with a semicolon.
c	If no commands are left in the current line, prompt for and read in
c	a new command line.
c
c       This routine is used by the "2dplot" subpackage.
c
	save
	character*(*) string,prompt
	character*200 line
c
	data iscolon/0/nl/-1/line/' '/
c
	if (ns.lt.0) line = ' '
c
10	if (iscolon.gt.nl) then
c
            call devoff
	    call getstring(prompt,1,len(prompt),line)
c
c           Strip off trailing blanks and non-significant semicolons.
c
	    nl = len(line)
	    iscolon = nl+1
	    call stripbl(line,nl,*10,*10)
	    iscolon = 0
c
  	end if
c
c 	Split the line into individual commands.
c
	string = ' '
  	call getnext(line,nl,iscolon,ic0,ic1,string,ns)
c
c 	Strip leading blanks, convert the current command to lowercase,
c 	and locate the start of the argument list.
c
  	call cleanup(string,ns,istart,*10)
c
c 	Note:	The entire command is string(1:ns)
c		The argument list is string(istart:ns)
c
	end
