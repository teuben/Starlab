
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

	subroutine getnext(line,nl,iscolon,ic0,ic1,input,nin)
c
c	Extract the next semicolon-delimited piece from the input line.
c
c       Start at iscolon, which becomes ic0; on return, next semicolon
c       will be at ic1.  The string "input" has NO trailing semicolon.
c
	character*(*) line,input
c	
	ic0  =  iscolon
c
15      do ic1 = iscolon+1,nl
            if (line(ic1:ic1).eq.';') go to 20
        end do
        ic1 = nl+1
c
20      if (ic1.lt.nl.and.line(ic1+1:ic1+1).eq.';') then
c
c	    Map ";;" to ";" and continue.
c
            do i = ic1+1,nl-1
                line(i:i) = line(i+1:i+1)
            end do
            nl = nl- 1
            iscolon = ic1
            go to 15
        end if
c
        nin = ic1-ic0-1
        input(1:nin) = line(ic0+1:ic1-1)
        iscolon = ic1
c
	end
