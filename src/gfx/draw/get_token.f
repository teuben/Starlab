
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

c------------------------------------------------------------------------
c
c       String-decomposition and reading.
c
c       Contents:    gettokens  - split a string into pieces
c                    readrtoken - read a real number from a string
c                    readitoken - read an integer from a string
c
c------------------------------------------------------------------------

	subroutine gettokens(string,token,ntoken)
        save
c
c	Extract a list of separate "words" from a given string.
c
	character*(*) string,token(1)
c
	parameter (NSEPMAX = 3)
	character*1 sep(NSEPMAX),c
c
        common /token delim/ nsep
c
c       Present delimiters are tab, space, and comma, in that order.
c       Thus, setting nsep = 3 allows all three; nsep = 2 allows
c       "whitespace" only; nsep = 1 allows only tabs.
c
	data sep/'	',' ',','/nsep/3/
c
	ntoken = 0
	itoken = 0
	do 100 i=1,len(string)+1
	    if (i.le.len(string)) then
		c = string(i:i)
	    else
		c = ' '
	    end if
c
c	    Is this character a separator?
c
	    do 25 k=1,nsep
		if (c.eq.sep(k)) go to 50
25	    continue
c
	    if (itoken.eq.0) istart = i
	    itoken = 1
	    go to 100
c
50	    if (itoken.eq.1) then
		ntoken = ntoken + 1
		token(ntoken) = string(istart:i-1)
	    end if
	    itoken = 0
100	continue
c
	end


	subroutine readrtoken(token,xread,xx)
        save
	character*(*) token
        character*1 per
        character*1000 temp
c
        common /read_token_stat/ istat
c
c	Read a real quantity xread from the specified token, subject
c	to certain rules:
c
c		=	--> do nothing
c		number	--> use number
c		x	--> 1
c		y	--> 2
c		z	--> 3
c		. or *	--> use xx
c		>number	--> use xx, with number as a lower limit
c		<number	--> use xx, with number as an upper limit
c
c		anything else (for now) means "do nothing"
c
c	Assume it is OK to tack a period onto the string, following
c       the last nonblank character.  Note that we should make NO
c       assumptions about what follows the digits in the string.
c
c       Find the last legal character in the string.
c
        do iend=1,len(token)
            if (token(iend:iend).le.' ') go to 100
        end do
100     iend = iend - 1
c
	if (token(1:1).eq.'*'.and.iend.eq.1) token(1:1) = '.'
c
	isingledot = 0
	if (token(1:1).eq.'.'.and.iend.eq.1) isingledot = 1
	if (token(1:1).eq.'-'.and.iend.eq.1) isingledot = 1
c
c       Add a trailing period, if necessary (assume < 1000 characters):
c
        per = ' '
	if (index(token(1:iend),'.').eq.0) per = '.'
c
        istat = 0
c
c       Trying to read ". " or "- " will get NaN and io = 0...
c
        temp = token(1:iend)//per
	if (isingledot.eq.0) then
	    read(temp(1:iend+1),'(f32.16)',iostat=io)yy
	else
	    io = 1
	end if
c
	if (io.eq.0) then
	    xread = yy
	else if (token(1:1).eq.'x'.or.token(1:1).eq.'X') then
	    xread = 1.
	else if (token(1:1).eq.'y'.or.token(1:1).eq.'Y') then
	    xread = 2.
	else if (token(1:1).eq.'z'.or.token(1:1).eq.'Z') then
	    xread = 3.
	else if (token(1:1).eq.'=') then
	else if (token(1:1).eq.'.') then
	    xread = xx
	else if (token(1:1).eq.'>'.and.iend.gt.2) then
	    read(temp(2:iend+1),'(f32.16)',iostat=io)yy
	    if (io.eq.0) xread = max(xx,yy)
            istat = io
	else if (token(1:1).eq.'<'.and.iend.gt.2) then
	    read(token(2:iend+1),'(f32.16)',iostat=io)yy
	    if (io.eq.0) xread = min(xx,yy)
            istat = io
        else
            istat = io
	end if
c
	end


	subroutine readdtoken(token,xread,xx)

c       IDENTICAL to readrtoken, but for real*8 xread.

        save
	character*(*) token
        character*1 per
        character*1000 temp

	real*8 xread,yy
c
        common /read_token_stat/ istat
c
c	Read a real quantity xread from the specified token, subject
c	to certain rules:
c
c		=	--> do nothing
c		number	--> use number
c		x	--> 1
c		y	--> 2
c		z	--> 3
c		. or *	--> use xx
c		>number	--> use xx, with number as a lower limit
c		<number	--> use xx, with number as an upper limit
c
c		anything else (for now) means "do nothing"
c
c	Assume it is OK to tack a period onto the string, following
c       the last nonblank character.  Note that we should make NO
c       assumptions about what follows the digits in the string.
c
c       Find the last legal character in the string.
c
        do iend=1,len(token)
            if (token(iend:iend).le.' ') go to 100
        end do
100     iend = iend - 1

	if (token(1:1).eq.'*'.and.iend.eq.1) token(1:1) = '.'

	isingledot = 0
	if (token(1:1).eq.'.'.and.iend.eq.1) isingledot = 1
	if (token(1:1).eq.'-'.and.iend.eq.1) isingledot = 1
c
c       Add a trailing period, if necessary (assume < 1000 characters):
c
        per = ' '
	if (index(token(1:iend),'.').eq.0) per = '.'
c
        istat = 0
c
c       Trying to read ". " or "- " will get NaN and io = 0...
c
        temp = token(1:iend)//per
	if (isingledot.eq.0) then
	    read(temp(1:iend+1),'(f32.16)',iostat=io)yy
	else
	    io = 1
	end if
c
	if (io.eq.0) then
	    xread = yy
	else if (token(1:1).eq.'x'.or.token(1:1).eq.'X') then
	    xread = 1.
	else if (token(1:1).eq.'y'.or.token(1:1).eq.'Y') then
	    xread = 2.
	else if (token(1:1).eq.'z'.or.token(1:1).eq.'Z') then
	    xread = 3.
	else if (token(1:1).eq.'=') then
	else if (token(1:1).eq.'.') then
	    xread = xx
	else if (token(1:1).eq.'>'.and.iend.gt.2) then
	    read(temp(2:iend+1),'(f32.16)',iostat=io)yy
	    if (io.eq.0) xread = max(xx,yy)
            istat = io
	else if (token(1:1).eq.'<'.and.iend.gt.2) then
	    read(token(2:iend+1),'(f32.16)',iostat=io)yy
	    if (io.eq.0) xread = min(xx,yy)
            istat = io
        else
            istat = io
	end if
c
	end


	subroutine readitoken(token,iread,ii)
        save
	character*(*) token
c
c	Read an integer quantity iread from the specified token, subject
c	to the same rules as for real quantities.
c
	call readrtoken(token,xread,float(ii))
	iread = nint(xread)
c
	end
