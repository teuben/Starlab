
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
c
c-----------------------------------------------------------------------
c
        subroutine display text(text,ntext)
        save
        character*(*) text
c							    ******************
c	Print a message in standard form.		    *  DISPLAY TEXT  *
c							    ******************
        if(ntext.le.0)return
	call devoff
c
        if(itek.eq.1)then
            rr=rl
            ss=sl
            call plot(-ro,-so,3)
        end if
c
	lt=min(len(text),ntext)
        write(6,'(a)')text(1:lt)
        if(itek.eq.1)call plot(rr,ss,3)
c
        end
