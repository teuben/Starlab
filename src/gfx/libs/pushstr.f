
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

        subroutine pushstr
        save
c
c       Push current strpos settings on a stack. Don't alter settings.
c
        dimension stack(30)
        common/str posn/i pos set,frx,fry
c
        data ipoint/0/
c
        if(ipoint.gt.27)stop'str stack overflow.'
        ipoint=ipoint+1
        stack(ipoint)=i pos set
        ipoint=ipoint+1
        stack(ipoint)=frx
        ipoint=ipoint+1
        stack(ipoint)=fry
        return
c
        entry popstr
c
c       pop last strpos settings from the stack.
c
        if(ipoint.lt.2)stop'str stack underflow.'
        fry=stack(ipoint)
        ipoint=ipoint-1
        frx=stack(ipoint)
        ipoint=ipoint-1
        i pos set=stack(ipoint)
        ipoint=ipoint-1
c
        end
