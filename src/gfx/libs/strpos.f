
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

        subroutine strpos(xfr,yfr)
        save
c
        common/str posn/i pos set,frx,fry
        data i pos set/0/frx/0./fry/0./
c
        i pos set=1
        frx=xfr
        fry=yfr
c
        return
c
	entry getpos(xfr,yfr)
	entry getstr(xfr,yfr)
c
	xfr=frx
	yfr=fry
c
        end
