
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
c
        subroutine heart (x0,y0,r,t)
        save
c
        pi=3.1415926
        ct=cos(pi*t/180.)
        st=sin(pi*t/180.)
        ax=(180.+28.072487)*pi/180.
c
c---	right hump
        call plot(x0,y0,3)
        do 5 i=1,100
            a=i*.01*ax
            x=.25*r*(1.-cos(a))
            y=.25*r*sin(a)
            xx=x0+x*ct+y*st
            yy=y0+y*ct-x*st
            call plot(xx,yy,2)
5       continue
        call plot(x0-r*st,y0-r*ct,2)
c
c---	left hump
        call plot(x0,y0,3)
        do 10 i=1,100
            a=i*.01*ax
            x=.25*r*(cos(a)-1.)
            y=.25*r*sin(a)
            xx=x0+x*ct+y*st
            yy=y0+y*ct-x*st
            call plot(xx,yy,2)
10      continue
        call plot(x0-r*st,y0-r*ct,2)
c
        end
