
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
c---	subroutine:  spiral
c
c---	draws logarithmic spirals
c		r	radius of spiral in inches
c		a	number of spirals per inch
c		t	angle of starting curve (degrees cc from x-axis)
c		i	direction of spiral (+1 = ccl, -1 = cl)
c		x,y	location of center
c
        subroutine spiral (r,a,t,i,x,y)
        save
        pi=3.1415926
c
        c1=.02*pi
        c2=2.*pi*a
        c3=2.*pi*t/360.
        call plot(x,y,3)
c
        do 10 j=1,100*a*r
            th=c1*j
            ra=th/c2
            th=c3+i*th
            x1=x+ra*cos(th)
            y1=y+ra*sin(th)
            call plot(x1,y1,2)
10      continue
        end
