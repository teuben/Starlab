
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

        subroutine outp(x,y,ier)
        save
c
c	Decide on the visibility of point (x,y), and maintain the
c	"horizon" array.  This routine is used by nxtvu.
c
        parameter (nn = 2000, eps = 0.001)
        common /nxtv1/ xx(nn),yy(nn),kk,ll
c
        if (kk.eq.0) go to 10
        if (kk.eq.ll-1) go to 20
        if (abs(xx(kk)-x)+abs(yy(kk)-y).lt.eps) return
c
10      kk = kk+1
        xx(kk) = x
        yy(kk) = y
        return
c
20      ier = 1
        end
