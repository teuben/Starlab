
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
        subroutine mc sym lims
        save
c
c       Return the correct /fontc1/ info for symbl.
c
        real*4 lastinc
        common /fontc1/ offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
        common /numsym int/ rs,ss,dx,dy,sint,cost
        dimension x(4),y(4)
        x(1)=0.
        y(1)=0.
        x(2)=dx*cost
        y(2)=dx*sint
        x(4)=-dy*sint
        y(4)=dy*cost
        x(3)=x(2)+x(4)
        y(3)=y(2)+y(4)
        xmin=1.e6
        ymin=xmin
        xmax=-xmin
        ymax=xmax
        do 10 i=1,4
            if(x(i).gt.xmax)xmax=x(i)
            if(x(i).lt.xmin)xmin=x(i)
            if(y(i).gt.ymax)ymax=y(i)
            if(y(i).lt.ymin)ymin=y(i)
10      continue
        xmin=xmin+rs
        xmax=xmax+rs
        ymin=ymin+ss
        ymax=ymax+ss
        xp=rs+dx*cost
        yp=ss+dx*sint
        end
