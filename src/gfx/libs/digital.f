
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
c---	subroutine digital
c
c---	creates a digital map of an array
c
c---	created  12.aug.80
c		 27.mar.81  for xplot
c
c---	parameters
c	a	- data array
c	m,n	- dimensions of a to plot
c	isc	- icode for nomber
c	x0,x1,y0,y1 - delimeters of plot
c	msiz,nsiz - dimensions of array a
c
        subroutine digital (a,m,n,isc,x0,x1,y0,y1,msiz,nsiz)
        save
        real*4 a(msiz,nsiz)
        common/fontc1/offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
        dx=(x1-x0)/m
        dy=(y1-y0)/n
        ht=dy/2.
        do 10 j=1,n
            do 10 i=1,m
                il=alog10(abs(a(i,j)))+isc+1
                if (a(i,j).lt.0.) il=il+1
                h1=.5*dx/il
                ht=min(h1,ht)
10      continue
        do 20 j=1,n
            do 20 i=1,m
                call nomber (0.,0.,-ht,a(i,j),0.,isc)
                call nomber (x0+(i-.5)*dx-.5*(xmax-xmin),
     +                       y0+(j-.5)*dy-.5*(ymax-ymin)+.05,
     +                       ht,a(i,j),0.,isc)
20      continue
        end
