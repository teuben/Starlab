
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
      subroutine contor(a,m,n,v,iv,xmin,xmax,ymin,ymax)
        save
c
c     Easy-to-use contour plotter.
c
      dimension a(m,n),v(iv)
      dimension x(2000),y(2000)
      external plot
c
      if(m.le.1.or.n.le.1)return
c
      do 100 i=1,m
100   x(i)=xmin+(i-1.)*(xmax-xmin)/(m-1.)
      do 200 j=1,n
200   y(j)=ymin+(j-1.)*(ymax-ymin)/(n-1.)
c
c     Split large images into strips...
c
      npass=(m*n)/65536+1
      nstep=max(3,n/npass)
c
      do 500 j=1,n,nstep-1
500   call convec(a(1,j),x,y(j),m,min(nstep,n-j+1),v,iv,plot)
c
      end


      subroutine dcontor(a,m,n,v,iv,xmin,xmax,ymin,ymax)
        save
c
c     Same as contor, but uses dplot.
c
      dimension a(m,n),v(iv)
      dimension x(2000),y(2000)
      external plot,dplot
c
      if(m.le.1.or.n.le.1)return
c
      do 100 i=1,m
100   x(i)=xmin+(i-1.)*(xmax-xmin)/(m-1.)
      do 200 j=1,n
200   y(j)=ymin+(j-1.)*(ymax-ymin)/(n-1.)
c
c     Split large images into strips...
c
      npass=(m*n)/65536+1
      nstep=max(3,n/npass)
c
      do 500 j=1,n,nstep-1
500   call convec(a(1,j),x,y(j),m,min(nstep,n-j+1),v,iv,dplot)
c
      end
