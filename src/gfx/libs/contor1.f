
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


      subroutine contor1 (z,x,y,mm,nn,conv,userplot)
        save
      external userplot
      dimension       z(mm,nn),x(mm),y(nn)
c
c     this routine finds the beginnings of all contour lines at level conv.
c     first the edges are searched for lines intersecting the edge (open
c     lines) then the interior is searched for lines which do not intersect
c     the edge (closed lines).  beginnings are checked against icont to prevent
c     re-tracing of lines.
c
c     subroutine userplot converts from (x,y) coordinates to
c     the desired rectangular output grid. if the original grid
c     is already rectangular and correctly scaled (see mcinit),
c     just set userplot = plot in the call to contor.
c
      common /cinfo/ ix,iy,idx,idy,is,iss,cv,inx(8),iny(8) 
      common icont(65536)
      data inx(1),inx(2),inx(3),inx(4),inx(5),inx(6),inx(7),inx(8)/
     1        -1 ,   -1 ,    0 ,    1 ,    1 ,    1 ,    0 ,   -1 /
      data iny(1),iny(2),iny(3),iny(4),iny(5),iny(6),iny(7),iny(8)/
     1         0 ,    1 ,    1 ,    1 ,    0 ,   -1 ,   -1 ,   -1 /
c
      do 1 i=1,mm*nn+1
1     icont(i)=0
      l = ll
      m = mm
      n = nn
      cv = conv
      iss = 0
      do 102 ip1=2,m
         i = ip1-1
         if (z(i,1).ge.cv .or. z(ip1,1).lt.cv) go to 101
         ix = ip1
         iy = 1
         idx = -1
         idy = 0
         is = 1
         call drline (z,x,y,m,n,userplot)
  101    if (z(ip1,n).ge.cv .or. z(i,n).lt.cv) go to 102
         ix = i
         iy = n
         idx = 1
         idy = 0
         is = 5
         call drline (z,x,y,m,n,userplot)
  102 continue
      do 104 jp1=2,n
         j = jp1-1
         if (z(m,j).ge.cv .or. z(m,jp1).lt.cv) go to 103
         ix = m
         iy = jp1
         idx = 0
         idy = -1
         is = 7
         call drline (z,x,y,m,n,userplot)
  103    if (z(1,jp1).ge.cv .or. z(1,j).lt.cv) go to 104
         ix = 1
         iy = j
         idx = 0
         idy = 1
         is = 3
         call drline (z,x,y,m,n,userplot)
  104 continue
      iss = 1
      do 108 jp1=3,n
         j = jp1-1
         do 107 ip1=2,m
            i = ip1-1
            if (z(i,j).ge.cv .or. z(ip1,j).lt.cv) go to 107
            if(icont(i+m*(j-1)).ne.0)go to 107
            icont(i+m*(j-1))=1
            ix = ip1
            iy = j
            idx = -1
            idy = 0
            is = 1
c
c           note: start contour only if
c
c                 z(ix-1,iy) .lt. cv .le. z(ix,iy)
c
c           all closed contours must have this ordering somewhere
c           in the grid.
c
            call drline (z,x,y,m,n,userplot)
  107    continue
  108 continue
      end
