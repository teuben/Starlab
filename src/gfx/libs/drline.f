
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



      subroutine drline (z,xx,yy,mm,nn,userplot)
        save
      external userplot
      dimension       z(mm,nn),xx(mm),yy(nn)
c
c this routine traces a contour line when given the beginning by contor1.
c
      common /cinfo/ ix,iy,idx,idy,is,iss,cv,inx(8),iny(8) 
      common icont(65536)
c
      c(p1,p2) = (p1-cv)/(p1-p2)
c
      m = mm
      n = nn
c
c     store initial info for termination criterion.
c
      ix0 = ix
      iy0 = iy
      is0 = is
c
c     go to first point.
c
      if (idx.ne.0) then
          y = yy(iy)
          isub = ix+idx
          x = c(z(ix,iy),z(isub,iy))*(xx(isub)-xx(ix))+xx(ix)
      else
          x = xx(ix)
          isub = iy+idy
          y = c(z(ix,iy),z(ix,isub))*(yy(iy+idy)-yy(iy))+yy(iy)
      end if
      call userplot (x,y,3)
c
c     look for next crossing.
c
  106 is = is+1
      if (is .gt. 8) is = is-8
      idx = inx(is)
      idy = iny(is)
      ix2 = ix+idx
      iy2 = iy+idy
      if (iss .ne. 0) go to 107
      if (ix2.gt.m.or.iy2.gt.n.or.ix2.lt.1.or.iy2.lt.1) go to 120
  107 if (cv.gt.z(ix2,iy2)) go to 109
  108 is = is+4
      ix = ix2
      iy = iy2
      go to 106
  109 if (is/2*2 .eq. is) go to 106
c
c     draw next contour segment.
c
      if (idx.ne.0) then
          y = yy(iy)
          isub = ix+idx
          x = c(z(ix,iy),z(isub,iy))*(xx(isub)-xx(ix))+xx(ix)
      else
          x = xx(ix)
          isub = iy+idy
          y = c(z(ix,iy),z(ix,isub))*(yy(isub)-yy(iy))+yy(iy)
      end if
      call userplot (x,y,2)
      if (is .eq. 1) icont(ix-1+m*(iy-1))=1
c
c     mark presence of contour in this cell.
c
      if (iss .eq. 0) go to 106
      if (ix.ne.ix0 .or. iy.ne.iy0 .or. is.ne.is0) go to 106
c
c     end of line
c
  120 return
      end
