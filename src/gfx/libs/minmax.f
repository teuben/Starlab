
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


      subroutine minmax(a,n,amin,amax)
        save
      dimension a(1)
c
      amin = 1.e30
      amax = -1.e30
      call minmax1(a,n,amin,amax)
c
      end


      subroutine minmax1(a,n,amin,amax)
        save
c
c     Like minmax, but use the input amin, amax as initial values.
c
      dimension a(1)
c
      do 100 i=1,n
          amin = min(amin,a(i))
          amax = max(amax,a(i))
100   continue
c
      end


      subroutine iminmax(a,n,imin,imax)
        save
c
c     Like minmax, but return the locations of the minimum and the maximum.
c
      dimension a(1)
c
      imin = 0
      imax = 0
      amin=1.e30
      amax=-1.e30
      do 100 i=1,n
	  if (a(i).lt.amin) then
	      amin = a(i)
	      imin = i
	  end if
	  if (a(i).gt.amax) then
	      amax = a(i)
	      imax = i
	  end if
100   continue
c
      end
