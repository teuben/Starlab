
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

      subroutine fr lnfnc(x1, x2, dx, mode, first, npost)
      save
      data eps/0.001/

c     mode = 1:  Establish fenceposts between x1 and x2, at interval dx
c     mode = 2:  Move posts out to contain x1 and x2

      if (x2.lt.x1) then
         small = x2
         big = x1
      else
         small = x1
         big = x2
      end if

c     Note attempts to deal with real to integer conversion!

      dxg = abs(dx)
      fs = small/dxg
      if (fs .gt. 0) then
         ns = fs + eps
      else
         ns = fs - eps
      endif

      if (abs(fs-ns).gt.eps) then
         if (small.gt.0.) ns = ns + 1
         if (mode.eq.2) ns = ns - 1
      end if

      fb = big/dxg
      if (fb .gt. 0) then
         nb = fb + eps
      else
         nb = fb - eps
      endif

      if (abs(fb-nb).gt.eps) then
         if (big.lt.0.) nb = nb - 1
         if (mode.eq.2) nb = nb + 1
      end if

      if (x2.le.x1) then
         i = nb
         nb = ns
         ns = i
      end if

      first = ns*dxg
      npost = abs(nb-ns) + 1

      if (mode.eq.2) then
         x1 = first
         x2 = nb*dxg
      end if

      end
