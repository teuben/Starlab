
c     
c     Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c     by Steve McMillan, Drexel University, Philadelphia, PA.
c     
c     All rights reserved.
c     
c     Redistribution and use in source and binary forms are permitted
c     provided that the above copyright notice and this paragraph are
c     duplicated in all such forms and that any documentation,
c     advertising materials, and other materials related to such
c     distribution and use acknowledge that the software was developed
c     by the author named above.
c     
c     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c     IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c     WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c     


      subroutine fr labfrm(first,dx,nx,nlab,ndec,npow)
      save
c     
c     Given the (nx) values from (first) in increments of (dx)
c     to be written, find a suitable display format.
c     
c     ndec=-1   : integer format
c     ndec.gt.0 : number of places right of dec point
c     npow.ne.0 : e format
c     nlab      = total number of spaces for label
c     
c     nmax is the maximum number of characters allowed in a number
c     without moving to exponential format.
c
      data nmax /5/
c
      npow=0
c     
c     Decompose dx into mantissa f and exponent ldx.
c
      call compoz(dx,f,ldx)
c
c     Find largest label to set scale, and decompose it into f and lbig.
c
      final= first+(nx-1)*dx
      big=first
      if (abs(first).lt.abs(final)) big=final
      call compoz(big,f,lbig)
c
c     Set nsin = 1 if there are any minus signs involved.
c
      nsin=0
      if (first.lt.0. .or. final.lt.0.) nsin=1
c
c     *** Is integer format OK? (return with npow = 0, ndec = -1 if so) ***
c
      if (ldx.lt.0) go to 1
      nlab=nsin+1+lbig
c
c     Move to exponential format if there are too many digits.
c
      if (nlab.gt.nmax) go to 2
c
      ndec=-1
      return
c
c     *** Increment < 1, is F format OK? (ndec > 0, npow = 0) ***
c
    1 ndec=abs(ldx)
      nlab=nsin+1+ndec
      if (lbig.ge.0) nlab=nlab+lbig+1
      if (nlab.le.nmax) return
c
c     *** Use E format (ndec > 0, npow nonzero) ***
c
    2 nlab=nsin+5
      npow=lbig
      if (abs(npow).ge.10) nex=nex+1
      if (npow.lt.0) nlab=nlab+1
      ndec=npow-ldx
      if (ndec.lt.1) ndec=1
      nlab=nlab+ndec
c     
      end
