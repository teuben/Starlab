
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


      subroutine fr lgset(x1,x2,xlen,mode,dexps,dexpm,dexpl)
      save
c
      common/fr tik level/jtik level
c     
c     Find increments for small, major and labeled ticks for log plot.
c     
      x1=alog10(abs(x1))
      x2=alog10(abs(x2))
      call fr lndiv(x1,x2,xlen,dexps,dexpm,dexpl)
      if (abs(dexps).lt.1.) dexps=0.
      if (abs(dexpm).lt.1.) dexpm=dexpm/abs(dexpm)
      if (abs(dexpl).lt.1.) dexpl=dexpm
      if (mode.ne.-2) return
c
      dfnc=dexpm
      if (dexps.ne.0.) dfnc=dexps
      call fr lnfnc(x1,x2,dfnc,2,fexp,nexp)
c
c     dinch=xlen/abs(x2-x1)
c     if (jtik level.ne.1.or.dinch.le..5.or.abs(dexpm).ne.1.) then
c     
c     As in frlnset, force at least two numeric labels if mode = 2.
c     
      call fr lnfnc(x1,x2,abs(dexpl),1,fexp,nexp)
      if (fexp+abs(dexpl).gt.max(x1,x2)) dexpl=2.*dexpm
c
c     end if
c
      end
