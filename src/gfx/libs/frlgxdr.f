
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


      subroutine fr lgxdr(y,x1,x2,dexps,dexpm,dexpl,iax,ilab)
      save
c     
c     Draw x-axis, tick marks and numbers for a log plot.
c     (Numbers are handled by frlnxdr with ndec = 0.)
c     
      dimension tiks(8)
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/dev status/idevon,idevpen,idevwt
      common/fr tik level/jtik level
      common/fr hts/htl,htn/fr wts/iwts(4)/fr ticks/tikk(3),tikl
c
      data ntik,tiks/8,.301003,.4771213,.60206,.69897,.7781513,
     *        .845098,.90309,.9542425/
c     
      cinch(dumx,dor,dinch)=(dumx-dor)*dinch
c
c     Plot the numbers and major tick marks:
c
      nlab=2
      ndec=0
      lpow = 1
      call fr lnxdr(y,x1,x2,dexps,dexpm,dexpl,
     $              iax,ilab,nlab,ndec,lpow)
c
      if (abs(dexpm).ne.1.) return
      if (abs(dinchx).le..5) return
      if (jtik level.ne.1) return
c     
c     Add markers for integers
c     
      ds=1.
      if (iax.eq.2) ds=-ds
      sax=cinch(y,ybot,dinchy)
      stik=sax+ds*tikl
      xa=x1
      xb=x2
      if (xb.le.xa) then
          xa=x2
          xb=x1
      end if
c
      call fr lnfnc(xa,xb,1.,1,fexp,nexp)
c
      fexp=fexp-1.
      nexp=nexp+1
      if (iwts(1).gt.0) then
          jwt=idevwt
          call weight(iwts(1))
      end if
c
c     Add logarithmically spaced tick marks.
c
      do i=1,nexp
          do j=1,ntik
              exp=fexp+tiks(j)
              if (exp.ge.xa) then
                  if (exp.gt.xb) go to 100
                  r=cinch(exp,xl,dinchx)
                  call plot(r,sax,3)
                  st = stik
c
c                 Note test here is 4, not 3 because nexp was incremented above.
c
                  if (nexp.le.4.and.j.eq.2) st = sax+ds*tikk(3)
                  call plot(r,st,2)
              end if
          end do
          fexp=fexp+1.
      end do
c
100   if (iwts(1).gt.0) call weight(jwt)
c
      end
