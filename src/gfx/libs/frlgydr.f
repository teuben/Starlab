
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

      subroutine fr lgydr(x,y1,y2,dexps,dexpm,dexpl,iax,ilab)
      save
c     
c     As for logxdr, but for the y-axis.
c     
      dimension tiks(8)
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/dev status/idevon,idevpen,idevwt
      common/fr tik level/jtik level
      common/fr hts/htl,htn/fr wts/iwts(4)/fr ticks/tikk(3),tikl
      common/fr setax/kax,lax
c     
      data ntik,tiks/8,.301003,.4771213,.60206,.69897,.7781513,
     *        .845098,.90309,.9542425/
c     
      cinch(dumx,dor,dinch)=(dumx-dor)*dinch
c     
      if (lax.lt.0) return
c     
      nlab=3
      ndec=0
      lpow = 1
      call fr lnydr(x,y1,y2,dexps,dexpm,dexpl,
     $              iax,ilab,nlab,ndec,lpow)
c     
      if (abs(dexpm).ne.1.) return
      if (abs(dinchy).le..5) return
      if (jtik level.ne.1) return
      if (lax.eq.1.and.iax.eq.2) return
      if (lax.eq.2.and.iax.eq.1) return
c     
c     Add markers for integers.
c     
      dr=1.
      if (iax.eq.2) dr=-dr
      rax=cinch(x,xl,dinchx)
      rtik=rax+dr*tikl
      ya=y1
      yb=y2
      if (yb.le.ya) then
          ya=y2
          yb=y1
      end if
c     
      call fr lnfnc(ya,yb,1.,1,fexp,nexp)
c     
      fexp=fexp-1.
      nexp=nexp+1
      if (iwts(1).gt.0)then
          jwt=idevwt
          call weight(iwts(1))
      end if
c     
c     Add logarithmically spaced tick marks.
c     
      do i=1,nexp
          do j=1,ntik
              exp=fexp+tiks(j)
              if (exp.ge.ya) then
                  if (exp.gt.yb) go to 100
                  s=cinch(exp,ybot,dinchy)
                  call plot(rax,s,3)
                  rt = rtik
c
c                 Note test here is 4, not 3 because nexp was incremented above.
c
                  if (nexp.le.4.and.j.eq.2) rt = rax+dr*tikk(3)
                  call plot(rt,s,2)
              end if
          end do
          fexp=fexp+1.
      end do
c
100   if (iwts(1).gt.0) call weight(jwt)
c     
      end
