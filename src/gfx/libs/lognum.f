
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


      subroutine lognum(xi,yi,ht1,expt,theta,ll)
      save
c     
c     Plots the number 10**(expt)
c     Leaves xi,yi unchanged.
c     For integer values of the exponent, set ll = -1
c     
      common/dev status/idevon,idevpen,idevwt
      common/fr int/iframe/fr wts/iwts(4)/fr plain/iplain/fr bare/ibare
      common/fr pens/ipens(3)
c     
      if((iwts(2).lt.0.or.iwts(3).lt.0).and.iframe.eq.1)return
c
      if(iframe.eq.1.and.ipens(2).gt.0)then
          jpen=idevpen
          call pen(ipens(2))
      end if
c
      if(iplain.eq.0)then
          call simbol(xi,yi,ht1,'10',theta,2)
      else
          if(ibare.eq.0)then
              call simbol(xi,yi,ht1,'@1@0',theta,4)
          else
              jfr=iframe
              iframe=0
              call fr numbr(xi,yi,ht1,10.,theta,-1)
              iframe=jfr
          end if
      end if
c
      call simwhe(xp,yp)
      ht2=0.53333333*ht1
      thetr=0.01745329*theta
      xe=xp-.75*ht1*sin(thetr)
      ye=yp+.75*ht1*cos(thetr)
      if(iwts(3).gt.0.and.iframe.eq.1)then
          jwt=idevwt
          call weight(iwts(3))
      end if
      jframe=iframe
      iframe=0                  ! don't want centered symbol.
      call fr numbr(xe,ye,ht2,expt,theta,ll)
      iframe=jframe
      if(iwts(3).gt.0.and.iframe.eq.1)call weight(jwt)
      if(ipens(2).gt.0.and.iframe.eq.1)call pen(jpen)
c
      end
