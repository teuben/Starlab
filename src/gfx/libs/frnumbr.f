
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


      subroutine fr numbr(r,s,h,x,t,n)
        save
      common/fr plain/iplain/fr bare/ibare
      common/fr wts/iwts(4)/fr pens/ipens(3)
      common/dev status/idevon,idevpen,idevwt
      integer*4 n
      common/debug trace/itrace
c
      if(iwts(2).lt.0.or.iwts(3).lt.0)return
      if(itrace.eq.1)write(2,*)'fr numbr:',r,s,h,x
      if(ipens(2).gt.0)then
          jpen=idevpen
          call pen(ipens(2))
      end if
      if(ibare.eq.1)then
          call numbr(r,s,h,x,t,n)
      else
          if(iplain.eq.0)then
              call nomber(r,s,h,x,t,n)
          else
              call nombr(r,s,h,x,t,n)
          end if
      end if
      if(ipens(2).gt.0)call pen(jpen)
      end
