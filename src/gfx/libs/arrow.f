
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


      subroutine arrow(xtail,ytail,xhead,yhead,mode)
      save
c
c	Draw arrow from (x,y)-tail to (x,y)-head.
c	If mode = 0, x and y are user units.
c	If mode = 1, x and y are in inches.
c
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      data open,pi/.5236,3.141593/head_size/-1.0/

      if(mode.eq.0)then
          rtail=(xtail-xl)*dinchx
          stail=(ytail-ybot)*dinchy
          dr=(xhead-xtail)*dinchx
          ds=(yhead-ytail)*dinchy
      else
          rtail=xtail
          stail=ytail
          dr=xhead-xtail
          ds=yhead-ytail
      end if

c     Draw the shaft of the arrow.

      theta=atan2(ds,dr)
      call plot(rtail,stail,3)

      rhd=rtail+dr
      shd=stail+ds
      call plot(rhd,shd,2)

c     Draw the point, scaled to shaft size if not otherwise specified.

      alpha=theta+pi-open
      shaft=sqrt(dr*dr+ds*ds)

      if (head_size.le.0.0) then
          if(shaft.gt.1.)then
              point=.15*shaft
          else if(shaft.lt..5)then
              point=max(.05,.2*shaft)
          else
              point = 2.*(shaft-.5)*.15*shaft
     &                + 2.*(1.-shaft)*max(.05,.2*shaft)
          end if
      else
          point = head_size
      end if

      r=rhd+point*cos(alpha)
      s=shd+point*sin(alpha)
      call plot(r,s,2)
      call plot(rhd,shd,3)
      alpha=alpha+2.*open
      r=rhd+point*cos(alpha)
      s=shd+point*sin(alpha)
      call plot(r,s,2)

      return

      entry set_arrow_head(size)

c     Specify arrowhead size, in "inches".

      head_size = size

      end
