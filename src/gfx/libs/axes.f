
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


      subroutine axes(yax,xax,id)
        save
c
c     This routine draws an x-axis (y=yax) and/or a y-axis (x=xax).
c     It will not draw either axis outside the boundaries defined by
c     'frame'.  Axis is dashed for -ve id.
c
      dimension xa(2),ya(2)
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
c
      iboth=1
      go to 10
c
      entry xaxis(yax,id)
      iboth=-1
10    if(sign(1.,(ybot-yax)).eq.sign(1.,(ytop-yax)))go to 20
      xa(1)=xl
      ya(1)=yax
      xa(2)=xr
      ya(2)=yax
      if(id.ge.0)call mline(xa,ya,2,0,0,0.)
      if(id.lt.0)call dline(xa,ya,2,0,0,0.)
20    if(iboth.lt.0)return
      go to 25
c
      entry yaxis(xax,id)
25    if(sign(1.,(xl-xax)).eq.sign(1.,(xr-xax)))go to 30
      xa(1)=xax
      ya(1)=ybot
      xa(2)=xax
      ya(2)=ytop
      if(id.ge.0)call mline(xa,ya,2,0,0,0.)
      if(id.lt.0)call dline(xa,ya,2,0,0,0.)
c
30    return
      end
