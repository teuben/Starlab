
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


      subroutine aerrbars(xc,yc,dx,dy)
        save
c
c	Draw one-sided error bars from (xc,yc), of length(s)
c	dx and/or dy, all measured in user units.
c
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      data cap/.07/
c
      iboth=1
      go to 10
c
      entry axerr(xc,yc,dx)
c
      iboth=-1
10    rc=(xc-xl)*dinchx
      sc=(yc-ybot)*dinchy
      dr=dx*dinchx
c
c      *** horizontal bar ***
c
      call plotin(rc,sc,3)
      r=rc+dr
      call plotin(r,sc,2)
      call plotin(r,sc-cap,3)
      call plotin(r,sc+cap,2)
      if(iboth.lt.0)return
      go to 20
c
      entry ayerr(xc,yc,dy)
c
      rc=(xc-xl)*dinchx
      sc=(yc-ybot)*dinchy
20    ds=dy*dinchy
c
c      *** vertical bar ***
c
      call plotin(rc,sc,3)
      s=sc+ds
      call plotin(rc,s,2)
      call plotin(rc-cap,s,3)
      call plotin(rc+cap,s,2)
c
      end
