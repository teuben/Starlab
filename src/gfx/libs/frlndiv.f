
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

      subroutine fr lndiv(x1,x2,xlen,dxs,dxm,dxl)
        save
c
c     Figure out "suitable" tick spacings (neither too many nor
c     too few numbers or major ticks).
c
      dimension subt(4),subl(3)
      data subt/.1,.2,.25,.5/,subl/1.,2.,5./,ths,thl/.15,1.5/
c
      dxinch = (x2-x1)/xlen
      call compoz(dxinch,dxpi,lpow)
      fac = 10.**(lpow)
      if (x2.lt.x1) fac = -fac
      scale = 1.
c
c     Choose small tick spacing first, (try to make them 0.1 screen
c     units apart).
c
      dipx = 1./abs(dxpi)
    1 do 2 i = 1,4
          j = i
          if (subt(j)*dipx .gt. ths) go to 3
    2 continue
      dipx = dipx*10.
      scale = scale*10.
      go to 1
c
    3 xx = scale*fac
      dxs = subt(j)*xx
c
c     Now set the medium tick spacing to the next factor of ten up,
c     except that we put medium ticks at half that interval (i.e. at
c     the 0.5 marks) if dxs = 0.1.
c
      dxm = xx
      if (j.eq.1) dxm = .5*xx
c
c     Finally, choose the large tick spacing (aim for thl screen
c     units apart).
c
    4 do 5 i = 1,3
          j = i
          if (subl(j)*dipx .ge. thl) go to 6
    5 continue
      scale = scale*10.
      dipx = dipx*10.
      go to 4
c
    6 dxl = subl(j)*scale*fac
c
c     Make sure there are enough large ticks.
c
      if (abs(dxl).gt..5*abs(x2-x1)) then
	  if (j.gt.1) then
	      j = j - 1
	      dxl = subl(j)*scale*fac
	  else
	      dxl = .5*dxl
	  end if
      end if
c
      end
