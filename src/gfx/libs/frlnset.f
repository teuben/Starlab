
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


      subroutine fr lnset(x1,x2,xlen,mode,dxs,dxm,dxl,
     $                    labsp,labdp,labpow)
      save
c     
c     Determines increments for small, major and numbered ticks
c     for a linear plot
c     
      call fr lndiv(x1,x2,xlen,dxs,dxm,dxl)
      call fr lnfnc(x1,x2,dxm,mode,firstl,nl)
      call fr lnfnc(x1,x2,dxl,1,firstl,nl)
c     
      if (mode.eq.2) then
c         
c         Insist on at least two numeric labels.
c         
          if(firstl+dxl.gt.max(x1,x2))then
              dxl=2*dxm
              if(firstl-dxl.ge.min(x1,x2))firstl=firstl-dxl
          end if
      end if
c
      call fr labfrm(firstl,dxl,nl,labsp,labdp,labpow)
c     
      end
