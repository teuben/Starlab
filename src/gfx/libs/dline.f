
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

      subroutine dline(xarray,yarray,n,jth,jsymbl,htsym)
      save
c     
c     Like mline, but draws a dashed line.
c
c     ***** This appears to have been absorbed into mline!
c     
      dimension xarray(1),yarray(1)
      character sim*3
      common/scales/ xminim,xmax,dxinch,yminim,ymax,dyinch,rlen,slen
      common/dash/dpatrn(10),dpat,npatrn,ipat,lpen
      common/fr bnds/mode
c     
      common /ngon stars/ istar
c     
      cinch(x,x0,dxi)=(x-x0)*dxi
c     
      call routine id('dline')
c     
      if(jth.eq.0) go to 2
      hite=abs(htsym)
      if(jsymbl.lt.0)then
          iascii=-jsymbl
          if(iascii.le.127)then
              nsym=-1
              sim=char(iascii)
          else if(iascii.le.223)then
              nsym=-2
              sim='@'//char(iascii-96)
          else
              nsym=-2
              sim='%'//char(iascii-192)
          end if
          if(htsym.lt.0.)nsym=-nsym
      end if
c     
c     *** plot symbols ***
c     
      int=abs(jth)
      do i=1,n,int
          x1=cinch(xarray(i),xminim,dxinch)
          y1=cinch(yarray(i),yminim,dyinch)
          if(mode.ne.0.and.
     +            (x1.lt.0..or.x1.gt.rlen.or.y1.lt.0..or.y1.gt.slen))go to 1
          if(jsymbl.ge.0)then
              if(istar.eq.0)then
                  call ngon(x1,y1,.5*htsym,jsymbl,0.)
              else
                  call ngon(x1,y1,.5*htsym,-jsymbl,0.)
              end if
          else
              if (iascii.le.1) then
c                 
c                 Encode the number of the point in sim.
c                 
                  ii = i
                  do while (ii.gt.61)
                      ii = ii - 61
                  end do
                  if (ii.le.9) then
                      sim = char(48+ii)
                  else if (ii.le.35) then
                      sim = char(87+ii)
                  else
                      sim = char(55+ii)
                  end if
              end if
              call simbol(x1,y1,hite,sim,0.,nsym)
          end if
      end do
      if(jth.lt.0) return
c     
c     *** plot line ***
c     
    2 x1=cinch(xarray(1),xminim,dxinch)
      y1=cinch(yarray(1),yminim,dyinch)
      if(mode.eq.0)then
          call plot(x1,y1,3)
      else
          call plotin(x1,y1,3)
      end if
      do i=2,n
          x1=cinch(xarray(i),xminim,dxinch)
          y1=cinch(yarray(i),yminim,dyinch)
          if(mode.eq.0)then
              call dplot(x1,y1,2)
          else
              call dplotin(x1,y1,2)
          end if
      end do
c
      end
