
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


      subroutine dplot(x,y,ipen)
        save
c
c     Draw a setpat-defined dashed line to (r,s).
c
      common/dash/dpatrn(10),dpat,npatrn,ipat,lpen
      data tol/0.002/
c
      r = x
      s = y
      go to 10
c
      entry udplot(x,y,ipen)
      entry userdplot(x,y,ipen)
c
      call fr inches(x,y,r,s)
c
   10 if(abs(ipen).eq.3)then
          call plot(r,s,ipen)
          return
      end if
c
      if(kount.eq.0)then
          kount=1
          if(dpat.eq.0.)call setpat(0,0,0,0)
      end if
c
      call lastp(rl,sl)
      dr=r-rl
      ds=s-sl
      travl=sqrt(dr**2+ds**2)
      if(travl.gt.1.e6*dpat)stop 'warning: (d)line length >1.e6*dpat.'
c
      if(travl.eq.0.)then
        drdl=0.
        dsdl=0.
        go to 1
      end if
c
      drdl=dr/travl
      dsdl=ds/travl
    1 do 2 i=ipat,npatrn
      iplast=i
      step=min(travl,dpat)
      rl=rl+step*drdl
      sl=sl+step*dsdl
      call plot(rl,sl,lpen)
      travl=travl-step
      if(travl.le.tol) go to 4
      dpat=dpatrn(i+1)
    2 lpen=5-lpen
    3 dpat=dpatrn(1)
      ipat=1
      lpen=2
      if(travl.gt.tol) go to 1
      return
c
    4 ipat=iplast
      dpat=dpat-step
      if(dpat.gt.tol) return
      ipat=ipat+1
      if(ipat.gt.npatrn) go to 3
      dpat=dpatrn(ipat)
      lpen=5-lpen
c
      end
