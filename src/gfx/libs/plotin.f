
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

      subroutine plotin(rcall,scall,ipen)
        save
c
c     Move pen from the last point to a new point, drawing only that
c     segment of the line lying within the "eframe" boundaries.
c
      common /scales/ x1,x2,dix,y1,y2,diy,rlen,slen
      common /fr sord/ idash
      common /last point/ rl,sl
c
c     Note: a "plotin" session should start with a "plotin(.,.,3)"
c	    call (to initialize il, etc.), not with "plot(.,.,3)"
c
      integer MOVE, DRAW
      parameter (MOVE = 3, DRAW = 2)
c
      save il,jl,kl,inl
c
c-----------------------------------------------------------------------------
c
c     Alternate entry points:
c     ----------------------
c
      idash = 0
      r = rcall
      s = scall
      go to 1
c
      entry uplotin(rcall,scall,ipen)
      entry userplotin(rcall,scall,ipen)
      idash = 0
      call fr inches(rcall,scall,r,s)
      go to 1
c
      entry dplotin(rcall,scall,ipen)
      idash = 1
      r = rcall
      s = scall
      go to 1
c
      entry udplotin(rcall,scall,ipen)
      entry userdplotin(rcall,scall,ipen)
      idash = 1
      call fr inches(rcall,scall,r,s)
c
c-----------------------------------------------------------------------------
c
c     Locate the ends of the line relative to the current box:
c
c			!		!
c	  (1,3) --> 7	!  (2,3) --> 8	!  (3,3) --> 9
c			!		!
c	- - - - - - - - +---------------+ - - - - - - - -
c			|		|
c			|  EFRAME box	|
c	  (1,2) --> 4	|		|  (3,2) --> 6
c			|  (2,2) --> 5	|
c			|		|
c	- - - - - - - - +---------------+ - - - - - - - -
c			!		!
c	  (1,1) --> 1	!  (1,2) --> 2	!  (1,3) --> 3
c			!		!
c
    1 in = 0
      i = 2
      j = 2
      if (r.lt.0.)   i = 1
      if (r.gt.rlen) i = 3
      if (s.lt.0.)   j = 1
      if (s.gt.slen) j = 3
c
      k = i + 3*(j-1)
      if (k.eq.5) in = 1
      if (ipen.eq.MOVE) go to 99
c
      if (r.eq.rl.and.s.eq.sl) then
          if (in.eq.1) call plot(r,s,DRAW)
          go to 99
      end if
c
      rp = r
      sp = s
c
c     *** both in ***
c
      if (in+inl.eq.2) go to 98
c
c     *** both out ***
c
      if (in+inl.eq.0) go to 20
c
c     *** 1 in, 1 out ***
c
      if (in.eq.1) go to 3
c
c     Last point in.
c
      rin = rl
      sin = sl
    2 rout = r
      sout = s
      iout = i
      kout = k
      go to 5
c
c     Current point in.
c
    3 rin = r
      sin = s
    4 rout = rl
      sout = sl
      iout = il
      kout = kl
c
c     Move from (rin,sin) to (rout,sout).
c
    5 call fr plt(rin,sin,MOVE)
      if (kout.eq.2) go to 11
      if (kout.eq.8) go to 12
      if (kout.eq.6) go to 14
      if (kout.eq.4) go to 15
c
c     Exit to corner square; find exit side.
c
      rp = 0.
      if (iout.eq.3) rp = rlen
      sp = sin + (rp-rin)*(sout-sin)/(rout-rin)
      if (sp.lt.0.) go to 11
      if (sp.gt.slen) go to 12
      go to 98
c
c     Exit side known...
c
c     ...through bottom.
c
   11 sp = 0.
      go to 13
c
c     ...through top.
c
   12 sp = slen
   13 rp = rin + (sp-sin)*(rout-rin)/(sout-sin)
      go to 98
c
c     ...through right.
c
   14 rp = rlen
      go to 16
c
c     ...through left.
c
   15 rp = 0.
   16 sp = sin + (rp-rin)*(sout-sin)/(rout-rin)
      go to 98
c
c     Both points outside -- look for any intersection with the frame.
c
c     First dispose of points in same square, or in same outside row or column.
c
   20 if (k.eq.kl)   go to 99
      if (j+jl.eq.2) go to 99
      if (j+jl.eq.6) go to 99
      if (i+il.eq.2) go to 99
      if (i+il.eq.6) go to 99
c
c     Is any square on a side of the frame?
c
      kout = kl
      if (k.eq.2) go to 21
      if (k.eq.8) go to 22
      if (k.eq.6) go to 24
      if (k.eq.4) go to 25
c
      kout = k
      if (kl.eq.2) go to 21
      if (kl.eq.8) go to 22
      if (kl.eq.6) go to 24
      if (kl.eq.4) go to 25
      go to 27
c
c     Through bottom.
c
   21 sin = 0.
      go to 23
c
c     Through top.
c
   22 sin = slen
   23 rin = r + (sin-s)*(rl-r)/(sl-s)
      if (rin.le.0.) go to 99
      if (rin.ge.rlen) go to 99
      if (kout.eq.kl) go to 4
      go to 2
c
c     Through right side.
c
   24 rin = rlen
      go to 26
   25 rin = 0.
   26 sin = s + (rin-r)*(sl-s)/(rl-r)
      if (sin.le.0.) go to 99
      if (sin.ge.slen) go to 99
      if (kout.eq.kl) go to 4
      go to 2
c
c     Connect diagonally opposite corner squares.
c
   27 ip = MOVE
      sin = 0.
      rin = r + (sin-s)*(rl-r)/(sl-s)
      if (rin.le.rlen.and.rin.ge.0.) then
          call fr plt(rin,sin,ip)
          ip = DRAW
      end if
c
      sin = slen
      rin = r + (sin-s)*(rl-r)/(sl-s)
      if (rin.le.rlen.and.rin.ge.0.) then
          call fr plt(rin,sin,ip)
          if (ip.eq.DRAW) go to 99
          ip = DRAW
      end if
c
      rin = 0.
      sin = s + (rin-r)*(sl-s)/(rl-r)
      if (sin.le.slen.and.sin.ge.0.) then
          call fr plt(rin,sin,ip)
          if (ip.eq.DRAW) go to 99
          ip = DRAW
      end if
c
      rin = rlen
      sin = s + (rin-r)*(sl-s)/(rl-r)
      if (sin.le.slen.and.sin.ge.0.) call fr plt(rin,sin,ip)
      go to 99
c
c     Plot final segment, then terminate on last point.
c
   98 call fr plt(rp,sp,DRAW)
   99 call fr plt(r,s,MOVE)
      rl = r
      sl = s
      il = i
      jl = j
      kl = k
      inl = in
c
      end
