
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

	program energies
c
	parameter (XLEN = 5., YLEN = 5., YMAX = 12.)
	parameter (rc = 1., rh = 50., rt = 500.)
c
	parameter (N = 17, ALPHA = .4)
	dimension e(0:N+1),decm(N)
c
	dimension phi(500),rapo(500)
c
	e(0) = 0.
	e(1) = 1.
	do 10 i=2,N+1
10	e(i) = e(i-1)*(1.+ALPHA)
c
	do 20 i=1,N
20	decm(i) = (e(i+1) - e(i))/6.
c
	ddphi = .05
	dphi = -ddphi
	nph = 0
c
30	dphi = dphi + ddphi
	nph = nph + 1
	phi(nph) = dphi
	rapo(nph) = log10(apocenter(rc,rh,rt,dphi))
	if (rapo(nph).lt.log10(rt).and.nph.lt.500) go to 30
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	call mcinit
	call nobounds
	call sethts(.25,.2)
	call plot(4.25,1.,-3)
c
	call strpos(1.,1.)
	call simbol(XLEN,YLEN+1.2,.15,'Figure 1%%',0.,999)
	call clrstr
c
	call weight(10)
	call eframe(rc,rt,XLEN,-1,'!1.1!r\\a',0.,YMAX,YLEN,2,' ')
c
	call simsize(.25,'@D@f',4,dr1,ds)
	call simsize(.25,'@s^2',4,dr2,ds)	
	dr = max(dr1,dr2)
	s = 7.5*YLEN/YMAX
	call plot(-.25-dr,s,3)
	call plot(-.25,s,2)
c
	call strpos(1.,0.)
	call simbol(-.25,s+.05,.25,'@D@f',0.,4)
	call strpos(1.,1.)
	call simbol(-.25,s-.05,.25,'@s^2',0.,4)
	call clrstr
	call mline(rapo,phi,nph,0,0,0.)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	do 60 i=7,N
	    if (i.eq.8) go to 60
	    s = decm(i)*YLEN/YMAX
	    if (s.le.YLEN) then
	        call dplot(0.,s,3)
	        do 55 j=1,nph
		    if (phi(j).gt.decm(i)) then
		        r = rapo(j-1) + (rapo(j)-rapo(j-1))
     &					    *(decm(i)-phi(j-1))
     &					        /(phi(j)-phi(j-1))
		        call frinches(r,ss,r,ss)
		        go to 56
		    end if
55	        continue
56	        call setpat(1,2,1,2)
		call dplot(r,s,2)
	        call setpat(1,2,1,2)
	        call dplot(r,0.,3)
	        call dplot(r,s,2)
	    end if
60	continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	call plot(-2.25,0.,-3)
	call strpos(0.,.5)
	do 40 i=1,N
	    s = decm(i)*YLEN/YMAX
	    if (s.le.YLEN) then
	        call plot(0.,s,3)
	        call plot(.6,s,2)
	        if (i.eq.7.or.i.ge.9) then
		    call nomber(.65,s,.1,float(i),0.,-1)
	    	    call setpat(1,2,1,2)
		    call dplot(.925,s,3)
		    call dplot(1.325,s,2)
		end if
	    end if
40	continue
c
	call strpos(.5,0.)
	call simbol(.6,YLEN,.225,'%E\\C\\M/@s^2%%',0.,999)
c
	r1 = 1.5
	s1 = decm(7)*YLEN/YMAX
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	call plot(-r1,0.,-3)
	call strpos(0.,.5)
	do 50 i=0,N
	    s = e(i)*YLEN/YMAX
	    if (s.le.YLEN) then
		call plot(0.,s,3)
		call plot(.6,s,2)
		if (i.gt.0) then
		    call nomber(.65,s,.1,float(i),0.,-1)
		else
		    call weight(1)
		    dr = .6/15
		    do 45 r = dr,.6,dr
			call plot(r-dr,s-2.5*dr,3)
			call plot(r,s,2)
45		    continue
		    call weight(10)
		end if
	    end if
50	continue
c
	call strpos(.5,0.)
	call simbol(.34,YLEN,.2,'E\\b\\i\\n/kT%%',0.,999)
c
	r0 = .8
c
	ss8 = e(8)*YLEN/YMAX
	call plot(r0,ss8,3)
	call plot(r0+.05,ss8,2)
	ss7 = e(7)*YLEN/YMAX
	call plot(r0+.05,ss7,2)
	call plot(r0,ss7,2)
c
	r0 = r0 + .05
	s0 = .5*(ss7+ss8)
	r1 = r1 - .05
c
	frac = .985
	r1 = r0 + frac*(r1-r0)
	s1 = s0 + frac*(s1-s0)
c
	call plot(r0,s0,3)
	call plot(r1,s1,2)
	frac = .1
	r0 = r1 - frac*(r1 - r0)
	s0 = s1 - frac*(s1 - s0)
	call arrow(r0,s0,r1,s1,1)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	call mcquit
c
	end
	
