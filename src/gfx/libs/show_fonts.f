
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

	program fonts
c
        parameter (BOTTOM = 0.25, TOP = 6.9, NROWS = 24)
c
	data r0/.5/dr/2.25/
c
        ds = (TOP - BOTTOM) / (NROWS - 1)
        s0 = TOP + ds
	ds2=.5*ds
	dr1=2.*ds
c
        call mcinita(' ')
c        call plot(1., 1., 3)
c        call plot(10., 10., 2)
c        call mcdxidle
c        call mcdxnobuffer
c        write(6,*)'After mcinit'
	call nobounds
        call weight(7)
	call simbol(r0,s0+.75*ds,ds2,
     &		    'Here is the present complete SIMBOL font set'//
     &		    ' ("@Aa"\\ _\\ @a, "@Ra"\\ _\\ %a):%%',
     &		    0.,999)
c
	i = 31
	r = r0 - dr
	do 100 jcol=1,4
	    r = r + dr
	    s = s0
	    do 75 irow=1,24
		i = i + 1
		s = s - ds
		call nombr(r,s,ds2,float(i),0.,-1)
		call simbol(r+dr1,s,ds2,char(i),0.,1)
		call simbol(r+dr1+ds,s,ds2,'@'//char(i),0.,2)
		call simbol(r+dr1+2.*ds,s,ds2,'%'//char(i),0.,2)
75	    continue
100	continue
c
	call mcquit1
	end
