
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
	program example
c
	dimension x(100),y(100),z(100)
c
	do 10 i=1,100
	    x(i) = 0.1*i
	    y(i) = x(i)*sin(.1*x(i)**2)
	    z(i) = 0.03*x(i) - 1.0
10	continue
c
	call mcinit
c
c	Move and re-origin:
c	
	call plot(2.0,1.5,-3)
c
	call minmax(x,100,xmin,xmax)
	call minmax(y,100,ymin,ymax)
	xlen = 7.0
	ylen = 4.0
	modex = 1
	modey = 2
c
	call eframe(xmin,xmax,xlen,modex,'x-axis',
     &		    ymin,ymax,ylen,modey,'y-axis')
c     &		    ymin,ymax,ylen,modey,' y-axis')
c
	call mline(x,y,100,0,0,0.0)
c
	call simbol(4.0,5.0,0.2,'Hi, Paul!',0.0,9)
	call strpos(0.75,0.5)
	call simbol(4.0,4.5,0.2,'Hi, Paul!%%',0.0,999)
	call clrstr
	call simbol(1.0,4.0,0.2,'@H@i, @P%a%u%l%!%%',45.0,99)
c
	call newplot
c
	modey = -1
	zmin = 0.1
	zmax = 100.0
	call sethts(.2,.3)
	call eframe(xmin,xmax,xlen,modex,'x-axis',
     &		    zmin,zmax,ylen,modey,'A much longer y-axis label')
c
	call setpat(1,2,1,3)
	call dline(x,z,100,10,5,0.1)
	call weight(10)
	call mline(x,y,100,-5,-(96+ichar('p')),0.5)
c
	call weight(1)
	call simbol(0.5*xlen,ylen+0.25,0.2,'Happy hacking...%%',0.0,-999)
c
	call mcquit
c
	end
