
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
c	
c	A program plots a 3-dim picture
c	It calls the plt3d subroutine
c	a - the 2D data array, x direction increases faster.
c	m,n  -  dims of x and y directions
c	x and ywork - real working vector of length lwork >
c          2*min(m,n)
c	alt, azim - viewing altitute and azimuthal angles in
c			degrees
c	xlen,ylen - length of unprojected axis in device units
c	xoff,yoff - offset of plot origin in device units
c	zoff - data value corresponding to height of zero unit.
c	zfac - scale of z axis from data unit to unprojected
c		vertical device unit
c	ier - returns 0 for successful plot
c
 	integer*4 m,n,lwork,ier 
c
	real*4	a(500000),b(50000),xwork(1000),ywork(1000)
c
	read(1,*)m,n,etam,(a(k),k=1,(m+1)*(n+1)),(b(k),k=1,(m+1)*(n+1))
	m = m + 1
	n = n + 1
c
	amax=0.
	amin=1.e6
	bmax=0.
	bmin=1.e6
c
	do 20 k=1,m*n
		aa = max(-2., min(2., a(k)))
		amax = max(amax, aa)
		amin = min(amin, aa)
		a(k) = aa
c
		bb = max(-2., min(2., b(k)))
		bmax = max(bmax, bb)
		bmin = min(bmin, bb)
		b(k) = bb
10	    continue
20	continue
c
c	Rescale the data...
c
	abar = .5*(amin+amax)
	afac = 1./(amax-amin)
	bbar = .5*(bmin+bmax)
	bfac = 1./(bmax-bmin)
c
	do 110 k=1,m*n
	    a(k) = (a(k)-abar)*afac
	    b(k) = (b(k)-bbar)*bfac
110	continue
c
	call mcinit
c
	lwork=2*m
	xlen=6.
	ylen=6.
	xoff=5.
	yoff=5.
c
	zoff = 0.
	zfac = 6.
c
200	write(*,201)
201	format('Choice (1 or 2), Altitute, azimuthal viewing angles: '$)
	read(*,*,iostat=io)i,alt,azim
c
	if (io.eq.0) then
	    call clear
	    if (i.eq.1) then
		call plt3d(a,m,n,xwork,ywork,lwork,alt,azim,
     &		           xlen,xoff,ylen,yoff,zfac,zoff,ier)
	    else
		call plt3d(b,m,n,xwork,ywork,lwork,alt,azim,
     &		           xlen,xoff,ylen,yoff,zfac,zoff,ier)
	    end if
	    go to 200
	end if
c
99999	call mcquit
	end
