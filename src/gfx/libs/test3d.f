
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
	program test3d
c	
c	A program to plot a 3-D picture, using the plt3d subroutine.
c
c	a - the 2D data arrary, x direction increases faster.
c	m,n  -  dimensionss of x and y directions
c	alt, azim - viewing altitute and azimuthal angles in
c			degrees
c	xlen,ylen - length of unprojected axis in device units
c	xoff,yoff - offset of plot origin in device units
c	zoff - data value corresponding to height of zero unit.
c	zfac - scale of z axis from data unit to unprojected
c		vertical device unit
c	ier - returns 0 for successful plot
c
 	integer*4 m,n,ier 
	parameter (m = 50, n = 75)
	real*4	a(m,n),alt,azim,xlen,xoff,ylen,yoff,zfac,zoff
c
	amax=0.
	amin=1.e30
c
	do 5 j=1,n
	    y = float(j)/n
	    do 1 i=1,m
		x = float(i)/m
c
		a(i,j) = x*y*sin(12.56*y)
c
		amax = max(amax,a(i,j))
		amin = min(amin,a(i,j))
1	    continue
5	continue
c
	call mcinit
c
	xlen=6.
	ylen=6.
	xoff=5.
	yoff=5.
	zfac=7./(amax-amin)
	zoff=.5*(amin+amax)
c
2	write(*,'(''Altitute, azimuthal viewing angles: ''$)')
	read(*,*,err=99999,end=99999)alt,azim
	call clear
c
	call plt3d(a,m,n,alt,azim,xlen,xoff,ylen,yoff,zfac,zoff,ier)
	if (ier.ne.0) write(6,*)'Return status = ',ier
c
	call simbol(xplot3d(float(m+1),-1.,0.),
     &		    yplot3d(float(m+1),-1.,0.),
     &		   .225,'X',0.,1)
	call simbol(xplot3d(-1.,float(n+1),0.),
     &		    yplot3d(-1.,float(n+1),0.),
     &		   .225,'Y',0.,1)
c
	go to 2
c
99999	call mcquit
	end
