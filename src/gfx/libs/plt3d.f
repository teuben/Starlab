
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

        subroutine plt3d(a,m,n,alt,az,xlen,xoff,ylen,yoff,zfac,zoff,ier)
        save
c
c---	3d plot routine, adapted from the SAS-3 routine of the same name.
c
c---	calling sequence:
c	call plt3d (a,m,n,alt,az,xlen,xoff,ylen,yoff,zfac,zoff,ier)
c
c	a	- real data array, represents height of surface as
c		  function of location in plane
c	m,n	- dimensions of data array a
c	alt, az	- altitude, azimuth viewing angles in degrees
c	xlen,ylen - length of unprojected x,y axes in inches
c	xoff,yoff - offset of origin in inches
c	zfac    - scaling of z-axis from data units to unprojected inches
c	zoff	- offset of z-origin in user units
c	ier	- returns 0 for no error
c
c	Altitude is measured from the x-y plane toward the positive z-axis.
c	Azimuth is measured around the z-axis (as seen from above), 
c	counterclockwise from the positive x-axis.
c
c	The origin (xoff, yoff), in screen coordinates, is the image of the
c	point at the CENTER of the x-y grid and at a height zoff above it
c	(in the units used for a).
c
        dimension a(m,n)
c
c	Note that the original user-provided work arrays are now internal
c	and static.
c
	parameter (LDIM = 5000)
	dimension x(LDIM),y(LDIM)
c
        common /plt3b/ a1,a2,a3,b1,b2,b3,b4,zoff1
c
	parameter (PI = 3.141592654, PIDEG = PI/180.)
c
c	Coordinate to plot conversions:
c
c		(1,1) <--> (xmin,ymin)
c		(1,n) <--> (xmin,ymax)
c		(m,1) <--> (xmax,ymin)
c		(m,n) <--> (xmax,ymax)
c
        xplot(i,j) = a1*i + a2*j + a3
        yplot(i,j) = b1*i + b2*j + b3*(a(i,j)-zoff1) + b4
c
        lmax = 2*min(m,n)
        if (LDIM.lt.lmax) then
	    write(0,*)'PLT3D: Array too big!'
	    ier = 2
	    return
	end if
c
        taz = az*PIDEG
        talt = alt*PIDEG
c
        saz = sin(taz)
        caz = cos(taz)
        sal = sin(talt)
        cal = cos(talt)
c
        xsc = xlen/float(m-1)
        ysc = ylen/float(n-1)
c
c	The viewing plane is such that the "horizontal" vector is
c
c	    x' = (-saz, caz, 0),
c
c	and the "vertical" vector is
c
c	    y' = (-sal*caz, -sal*saz, cal).
c
        a1 = -saz*xsc
        a2 = caz*ysc
        a3 = xoff - 0.5*(a1*float(m+1)+a2*float(n+1))
        b1 = -caz*sal*xsc
        b2 = -saz*sal*ysc
        b3 = zfac*cal
        b4 = yoff - 0.5*(b1*float(m+1)+b2*float(n+1))
	zoff1 = zoff
c
c	Find which quadrant we are viewing from.
c
        iaz = 1
        if (caz.le.0.0) iaz = iaz + 1
        if (saz.le.0.0) iaz = 5 - iaz
c
c	Choose the direction of traversal of the plot based on view angle.
c
	if (iaz.eq.2.or.iaz.eq.3) then
            ifirst = 1
	else
            ifirst = m
	end if
	ilast = m + 1 - ifirst
	istep = sign(1,ilast-ifirst)
c
	if (iaz.ge.3) then
            jfirst = 1
	else
            jfirst = n
	end if
	jlast = n + 1 - jfirst
	jstep = sign(1,jlast-jfirst)
c
c	Make sure the array sent to nxtvu is traversed from left to right.
c
	if (iaz.eq.1.or.iaz.eq.3) then
	    lli = -1
	else
	    lli = 1
	end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        ic = 0
        ibeg = ifirst+istep
c
c	Traverse the plot region in a zig-zag motion, transverse to the
c	viewing direction (left to right), and from front to back.
c	The plotting routine "nxtvu" takes care of hidden line removal.
c
c	First, deal with the "front" portion of the plot (i.e. nearest
c	the viewer).
c
 70     lnth = min(2*iabs(ibeg-ifirst)+1,lmax)
	if (lli.eq.-1) then
	    ll = lnth + 1
	else
	    ll = 0
	end if
c
        i = ibeg
        j = jfirst
        ll = ll+lli
        x(ll) = xplot(i,j)
        y(ll) = yplot(i,j)
c
 80     i = i-istep
        ll = ll+lli
        x(ll) = xplot(i,j)
        y(ll) = yplot(i,j)
        if (j.ne.jlast) then
            j = j+jstep
            ll = ll+lli
            x(ll) = xplot(i,j)
            y(ll) = yplot(i,j)
            if (i.ne.ifirst) go to 80
	end if
c
        call nxtvu(ic,x,y,lnth,ier)
        if (ier.ne.0) return
        ic = 1
        if (ibeg.ne.ilast) then
            ibeg = ibeg+istep
            go to 70
	end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	Now do the "rear" half of the figure in the same manner.
c
        jbeg = jfirst
c
 100    lnth = min(2*iabs(jbeg-jlast)+1,lmax)
	if (lli.eq.-1) then
	    ll = lnth + 1
	else
	    ll = 0
	end if
c
        i = ilast
        j = jbeg
        ll = ll+lli
        x(ll) = xplot(i,j)
        y(ll) = yplot(i,j)
c
 110    j = j+jstep
        ll = ll+lli
        x(ll) = xplot(i,j)
        y(ll) = yplot(i,j)
        if (i.ne.ifirst) then
            i = i-istep
            ll = ll+lli
            x(ll) = xplot(i,j)
            y(ll) = yplot(i,j)
            if (j.ne.jlast) go to 110
	end if
c
        call nxtvu(1,x,y,lnth,ier)
        if (ier.ne.0) return
        jbeg = jbeg+jstep
        if (jbeg.eq.jlast) return
	go to 100
c
        end


	function xplot3d(x,y,z)
        save
c
c	Externally accessible version of the statement function in plt3d.
c	This will work only AFTER plt3d has been used, as that routine
c	sets up the common array.
c
c	The coordinates x and y are expected to range from 1 to m and
c	1 to n, respectively.
c
c	Typical usage:  xx = xplot3d(float(i),float(j),array(i,j))
c
c	Note that z is not used.
c
        common /plt3b/ a1,a2,a3,b1,b2,b3,b4,zoff1
c
        xplot3d = a1*x + a2*y + a3
c
	end


	function yplot3d(x,y,z)
        save
c
c	Externally accessible version of the statement function in plt3d.
c	This will work only AFTER plt3d has been used, as that routine
c	sets up the common array.
c
c	The coordinates x and y are expected to range from 1 to m and
c	1 to n, respectively.
c
        common /plt3b/ a1,a2,a3,b1,b2,b3,b4,zoff1
c
c	Typical usage:  yy = yplot3d(float(i),float(j),array(i,j))
c
        yplot3d = b1*x + b2*y + b3*(z-zoff1) + b4
c
	end
