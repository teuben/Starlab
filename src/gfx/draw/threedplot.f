
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

	subroutine ez3dplot(a,na,iprompt,*)
        save
c	
c	Plot a 3-D picture, using the array a.
c
	dimension a(1)
	character string*80,c*1
        common /mcd color/ icolor
c
	save m,n
	data m/0/n/0/
c
	call devoff
	if (iprompt.eq.1) write(6,'(''*** Three-dimensional plot'')')
c
        call get_dimensions(m,n,na,iprompt,*99999)
c
	call minmax(a,m*n,amin,amax)
	if (iprompt.eq.1) write(6,'(a,f,a,f)')
     &          'Array minimum = ',amin,', maximum = ',amax
	if (amin.ge.amax) return
c
	xlen=6.
	ylen=6.*aspectratio()
	xoff=5.
	yoff=5.25*aspectratio()
	zfac=6./(amax-amin)
	zoff=.5*(amin+amax)
c
10	if (iprompt.eq.1) then
	    string = ' '
	    call getstring('Altitude, azimuthal viewing angles, '//
     $	                   '[color] (h = help): ',1,56,string)
	else
	    call getstring(' ',1,1,string)
	end if
c
	do 20 i=1,len(string)
	    if (string(i:i).gt.' ') then
		c = string(i:i)
		if (c.ge.'A'.and.c.le.'Z') c = char(ichar(c)+32)
	 	go to 30
	    end if
20	continue
	go to 10
c
30	if (c.eq.'q'.or.c.eq.'x'.or.c.eq.'e') then
	    go to 9999
	else if (c.eq.'h') then
	    write(6,*)
	    write(6,*)'e or x = back to mcdraw prompt'
	    write(6,*)'q = quit mcdraw'
	    write(6,*)'altitude = angle above the x-y plane'
	    write(6,*)'azimuth  = angle around the z-axis, '//
     &		      'measured from the x-axis toward the y-axis'
	    write(6,*)'color = "underside" color, if specified'
	    write(6,*)
	    go to 10
	end if
c
	jcolor = 0
        call readrq(string,3,alt,azim,xjcolor,xdum,*99)
99	jcolor = nint(xjcolor)
c
	call devon
	call clear
c
	if (jcolor.gt.0) then
c
c	    Draw the "underside" first:
c
	    call setvu(2)
	    call color(jcolor)
	    call plt3d(a,m,n,alt,azim,xlen,xoff,ylen,yoff,zfac,zoff,ier)
	    call color(icolor)
	    call setvu(1)
	end if
c
	call plt3d(a,m,n,alt,azim,xlen,xoff,ylen,yoff,zfac,zoff,ier)
c
	ro = xplot3d(1.,1.,amin)
	so = yplot3d(1.,1.,amin)
	rx = xplot3d(m+1.,1.,amin)
	sx = yplot3d(m+1.,1.,amin)
	ry = xplot3d(1.,n+1.,amin)
	sy = yplot3d(1.,n+1.,amin)
	rz = xplot3d(1.,1.,amin+1.1*(amax-amin))
	sz = yplot3d(1.,1.,amin+1.1*(amax-amin))
c
c	Which axes should be visible?
c
c	call plot(ro,so,3)
c	call plot(rx,sx,2)
c	call plot(ro,so,3)
c	call plot(ry,sy,2)
c	call plot(ro,so,3)
c	call plot(rz,sz,2)
c
	call devoff
	if (ier.ne.0) write(6,'(a,i)')'Return status = ',ier
c
	go to 10
c
9999	return
99999	return 1
	end
