
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

	subroutine genfit(input,istart,nin,x,y,z,n,nz,iprompt,*)
        save
c
c	Perform a general least-squares fit to y(x).
c
	character*(*) input
c
	parameter (MMAX = 21)
	dimension x(1),y(1),z(1),a(MMAX)
c
	external poly,trig
c
c	Equal weights for all points if nz <= 0.
c	Otherwise, specify the weights in the z array...
c
	nw = nz
	if (nz.le.0) nw = 1
c
c	Determine the array ranges.  Note that the arrays MUST be sorted
c	in x for this to work.
c
	call sort3(n,x,y,z)
c
10      if (istart.gt.nin) then
            i1 = 1
            nn = n
	    ma = 2
        else
	    ma = 1
	    x1 = -1.e20
	    x2 = 1.e20
            call readrq(input(istart:nin),3,
     &                  xma,x1,x2,dum,*1001)
            ma = nint(xma)
	    if (ma.le.0) ma = 1
c
	    ma = ma + 1
	    if (ma.gt.MMAX) then
		if (iprompt.eq.1) then
		    call devoff
		    write(6,*)' Must specify order <= ',MMAX-1
		    go to 1001
		end if
	    end if
c
            i1=0
            do 150 i=1,n
                if (x(i).gt.x1.and.i1.eq.0) i1=i
                if (x(i).gt.x2) go to 160
150         continue
            i = n+1
160         i1 = max(1,i1)
            nn = i-i1
            if (nn.le.1) then
                write(6,*)'No points in specified range.'
                go to 1001
            end if
        end if
c
c	Note: Input ma is the order of the polynomial fit (<=20).
c	      Internal ma is one greater than this number.
c
c	Rescale the x and y data to lie in [-1,1], to reduce the effects of
c	finite-precision arithmetic.  This appears to do a better job than 
c	scaling to [0,1].
c
	call minmax(x,n,xmin,xmax)
	xscale = 2./(xmax - xmin)
	xref = .5*(xmin+xmax)
c
	call minmax(y,n,ymin,ymax)
	yscale = 2./(ymax - ymin)
	yref = .5*(ymin+ymax)
c
	do 200 i=1,n
	    x(i) = (x(i) - xref)*xscale
	    y(i) = (y(i) - yref)*yscale
200	continue
c
	xscale = 1./xscale
	yscale = 1./yscale
c
	x1 = x(i1)
	x2 = x(i1+nn-1)
c
c	Use an intermediate C routine to allocate workspace space.
c	(Pointers are NOT standard FORTRAN yet!)
c
	if (input(2:2).eq.'t') then
	    call cfit(x(i1),y(i1),z(i1),nn,a,ma,trig,nw,chisq,iret)
	else
	    call cfit(x(i1),y(i1),z(i1),nn,a,ma,poly,nw,chisq,iret)
	end if
c
c	Print out some basic data on the fit obtained.
c
        if (iprompt.eq.1) then
            call devoff
	    write(6,250)(a(i),i=1,ma)
250         format(' Coefficients: ',1p,5e11.3,:/(15x,5e11.3))
	    write(6,*)' Chisq = ',chisq
	end if
c
        if (input(3:3).eq.'p') then
c
c	    Plot the fit.
c
	    call uplotin(x1*xscale+xref,
     &			 afunc(x1,a,ma,input(2:2))*yscale+yref,3)
c
	    do 300 i=1,n
300	    if (x(i).gt.x1.and.x(i).le.x2)
     &		call uplotin(x(i)*xscale+xref,
     &			     afunc(x(i),a,ma,input(2:2))*yscale+yref,2)
        end if
c
c	Unscale x and y before returning.
c
	do 350 i=1,n
	    x(i) = xscale*x(i) + xref
	    y(i) = yscale*y(i) + yref
350	continue
c
	return
1001	return 1
c
	end


	function afunc(x,a,ma,which)
	save
c
	parameter (MMAX = 21)
	dimension a(ma),f(MMAX)
	character*1 which
c
	if (which.eq.'t') then
	    call trig(x,f,ma)
	else
	    call poly(x,f,ma)
	end if
c
	sum = 0.
	do 10 i=1,ma
10	sum = sum + a(i)*f(i)
c
	afunc = sum
c
	return
	end
