
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

	subroutine smooth(t,x,xsmooth,n,dts,iopt,iw)
        save
c
c	Smooth the input data, using one of a variety of methods.
c
c       Options:
c
c           iopt = 0:  median
c                  1:  arithmetic mean
c                  2:  harmonic mean
c                  3:  geometric mean
c                  11: clipped arithmetic mean
c                  12: clipped harmonic mean
c                  13: clipped geometric mean
c
c           iwindow = 0:  flat window
c                     1:  triangular window
c      
	parameter (NMAX = 10000, FCLIP = 0.1)
c
	dimension t(1),x(1),xsmooth(1),work(NMAX),tt(NMAX)
c
	external comp
c
c	Smooth the input data over a window of width dts.
c	Brute-force approach recalculates the sum at each data point!
c
c	Window function:
c
	window(xx) = max(0., 1. - iwindow*abs(xx/dt2))
c
	if (dts.le.0.) return
        iwindow = iw
        if (iwindow.ne.0) iwindow = 1
	dt2 = .5*dts
c
	if (iopt.gt.0.and.iopt.le.10.and.iwindow.eq.0) then
c
c	    Flat window, straight average only:
c
	    sum = 0.
	    ns = 0
	    j0 = 1
	    j1 = 0
c
	    do 190 i=1,n
c
120             continue
c
c	        if (j1.lt.n
c     &                  .and.t(j1+1).le.t(i)+min(dt2,t(i)-t(1))) then
c
	        if (j1.lt.n.and.t(j1+1).le.t(i)+dt2) then
		    j1 = j1 + 1
		    sum = sum + f(x(j1),iopt)
		    ns = ns + 1
		    go to 120
	        end if
c
130	        continue
c               if (t(j0).lt.t(i)-min(dt2,t(n)-t(i))) then
                if (t(j0).lt.t(i)-dt2) then
		    sum = sum - f(x(j0),iopt)
		    ns = ns - 1
		    j0 = j0 + 1
		    go to 130
	        end if
c
	        xsmooth(i) = finv(sum/max(ns,1),iopt)
c
190	    continue
	else
c
c	    All other options:
c
	    j0 = 1
	    j1 = 0
c
c           Locate the range containing the desired interval.
c
	    do 290 i=1,n
c
220	        if (j1.lt.n
     &                  .and.t(j1+1).le.t(i)+min(dt2,t(i)-t(1))) then
		    j1 = j1 + 1
		    go to 220
	        end if
c
230	        if (t(j0).lt.t(i)-min(dt2,t(n)-t(i))) then
		    j0 = j0 + 1
		    go to 230
	        end if
c
c               Range is j0 to j1.
c
                do 240 j=j0,j1
	    	    work(j-j0+1) = x(j)
                    tt(j-j0+1) = t(j)
240             continue
                nw = j1-j0+1
c
                if (iopt.eq.0.or.iopt.gt.10) then
c
c                   Sort the data.
c
	    	    if (nw.gt.1) call sort2(nw,work,tt)
c
                    if (iopt.eq.0) then
c
c                       Interpolate to the median.
c
                        xj = .5*(nw+1)
                        jj = xj
                        xsmooth(i) = work(jj)
     &                               + (xj-jj)*(work(jj+1) - work(jj))
                    else
c
c                       Clip the data.
c
                        xj1 = FCLIP*(nw+1)
                        jj1 = max(1,nint(xj1))
                        xj2 = (1.-FCLIP)*(nw+1)
                        jj2 = min(nw,nint(xj2))
c
c                       New range is jj1 to jj2.
c
                        do 250 j=jj1,jj2
                            work(j-jj1+1) = work(j)
                            tt(j-jj1+1) = tt(j)
250                     continue
                        nw = jj2-jj1+1
                    end if
                end if
c
		if (iopt.gt.0) then
c
c                   Perform averaging on the work data.
c
	            sum = 0.
	            wsum = 0.
c
	            do 260 j=1,nw
		        weight = window(tt(j)-t(i))
		        wsum = wsum + weight
		        sum = sum + weight*f(work(j),iopt)
260	            continue
c
	            if (wsum.gt.0.) then
		        xsmooth(i) = finv(sum/wsum,iopt)
	            else
		        xsmooth(i) = x(i)
	            end if
		else
		end if
c
290	    continue
	end if
c
	end



	function f(x,iopt)
	save
c
	if (iopt.eq.1.or.iopt.eq.11) then
	    f = x
	else if (iopt.eq.2.or.iopt.eq.12) then
	    if (x.ne.0.) then
		f = 1./x
	    else
		f = 0.
	    end if
	else if (iopt.eq.3.or.iopt.eq.13) then
	    if (x.gt.1.e-20) then
	 	f = log(x)
	    else
		f = -46.05
	    end if
	end if
c
	end



	function finv(x,iopt)
	save
c
c	Compute the inverse function to f above.
c
	if (iopt.eq.1.or.iopt.eq.11) then
	    finv = x
	else if (iopt.eq.2.or.iopt.eq.12) then
	    if (x.ne.0.) then
		finv = 1./x
	    else
		finv = 0.
	    end if
	else if (iopt.eq.3.or.iopt.eq.13) then
	    if (x.gt.-46.05) then
	 	finv = exp(x)
	    else
		finv = 1.e-20
	    end if
	end if
c
	end



	integer*2 function comp(x,y)
	save
c
c	Comparison function for use with qsort.
c
	if (x.gt.y) then
	    comp = 1
	else if (x.lt.y) then
	    comp = -1
	else
	    comp = 0
	end if
c
	end
