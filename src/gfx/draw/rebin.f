
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

	subroutine rebin(string,array,nmax,narr,istat,iout,temp,list)
        save
c
c	Rebin an array, carrying the others along (if same length).
c
	character*(*) string
	dimension array(nmax,3),narr(3),temp(nmax),list(nmax)
	character*40 token(2)
c
	call gettokens(string,token,nt)
	if (nt.ne.2) go to 1001
c
	call sdecode(token(1),1,iarr,*1001)
	if (narr(iarr).le.0) go to 1001
c
	call readrtoken(token(2),dx,0.)
	if (dx.eq.0.) go to 1001
c
	if (array(narr(iarr),iarr).lt.array(1,iarr)) then
	    dx = -abs(dx)
	else
	    dx = abs(dx)
	end if
c
c	Linearly rebin the specified array.  Assume the elements are
c	ordered, and make a list of elements to save.
c
	nl = 1
	list(1) = 1
	xfirst = array(1,iarr)
	x = xfirst
	xnext = xfirst + dx
	do 100 i=2,narr(iarr)
	    xprev = x
	    x = array(i,iarr)
	    if ((x-xfirst)*(x-xprev).lt.0.) then
		if (iout.eq.1) write(6,*)'Array not ordered'
		go to 1001
	    end if
	    if ((x-xnext)*(x-xfirst).ge.0.) then
		nl = nl + 1
		if (nl.gt.nmax) go to 1001
		list(nl) = i
90		xnext = xnext + dx
		if ((x-xnext)*(x-xfirst).ge.0.) go to 90
	    end if
100	continue
c
	if (list(nl).ne.narr(iarr)) then
	    nl = nl + 1
	    list(nl) = narr(iarr)
	end if
c
c	Rebin the array(s).
c
	nold = narr(iarr)
	do 200 k=1,3
	    if (narr(k).eq.nold) then
		do 150 i=1,nl
150		temp(i) = array(list(i),k)
		do 160 i=1,nl
160		array(i,k) = temp(i)
		narr(k) = nl
	    end if
200	continue
c
	if (iout.eq.1) write(6,*)'New number of points = ',nl
	return
c
1001	istat = 3
	end
