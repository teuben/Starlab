
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

	subroutine errors(input,istart,nin,x,y,z,n,h,*)
        save
c
c	Add error bars from z:
c
	character*(*) input
	dimension x(1),y(1),z(1)
c
	common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
c
c	Decide on the direction of the bars:
c
	idir = 0
	do 162 i=istart,nin
	    if (input(i:i).gt.' ') then
		 if (input(i:i).eq.'x'
     &			.or.input(i:i).eq.'1') then
		    idir = 1
		else if (input(i:i).eq.'y'
     &			.or.input(i:i).eq.'2') then
		    idir = 2
		end if
		go to 10162
	    end if
162	continue
c
10162	if (idir.eq.0) go to 1001
	istart = i + 1
c
c	Draw the bars.
c
        call readiq(input(istart:nin),1,
     &              iside,idum,idum,idum,*1001)
c
c       Point size h is specified in screen units.
c       Allow an additional 10% to improve appearance.
c
	if (idir.eq.1) then
	    pt_size = 0.55*h/dinchx
	else
	    pt_size = 0.55*h/dinchy
	end if
c
	do 10562 i=1,n
	    if (abs(z(i)).gt.pt_size) then

		if (abs(iside).eq.1) then
c
c                   One-sided:
c
		    if (idir.eq.1) then
			call axerr(x(i),y(i),iside*z(i))
		    else
			call ayerr(x(i),y(i),iside*z(i))
		    end if
		else if (iside.eq.2) then
c
c                   Two-sided:
c
		    if (idir.eq.1) then
			call xerr(x(i),y(i),z(i))
		    else
			call yerr(x(i),y(i),z(i))
		    end if
		else
		    go to 1001
		end if

	    end if
10562	continue
c
	return
1001	return 1
c
	end
