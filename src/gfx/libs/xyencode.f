
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
c
	subroutine xyencode(r,s,vec)
        save
c
c	Encode a Tektronix-type vector.
c
	character*1 vec(*)
	common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
	common /plot origin/ ro,so
c
	i=nx0+xfac*(r+ro)
	j=ny0+yfac*(s+so)
	if(nxpix.gt.1023)then
	    i4=i/4
	    ii=i-4*i4
	    j4=j/4
	    jj=j-4*j4
	    vec(2)=char(96+ii+4*jj)
	    i=i4
	    j=j4
	else
	    vec(2)='`'
	end if
	j32=j/32
	vec(1)=char(32+j32)
	vec(3)=char(96+(j-32*j32))
	i32=i/32
	vec(4)=char(32+i32)
	vec(5)=char(64+(i-32*i32))

	end
