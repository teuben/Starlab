
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


	subroutine boxf(x1,x2,y1,y2,icolor)
        save
c
c	Like box, but fills the rectangle with a specified color.
c
	dimension x(4),y(4)
c
	x(1)=x1
	y(1)=y1
	x(2)=x2
	y(2)=y1
	x(3)=x2
	y(3)=y2
	x(4)=x1
	y(4)=y2
	call setfill(icolor)
	call polyfill(x,y,4)
	call unsetfill
	call box(x1,x2,y1,y2)
c
	return
	end
