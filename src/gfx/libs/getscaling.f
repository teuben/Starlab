
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


 	subroutine get scaling(x1,x2,dx,xl,y1,y2,dy,yl)
        save
 	common/scales/ xmin,xmax,dxinch,ymin,ymax,dyinch,rlen,slen
c
 	entry get scales(x1,x2,dx,xl,y1,y2,dy,yl)
c
 	x1=xmin
 	x2=xmax
 	dx=dxinch
 	xl=rlen
 	y1=ymin
 	y2=ymax
 	dy=dyinch
 	yl=slen
c
	return
c
	entry set scaling(xref,fx,yref,fy)
	entry set scales(xref,fx,yref,fy)
c
c	Set essential scaling without a call to eframe, so uplot,
c	mline, etc. can be used.
c
	xmin=xref
	dinchx=fx
	ymin=yref
	dinchy=fy
c
 	return
 	end
