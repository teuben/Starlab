
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
c       X-window calls used by mcdraw:
c       -----------------------------
c
c       Initialization (ierr = 0 on success):
c
        subroutine mcdxinit(ierr)
        save
c
c       Set line width:
c
        entry mcdxlinew(iw)
c
c       Set line color:
c
        entry mcdxcolor(ic)
c
c       Set background color:
c
        entry mcdxbackg(ic)
c
c       Move "cursor" to (r,s):
c
        entry mcdxmove(r,s)
c
c       Draw from current cursor location to (r,s):
c
        entry mcdxdraw(r,s)
c
c       Plot a point (pixel) at (r,s):
c
        entry mcdxpoint(r,s)
c
c       Get graphics input -- mouse clicked at (r,s):
c
        entry mcdxgin(r,s)
c
c       Draw a filled polygon (points {(r(i),s(i)), i=1,..,n}, fill color if):
c
        entry mcdxpolyf(r,s,n,if)
c
c       Erase a polygon (points {(r(i),s(i)), i=1,..,n}):
c
        entry mcdxpolyc(r,s,n)
c
c       Draw a text string at (r,s), height = h, angle = a:
c
        entry mcdxtext(r,s,h,a,string)
c
c       Clear the display:
c
        entry mcdxclear
c
c       Reset (reinitialize) the display:
c
        entry mcdxreset
c
c       Quit the X-display:
c
        entry mcdxquit
c
c       X-idle modes:
c
        entry mcdxidle
        entry mcdxread_line
c
        entry mcdxnopen
c
c       Return number of open X windows.
c
        entry mcdxcurrwin
c
c       Return ID of current X window.
c
        entry mcdxsetwin
c
c       Set current X window.
c
        entry mcdxkillwin
c
c       Kill an X window.
c
        end

