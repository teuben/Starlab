
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

	subroutine sundummies
c
c	Dummy entry points to avoid having the Sun graphics libraries
c	linked in to a code that doesn't want them.
c
        character*(*) string
c
	entry plinit(string,aspect,icolor,igp,ncolor,ierr)
c
c       Print a message and return an error status.
c
        write(6,*)'SunCore not available'
        ierr = 1
        return
c
	entry awtbuttongetloc2
	entry inqtextextent2
	entry lineabs2
	entry moveabs2
	entry plframe
	entry plstop
	entry polygonabs2
	entry setcharsize
	entry setfillindex
	entry setlineindex
	entry setlinewidth
	entry setpolyintrstyle
	entry settextindex
	entry text
        entry corebackg
        entry readmap
        entry sunwait
c
	end

