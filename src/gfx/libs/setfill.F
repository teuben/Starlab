
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

	subroutine setfill(icolor)
        save
c
c	Set the SunCore/X/PostScript fill index.
c
        character*80 device
        common /plot device/ device,aspect,idev
	common /findex/ index
c
#ifdef HAS_SUNCORE
c	include "/usr/include/f77/usercore77.h"
	include "/usr/old/usercore77.h"
#endif
	data index/-1/
c
	if (idev.ne.15.and.idev.ne.16.and.idev.ne.17) return
c
	index = icolor
c
        if (idev.eq.15) then
            if(icolor.ge.0)then
                call setfillindex(icolor)
                call setpolyintrstyle(PLAIN)
            end if
        end if
c
	return
c
	entry unsetfill
	index=-1
	return
c
        entry getfill(j)
c
c       Return current fill index.
c
        j = index
c
	end
