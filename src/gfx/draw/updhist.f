
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

	subroutine updhist(line,nl,nzoom,nbox)
        save
c
c	Update the history list.
c
	character*(*) line
c
	parameter (NHMAX = 500)
c
        character*200 history(NHMAX)
        common/histchars/history
        common/histnums/lhist(NHMAX),nhist,ishist(NHMAX),
     &		        izhist(NHMAX),ibhist(NHMAX)
c
        nhist=nhist+1
        jh=nhist
10015   if (jh.gt.NHMAX) then
            jh=jh-NHMAX
            go to 10015
        end if
c
        lhist(jh)=min(200,nl)
        history(jh)=' '
        history(jh)=line(1:lhist(jh))
c
	call getnsave(nsave)
	ishist(jh) = nsave
	izhist(jh) = nzoom
	ibhist(jh) = nbox
c
	return
	end
