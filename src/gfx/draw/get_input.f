
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

	subroutine getinput
        save
c
c	General means of getting input from stdin, the graphics cursor,
c	or the repeat stack, depending on mode.  This is fairly wasteful
c	of space, but...
c
	common /replay/ ireplay
	common /prompt/ iprompt
	common /x input/ interact
	logical gfxin,first_xuse
c
        character*80 device
        common /plot device/ device,aspect,idev
c
	common/draw params/roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &			   idevset,jbox,iorig
	parameter (NSAVEMAX = 500)
c
	common /instack1/ nsave,isave,rsave(NSAVEMAX),ssave(NSAVEMAX)
	character*200 strsave(NSAVEMAX)
	common /instack2/ strsave
c
	save /instack1/,/instack2/
c
	character*(*) input,string
	data first_xuse/.true./
c
	entry getnsave(i)
	i = nsave
	return
c
	entry setisave(i)
	isave = i
	return
c
	entry getgfx(r,s)
c
c	Graphics input.
c
	xl = xlen
	if (xl.le.0.) xl = 1.
	yl = ylen
	if (yl.le.0.) yl = 1.
c
	if (ireplay.eq.0) then
	    if (gfxin()) then
c
		if (idev.eq.17.and.first_xuse.and.iprompt.eq.1)
     &          write(6,'(a)')
     &             'Use right-hand mouse button to indicate position.'
c
		call graphin(r,s)
c
		if (idev.eq.17) first_xuse = .false.
	    else
		if (iprompt.eq.0) return
		call devoff
		write(6,'(''No graphics input.  Enter r, s: ''$)')
		read(5,*,end=99,err=99)r,s
	    end if
c
	    if (nsave.lt.NSAVEMAX) then
		nsave = nsave + 1
		rsave(nsave) = r/xl
		ssave(nsave) = s/yl
	    end if
	else
	    if (isave.lt.nsave) then
		isave = isave + 1
		r = xl*rsave(isave)
		s = yl*ssave(isave)
	    end if
	end if
c
99	return
c
	entry getstring(input,istart,nin,string)
c
c	General keyboard input, with optional prompt.
c
	if (ireplay.eq.0) then
c
	    call devoff
	    if (istart.le.nin) then
		write(6,'(a,'' '',$)')input(istart:nin)
	    else
		write(6,'(''Input string: ''$)')
	    end if
	    string = ' '
	    if (interact.eq.1.and.num_win(17).gt.0) then
c
c               Get input via X, keeping screen up to date.
c
		call myflush(6)
		string = '\0'
c
		call win_read_line(string)
c
	    else
		read(5,'(a)',err=999,end=999)string
	    end if
c
	    if (nsave.lt.NSAVEMAX) then
		nsave = nsave + 1
		strsave(nsave) = string(1:min(200,len(string)))
	    end if
	else
	    if (isave.lt.nsave) then
		isave = isave + 1
		string = strsave(isave)
	    end if
	end if
c
999	return
c
	end
