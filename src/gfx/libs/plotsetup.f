
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

        block data plot set init
        save
c
        common /dev status/ idevon,idevpen,idevwt
        common /plot sizes/ xsize,ysize
	common /plot invert/ inv,ipensto
        common /numbr on/ inumbr
c
        common /dev details/ itek,ivers
        common /hp plot/ ivdef
	common /sunscreen/ isun
	common /findex/ index
c
        data idevon/-1/idevpen/-1/idevwt/-1/xsize/-1./inv/-1/
        data ipensto/-1/inumbr/-1/itek/-1/ivers/-1/ivdef/5/isun/-1/
c
        end


	subroutine plot setup
        save
c
c	Initialize ALL package variables, so mcinit always (re)starts in
c	a standard state.
c
        common /dev status/ idevon,idevpen,idevwt
        common /plot sizes/ xsize,ysize
	common /plot invert/ inv,ipensto
        common /numbr on/ inumbr
c
        common /dev details/ itek,ivers
        common /hp plot/ ivdef
	common /sunscreen/ isun
	common /findex/ index
c
        logical set
c
c	Set/reset plot parameters...
c
	if (idevon.lt.0) idevon = 0
	if (idevpen.lt.0) idevpen = 1
	if (idevwt.lt.0) idevwt = 1
	if (xsize.lt.0.) xsize = 10.
	if (inv.lt.0) inv = 0
	if (ipensto.lt.0) ipensto = 1
	if (inumbr.lt.0) inumbr = 0
c
	call setorigin(0.,0.)
	call setlhe(0.)
	call setbot(0.)
	call bounded
c
	if (itek.lt.0) itek = 0
	if (ivers.lt.0) ivers = 0
	ivdef = 5
	if (isun.lt.0) isun = 0
	index = -1
c
c	...and frame parameters...
c
c	BE CAREFUL setting up default "frame" options, as these
c       may reset the graphics in an unwanted manner.
c
        call gethset(set)
        if (.not.set) call sethts(.25,.2)
c
        call getmset(set)
        if (.not.set) call setmod(0,1,0,0,0)
c
	call getwset(set)
        if (.not.set) call setwts(0,0,0,0)
c
	call getpset(set)
        if (.not.set) call setpens(0,0,0)
c
	call clrstars
c
	end
