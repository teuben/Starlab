
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
c         *************************************************************
c         *                                                           *
c         *    Initialization of mcpak and interface between it       *
c         *    and any "local" plotting package:                      *
c	  *							      *
c         *    ROUTINES mcinit(a), dev select(a), get device, weight, *
c         *    pen (color), dev on, dev off, graphin, erase, clear,   *
c	  *    invert, polyfill, bounded, eop, mcquit, display text   *
c	  *    and plot.				              *
c         *                                                           *
c         *    Also included are numbr and symbl, which use the       *
c         *    external number and symbol routines (if any) to mimic  *
c         *    the nomber and simbol calling conventions.             *
c         *                                                           *
c         *    For portability, no other plotting routine should      *
c         *    ever call the external routines directly. instead,     *
c         *    use plot, numbr, symbl, clear, etc.                    *
c         *                                                           *
c         *                                                           *
c         *    (The above-named routines are now split among          *
c         *     several source file.)                                 *
c         *                                                           *
c         *************************************************************
c
        subroutine mcinit
        save
c
c       Initialize device, clear screen, establish plotting
c       area, aspect ratio, etc.
c
        character*80 device
        character *(*) dev
        common /plot sizes/ xsize,ysize
        common /plot device/ device,aspect,idev
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /ncar/ nxpix1,nypix1,nx01,ny01,xfac1,yfac1
        common /plain font/ wid
        common /dev status/ idevon,idevpen,idevwt
        common /dev details/ itek,ivers
	common /sunscreen/ isun
	common /plot offset/ iin
c
        logical set
c
        common /dev init/ init
        data init/0/
c
c
c       Entry points:
c
c		mcinit:	 set up graphics, prompt for device
c					  (or get from environment)
c		mcinita: set up graphics, take device as argument
c
c		devselect:  don't change graphic state, prompt for device
c						  (or get from environment)
c		devselecta: don't change graphic state, take device as
c  								argument
c
c-----------------------------------------------------------------------
c
        device=' '
        go to 50
c
c-----------------------------------------------------------------------
c
        entry mcinita(dev)
        device=dev
c
50      if (init.gt.0) return
	call plot setup
c
        iclr=1
        if (init.lt.0) iclr=0
        init=1
        iin=1
        call init chars
        go to 100
c
c-----------------------------------------------------------------------
c
        entry dev selecta(dev)
        device=dev
        go to 75
c
c-----------------------------------------------------------------------
c
        entry dev select
        device=' '
c
c-----------------------------------------------------------------------
c
75      iin=2
	init = 1
        call devoff
c
c       Most devices are simply "on" or "off" and have only one display.
c       Some combinations of devices are allowed, others are forbidden.
c
c       Presently, allow an arbitrary number of X windows,
c       	   only have one PostScript file open at a time (but
c			do allow it to remain open when X is selected),
c       	   on return to an open PostScript file, append to it,
c       	   don't allow any other output channel if "Sun" is used.
c
c       These combinations are mainly handled in getdevice.
c
        if (idev.eq.5.or.idev.eq.6) then
c
c           Clean up the HP plotter.
c
            write(6,80)27,27
80          format(1x,a1,'.Y',a1,'.L')
            read(5,*)idummy
            call devoff
c
        end if
c
100     idevwt = 1
	if (idevpen.lt.0) idevpen = 0
c
c       Determine the new device.
c
	call get device
c
c       Restore old settings:
c
	call weight(idevwt)
	if (idevpen.gt.0) call pen(idevpen)
c
c	(These is necessary because idevwt, idevpen may be set by getdevice.)
c
        if (itek.eq.1.and.iclr.eq.1) then
            call devon
            call clear
            call devoff
        end if
c
c	Force ysize/xsize = aspect ratio, so coordinate angles = physical ones.
c
        ysize=xsize*aspect
        xfac=nxpix/xsize
        yfac=nypix/ysize
        if (idev.eq.2) then
            xfac1=nxpix1/xsize
            yfac1=nypix1/ysize
        end if
c
c	Be careful setting up default "frame" options, as these
c       may reset the graphics in an unwanted manner.
c
        if (iin.ne.1) return
c
        call gethset(set)
        if (.not.set) call sethts(.25,.2)
c
        call getmset(set)
        if (.not.set) call setmod(0,1,0,0,0)
c
        call setlhe(0.)
c
        call getfset(set)
        if (.not.set) then
            call setxtf
            if (itek.eq.1) call setpln
        end if
c
        end
