
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

        subroutine plot details
        save
c
c	*****************************************************************
c	*								*
c	*	PLOT DETAILS:   Translate plot commands for		*
c	*			specific output devices.		*
c	*								*
c	*								*
c	*	"Essential" commands are:				*
c	*								*
c	*		weight		set pen weight			*
c	*		color		set pen color			*
c       *               background      set background color            *
c	*		devon/off	enable/disable output		*
c	*				to graphics device		*
c	*		graphin		get graphic input		*
c       *               gfxin           check if GIN is enabled         *
c	*		clear/erase	erase all or part of screen	*
c	*		invert		set foreground = background	*
c	*		point		plot a single point/pixel	*
c	*		polydraw	plot a polygon		        *
c	*		polyfill	plot a filled polygon		*
c	*		eop/newpage	take new page			*
c	*		mcquit		clean up and terminate graphics	*
c	*								*
c	*		num_win		number of open windows		*
c	*		set_win		set (X) window			*
c	*		kill_win	kill (X) window			*
c	*		curr_win	current (X) window		*
c	*		iconify_win	iconify current (X) window	*
c	*		iconify_all	iconify all (X) windows		*
c	*		win_read_line	get input from the current	*
c	*						   (X) window	*
c       *		win_flush	flush plotting output		*
c	*								*
c	*	Other routines (stored elsewhere):			*
c	*								*
c	*		mcinit		initialize graphics		*
c	*		plot		move/draw a line		*
c	*		setfill		set color for polyfill		*
c	*								*
c	*****************************************************************
c
	dimension rpoly(1),spoly(1)
c
        character*80 device
        character*10 chtekin,chout
        character*1 leadin
        character*1 tekin(10),tekout(10)
        character*1 vec(4)
        equivalence (chtekin,tekin),(chout,tekout)

	character*200 temp
c
        character*1 ctrl(0:31),
     &              null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        common /ctrlch/ ctrl,
     &                  null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        character*1 tekc(4),hpc(3),hpon(3),hpoff(6)
        common /dev controls/ tekc,hpc,hpon,hpoff
c
c	Internal plotting variables:
c
        common /dev status/ idevon,idevpen,idevwt
        common /dev details/ itek,ivers
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /last point/ rl,sl
        common /pen posn/ npo
        common /plain font/ wid
        common /plot device/ device,aspect,idev
        common /plot origin/ ro,so
        common /plot sizes/ xsize,ysize
        common /dev init/ init
	common /plot offset/ iin
	common /plot invert/ inv,ipensto
	save /plot invert/
        common /backgcolor/ ibackgcolor
        save /backgcolor/
c
c	HP options:
c       ----------
c
        common /hp plot/ ivdef
c
c	TEK options:
c       -----------
c
        common /input character/ leadin /input posn/ iposn,jposn
c
c	PS option:
c       ---------
c
	parameter (NSMAX = 100)
	common /ps strokes/ nstroke,nstrpage,nstroketot
	common /ps enforced/ ibounds,ps rmax,ps smax
	common /ps bounding box/ ixleft,iybot,ixright,iytop
	character*(*) id string
c
c       SUN options:
c       -----------
c
#ifdef HAS_SUNCORE
c        include "/usr/include/f77/usercore77.h"
        include "/usr/old/usercore77.h"
#endif
c
        common /graphproc/ igp
        common /mcpak_colormap/ ncolor,red(0:255),green(0:255),
     &                          blue(0:255)
	common /sunscreen/ isun
	common /findex/ index
        dimension xbox(4),ybox(4)
c
c       X options:
c       ---------
c
        logical gfxin
	character*(*) filename
c
	data ibounds/1/
c
c-----------------------------------------------------------------------
c
        entry getstatus(jgdevice,jgdevon,jgdevpen,jgdevwt)
c							    ***************
c	Return overall device status.			    *  GETSTATUS  *
c							    ***************
        jgdevice = idev
        jgdevon = idevon
        jgdevpen = idevpen
        jgdevwt = idevwt
c
        return
c
c-----------------------------------------------------------------------
c
        entry weight(iwt)
c							 	************
c       Adjust pen width/weight, if appropriate.	 	*  WEIGHT  *
c							 	************
c
c	CONVENTION:  Fiddle with details so that 1 is a "thin" line,
c	             10 is a "double weight" line, 20 is a "heavy"
c	             line, etc., independent of the device.
c
        if(idev.eq.1)then
            do while (iwt .ge.64)
                iwt = iwt - 64
            end do
            write(60,500)4,char(iwt+32)
500         format(i1,a1)
        else if((idev.eq.5.or.idev.eq.6).and.iwt.ne.idevwt)then
c
c	    HP GL (what does "weight" mean here?)
c
            if(idevon.le.0)write(6,510)hpon
510         format(1x,10a1)
            jpen = max(0,min(7,iwt))
            write(6,520)jpen+1,ivdef*(11-jpen)
520         format(' FS',i1,';VS',i2,';')
            if(idevon.le.0)write(6,510)hpoff
c
        else if((idev.eq.7.or.idev.eq.8.or.idev.eq.11.or.idev.eq.12)
     &          .and.iwt.gt.0)then
c
c	    Tek/Versaterm -- pensize specifies pen size in pixels.
c
	    iw = (iwt+4)/10+1
            call pensize(iw,iw)
c
        else if((idev.eq.15))then
c
c	    SunCore ideosyncrasies:
c
c		setlinewidth(x) -->	single for 0.0 to 0.3
c					double for 0.4, 0.5
c					heavy for 0.6
c					heavier for 0.7
c					heavier still for 0.8, 0.9
c					...
c
	    if (iwt.le.5) then
		ww = 0.1
	    else if (iwt.le.15) then
		ww = 0.4 + 0.01*(iwt-5)
	    else
		ww = 0.6 + 0.01*(iwt-15)
	    end if
            call setlinewidth(ww)
c
	else if(idev.eq.16)then
c
c	    PostScript.
c
	    call ps stroke
	    write(42,'(f5.1,'' setlinewidth'')').1*iwt
	    write(42,'(''2 setlinejoin 2 setlinecap'')')
c
c	    N.B. PS weight applies to SUBSEQUENT lines only.
c
        else if (idev.eq.17) then
            jwt = max(0, nint(.1*iwt+.5))
            call mcdxlinew(jwt)
        end if
c
        idevwt = iwt
c
        return
c
c-----------------------------------------------------------------------
c
        entry pen(ipenin)
c								***************
c       Set pen characteristics ("color").			*  PEN/COLOR  *
c								***************
        ipen = ipenin
525     if(idev.eq.1)then
            do while (ipen .ge.64)
                ipen = ipen - 64
            end do
            write(60,500)5,char(ipen+32)
        else if((idev.eq.5.or.idev.eq.6).and.ipen.ne.idevpen)then
            if(idevon.le.0)write(6,510)hpon
            write(6,530)max(0,min(8,ipen))
530         format(' SP',i1,';')
            if(idevon.le.0)write(6,510)hpoff
            idevpen = ipen
        else if(idev.eq.13.or.idev.eq.14)then
            call output char(esc)
            write(chout,531)min(7,max(0,ipen))
531         format('ML',i1)
            call type string(tekout,3)
        else if(idev.eq.15)then
            call setlineindex(ipen)
            call settextindex(ipen)
c           Note that ipen = 0 erases
	else if(idev.eq.16)then
	    call ps color(ipen)
c
c	    N.B. PS pen color applies to SUBSEQUENT lines only.
c
        else if (idev.eq.17) then
c
c           For now, just wrap the colormap...
c
            jpen = ipen/ncolor
            if (jpen.gt.0) ipen = ipen - ncolor*jpen
            call mcdxcolor(ipen)
        end if
        idevpen = ipen
c
        return
c
        entry color(icolor)
c
        ipen = icolor
        go to 525
c
	entry getcolor(jcolor)
	jcolor = idevpen
	return
c
c-----------------------------------------------------------------------
c
        entry background(iback)
c							       ****************
c       Set background color.			               *  BACKGROUND  *
c							       ****************
        if (idev.eq.15) then
            call corebackg(iback)
        else if (idev.eq.17) then
c
c           For now, just wrap the colormap...
c
            jback = iback/ncolor
            if (jback.gt.0) iback = iback - ncolor*jback
            call mcdxbackg(iback)
        end if
c
        ibackgcolor = iback
c
        return
c
c-----------------------------------------------------------------------
c
        entry devon
c								***********
c       Explicitly initiate output to graphics device.		*  DEVON  *
c								***********
        if (idevon.ne.1) then
            if(itek.eq.1)then
                if(ivers.eq.0)then
                    if (idev.eq.9.or.idev.eq.10) then
c                       
c                       Xterm escape sequences:
c                       
                        call output char(esc)
                        call output char('[')
                        call output char('?')
                        call output char('3')
                        call output char('8')
                        call output char('h')
                    else
                        call output char(esc)
                        call output char(';')
                        write(6,*)
                    end if
                end if
                call output char(gs)
            else if(idev.eq.5.or.idev.eq.6)then
                write(6,510)hpon
            end if
            idevon = 1
        end if
        return
c
c-----------------------------------------------------------------------
c
        entry devoff
c								************
c       Explicitly terminate output to graphics device.		*  DEVOFF  *
c								************
        iquit = 0
c
200     if (idevon.ne.0) then
c
            if (itek.eq.1) then
                if (ivers.eq.0) then
                    if (idev.eq.9.or.idev.eq.10) then
c                       
c                       Xterm escape sequences:
c                       
c                       Flush the buffer (100 may be overkill, but...).
c                       
                        call plot(-10.,-10.,3)
                        do 210 i=1,100
                            call plot(-10.,-10.,2)
210                     continue
c                       
                        call output char(esc)
                        call output char(char(3))
c                       
                    else
                        call plot(-10.,-10.,3)
                        call plot(-10.,-10.,2)
                    end if
                end if
                call output char(esc)
                call output char('2')
c
            else if(idev.eq.5.or.idev.eq.6)then
                write(6,510)hpoff
                npo = 3
            end if
c
            idevon = 0
c
        end if
c
        if (iquit.eq.1) then
c
c           Complete MCQUIT operations.
c
            init = 0
            idev = 0
            device = ' '
            iin = 0
        end if
c
        return
c
c-----------------------------------------------------------------------
c
        entry graphin(r,s)
c
c       Return the location of the graphics pointer relative	*************
c	to the current origin (implemented for Tektronix, X,	*  GRAPHIN  *
c	and SUN options only).					*************
c
        if (.not.gfxin()) return
c
        if(idev.eq.15)then
            kount = 0
400         call awtbuttongetloc2(0,1,ibutton,xx,yy)
            if(ibutton.eq.0)then
                kount = kount+1
                if(kount.lt.12000)then
                    call uwait(5000)
                    go to 400
                end if
                write(6,*)'timed out...'
                xx = 0.
                yy = 0.
            end if
c
            r = xsize*xx-ro
            s = xsize*yy-so
c           Yes, this SHOULD be xsize!!
c
            iposn = 1000*xx
            jposn = 1000*yy
            return
        else if (idev.eq.17) then
            call mcdxgin(r,s)
            iposn = nxpix*r/xsize
            jposn = nypix*s/ysize
            r = r - ro
            s = s - so
            return
        end if
c
        if (itek.eq.0) return
c
c       The rest is for Tektronix mode(s).
c
        if (idevon.eq.0) then
            if (idev.eq.9.or.idev.eq.10) then
                call output char(esc)
                call output char('[')
                call output char('?')
                call output char('3')
                call output char('8')
                call output char('h')
            else
                call output char(gs)
            end if
            idevon = 1
        end if
c
c       enter gin mode...
c
        call output char(esc)
        call output char(ctrlz)
c
c       ...read the pointer address, preceded by a lead-in character...
c
        if(nxpix.gt.1023)then
            kmax = 6
        else
            kmax = 5
        end if
c        do 95 k=1,kmax
        do 95 k=1,kmax
            call input char(tekin(k))
            if((tekin(k).eq.lf.or.tekin(k).eq.cr).and.k.gt.1)go to 96
95      continue
        k = kmax+1
96      kmax = k-1
c
c       ...and exit gin mode.
c
        if (idev.ne.9.and.idev.ne.10) then
            call output char(cr)
        else
c
c           Xterm version -- leave GFX mode:
c
            call output char(esc)
            call output char(char(3))
            idevon = 0
        end if
c
c       Each use of graphin causes the tektronix to go out of graph
c       mode momentarily, so we must reposition the pointer:
c
        call plot(rl,sl,3)
c
c	This will also undo VersaTerm "inverted" mode, so reset, if necessary:
c
        if(inv.eq.1)then
            call output char(esc)
            if(idev.eq.14)then
                chout = 'ML0'
                call type string(tekout,3)
            else
                call output char(del)
            end if
	end if
c
c       Now decode the address from tekin(2),...,tekin(kmax).
c       note that the input string is encoded differently from the
c       corresponding output in plot or pixel.
c
        i32 = 32
c	(force integer*4 arithmetic below)
c
        if(kmax.eq.5)then
            i = (ichar(tekin(2))-32)*i32 + (ichar(tekin(3))-32)
            j = (ichar(tekin(4))-32)*i32 + (ichar(tekin(5))-32)
            if(nxpix.gt.1023)then
                i = 4*i
                j = 4*j
            end if
        else
c
c           This possibility might not actually be legal. It is
c           not yet tested.
c
            ii = ichar(tekin(5))-32
            jj = ii/4
            ii = ii-4*jj
            i = 4*((ichar(tekin(2))-32)*i32 + (ichar(tekin(3))-32)) + ii
            j = 4*((ichar(tekin(4))-32)*i32 + (ichar(tekin(6))-32)) + jj
        end if
        r = (i-nx0)/xfac-ro
        s = (j-ny0)/yfac-so
        iposn = i
        jposn = j
        leadin = chtekin(1:1)
        if(idevon.le.0)then
           call output char(esc)
           call output char('2')
        end if
        return
c
c-----------------------------------------------------------------------
c
        entry erase(rmin,rmax,smin,smax)
c								***********
c       Clear a rectangular region of the display (for		*  ERASE  *
c	SUN, X, PostScript and "Versaterm" options only)	***********
c
        if(idev.eq.15.or.idev.eq.16.or.idev.eq.17)then
            xbox(1) = rmin
            xbox(2) = rmax
            xbox(3) = rmax
            xbox(4) = rmin
            ybox(1) = smin
            ybox(2) = smin
            ybox(3) = smax
            ybox(4) = smax
	    call polyerase(xbox,ybox,4)
        end if
c
        if(ivers.eq.0)return
c
        if(idevon.le.0)then
            call output char(gs)
            idevon = 1
        end if
c
        if(1.eq.0.and.idev.ge.8.and.idev.le.12)then
c
c           if pensize doesnt work, use this...
c
            i1 = (rmin+ro)*xfac+nx0
            i2 = (rmax+ro)*xfac+nx0
            j1 = (smin+so)*yfac+ny0
            j2 = (smax+so)*yfac+ny0
            call output char(esc)
            call output char(del)
            if(i2-i1.gt.j2-j1)then
                do 11001 j=j1,j2
                    call pixel(i1,j,3)
                    call pixel(i2,j,2)
11001           continue
            else
                do 11002 i=i1,i2
                    call pixel(i,j1,3)
                    call pixel(i,j2,2)
11002           continue
            end if
            call output char(cr)
            return
        else if(idev.lt.13)then
            i = 0
            j = (smax-smin)*yfac+1
            call pensize(i,j)
c
c           turn on erase mode:
c
            call output char(esc)
            call output char(del)
c
c           draw the line:
c
            call plot(rmin,smax,3)
            call plot(rmax,smax,2)
c
c           reset:
c
            call output char(cr)
        else if(idev.eq.14)then
            call output char(esc)
            chout = 'MP0'
            call type string(tekout,3)
            call output char(esc)
            chout = 'LP'
            call type string(tekout,2)
            call xyencode(rmin,smin,vec)
            call type string(vec,5)
            call type string(tekout,2)
            call plot(rmin,smax,2)
            call plot(rmax,smax,2)
            call plot(rmax,smin,2)
            call plot(rmin,smin,2)
            call output char(esc)
            chout = 'LE'
            call type string(tekout,2)
        end if
        return
c
c-----------------------------------------------------------------------
c
        entry clear
c								***********
c       clear entire screen, if any,...				*  CLEAR  *
c								***********
        if (idev.eq.1) then
            write(60,500)6
        else if (idev.eq.15) then
            call plframe
	else if (idev.eq.16) then
c
c	    Old version:
c
c	    write(42,'(''erasepage'')')
c
c	    New version (saves unnecessary output to the printer):
c	    Close the output file, and copy it to a new file, omitting
c	    all "plot" commands.
c
	    write(42,'(/''%Re-writing...''/)')
	    call ps rewrite
            ixleft = 10000
            iybot = 10000
            ixright = 0
            iytop = 0
c
	    nstroke = 0
	    nstrpage = 0
	    nstroketot = 0
        else if (idev.eq.17) then
            call mcdxclear
        else if (itek.eq.1) then
            if(idev.eq.12.or.idev.eq.14)then
                call output char(esc)
                call output char(ff)
            else
                call type string(tekc,4)
            end if
            if(ivers.eq.0)then
                call plot(-10.,10.,3)
                call plot(-10.,10.,2)
            end if
        end if
        return
c
c-----------------------------------------------------------------------
c
        entry invert
c								************
c       Toggle between draw and undraw modes (SUN,		*  INVERT  *
c	PostScript, X and Versaterm modes only).		************
c
        if (idev.eq.15.or.idev.eq.16.or.idev.eq.17) then
            if(idevpen.gt.0)then
                ipensto = idevpen
                idevpen = 0
            else
                idevpen = ipensto
            end if
c
	    if(idev.eq.15)then
        	call setlineindex(idevpen)
	    else if (idev.eq.16) then
		call ps color(idevpen)
            else if (idev.eq.17) then
                call mcdxcolor(idevpen)
	    end if
            return
        end if
c
        if (idev.ne.7.and.idev.ne.8.and.idev.ne.10
     &       .and.idev.ne.12.and.idev.ne.14) return
c
        if(idevon.le.0)then
            call output char(gs)
            idevon = 1
        end if
        if(inv.eq.0)then
            call output char(esc)
            if(idev.eq.14)then
                chout = 'ML0'
                call type string(tekout,3)
            else
                call output char(del)
            end if
            inv = 1
        else
            if(idev.eq.14)then
                call output char(esc)
                chout = 'ML1'
                call type string(tekout,3)
            else
                call output char(cr)
            end if
            inv = 0
        end if
c
c       Note that we have probably lost track of the graphics pointer
c       through this action.
c
        return
c
c-----------------------------------------------------------------------
c
	entry point(r,s)
c								***********
c	Plot a single point.					*  POINT  *
c								***********
	if(idev.eq.16)then
	    rr = nx0+xfac*(r+ro)
	    ss = ny0+yfac*(s+so)
	    if(ibounds.ne.0)then
		rr = max(0.,min(ps rmax,rr))
		ss = max(0.,min(ps smax,ss))
	    end if
	    write(42,'(f7.3,f8.3,'' p%'')')rr,ss
	    if(nstroke.ge.NSMAX/2)then
		call ps stroke
		write(42,'(f7.3,f8.3,'' m%'')')rr,ss
	    end if
	    nstroke = nstroke+1
	    nstrpage = nstrpage+1
	    nstroketot = nstroketot+1
        else if (idev.eq.17) then
            call mcdxpoint(r+ro,s+so)
	else
	    call plot(r,s,3)
	    call plot(r,s,2)
	end if
c
	return
c
c-----------------------------------------------------------------------
c
	entry polydraw(rpoly,spoly,npoly)
c								**************
c	Draw a regular polygon.		          		*  POLYDRAW  *
c								**************
        call plot(rpoly(npoly),spoly(npoly),3)
        do 10 i = 1,npoly
	    call plot(rpoly(i),spoly(i),2)
10      continue
c
	return
c
c-----------------------------------------------------------------------
c
	entry polyfill(rpoly,spoly,npoly)
c								**************
c	Draw a filled regular polygon.				*  POLYFILL  *
c								**************
	if((idev.ne.15.and.idev.ne.16.and.idev.ne.17)
     &          .or.index.lt.0) return
c
        if (idev.eq.15.or.idev.eq.17) then
c
c           SunCore/X want absolute coordinates.
c
	    do 110 i=1,npoly
		rpoly(i) = rpoly(i) + ro
		spoly(i) = spoly(i) + so
110	    continue
        end if
c
	if(idev.eq.15)then
c
	    call polygonabs2(rpoly,spoly,npoly)
c
	else if (idev.eq.16) then
c
	    call ps color(index)
c
	    call plot(rpoly(npoly),spoly(npoly),3)
	    do 120 i=1,npoly
                call plot(rpoly(i),spoly(i),2)
120         continue
	    write(42,'(''fill'')')
	    call ps color(ipen)
c
        else if (idev.eq.17) then
c
            call mcdxpolyf(rpoly,spoly,npoly,index)
c
	end if
c
        if (idev.eq.15.or.idev.eq.17) then
	    do 130 i=1,npoly
		rpoly(i) = rpoly(i) - ro
		spoly(i) = spoly(i) - so
130         continue
        end if
c
        return
c
c-----------------------------------------------------------------------
c
	entry bounded
	entry bounds
c								*************
c	Enforce/dont enforce xsize x ysize bounds		*  BOUNDED  *
c	on (currently, only PostScript) output.			*************
c
	ibounds = 1
	return
c
	entry unbounded
	entry nobounds
	ibounds = 0
	return
c
c-----------------------------------------------------------------------
c
        entry eop
	entry newp
	entry newplot
	entry nextplot
	entry newpage
c								***************
c       End plot, prepare for new one (flush buffer, 		*  EOP, etc.  *
c	clear screen, etc.),...					***************
c
        if(idev.eq.1)then
            write(60,500)7
        else if(idev.eq.2)then
            call frame
        else if(itek.eq.1)then
            call display text(
     &      'Enter <CR> to clear screen.             ',40)
            read(5,350)idummy
350         format(a1)
            if(idev.ne.13.and.idev.ne.14)then
                call type string(tekc,4)
            else
                call output char(esc)
                call output char(ff)
            end if
        else if(idev.eq.5.or.idev.eq.6)then
            if(idevon.le.0)write(6,510)hpon
            write(6,510)hpc
            if(idevon.le.0)write(6,510)hpoff
        else if(idev.eq.15)then
	    call sunwait
            call plframe
	else if(idev.eq.16)then
	    call ps page(1)
        else if (idev.eq.17) then
            write(6,'(a)')
     &      'Click right-hand mouse button in window to exit X.'
            call mcdxidle
            call mcdxreset
        end if
        go to 200
c
c-----------------------------------------------------------------------
c
        entry set_win(iwin, iret)
c								***************
c	Choose output window.					*  SET_WIN    *
c								***************
c
	if (idev.eq.17) then
	    iret = mcdxsetwin(iwin)
c
c           Return status is -1 on error.
c
	end if
c
	return
c
c-----------------------------------------------------------------------
c
        entry kill_win(iwin)
c								***************
c	Kill output window.					*  KILL_WIN   *
c								***************
c
	if (idev.eq.17) then
	    call mcdxkillwin(iwin)
            if (mcdxnopen().le.0) then
c
c               No windows left!
c
                idev = 0
                device = ' '
                iin = 0
                idevon = 0
            end if
	end if
c
	return
c
c-----------------------------------------------------------------------
c
        entry iconify_win(iwin)
c								 ***************
c	Iconify an output X window.				 * ICONIFY_WIN *
c								 ***************
c
	call mcdxiconifywin(iwin)
c
	return
c
c-----------------------------------------------------------------------
c
        entry iconify_all
c								 ***************
c	Iconify all output X windows.				 * ICONIFY_ALL *
c								 ***************
c
	call mcdxiconifyall()
c
	return
c
c-----------------------------------------------------------------------
c
        entry set_colormap(filename,ier)
c								 ***************
c	Load a colormap from a file.				 * SET_COLORMAP*
c								 ***************
c
	nf = 0
	do i=len(filename),1,-1
	    if (nf.eq.0.and.filename(i:i).gt.' ') nf = i
	end do
c
c       HP fortran doesnt understand '\0'...
c
	do i=1,nf
	    temp(i:i) = filename(i:i)
        end do
        temp(nf+1:nf+1) = char(0)
	call mcdxsetwincolormap(temp(1:nf+1), ncolor, ier)
c
	return
c
c-----------------------------------------------------------------------
c
        entry win_read_line(string)
c							       *****************
c	Read a line of text from the (X) window system.	       * WIN_READ_LINE *
c							       *****************
c
	if (mcdxnopen().gt.0) then
	    call mcdxrdline(string)
	end if
c
	return
c
c-----------------------------------------------------------------------
c
        entry win_flush()
c							           *************
c	Flush graphics output in  the (X) window system.           * WIN_FLUSH *
c							           *************
c
	if (mcdxnopen().gt.0) then
	    call mcdxflush()
	end if
c
	return
c
c-----------------------------------------------------------------------
c
        entry routine id(id string)
c								****************
c	Identify key routines.					*  ROUTINE ID  *
c								****************
c
	if (idev.eq.16) then
	    do 450 i=len(id string),1,-1
		if (id string(i:i).gt.' ') go to 455
450	    continue
	    return
455	    write(42,'(''%%% Routine '',a,'':'')')id string(1:i)
	end if
c
	return
c
c-----------------------------------------------------------------------
c
        entry mcquit
        iquit = 0
        go to 1000
c
        entry mcquit1
        iquit = 1
c								************
c	Terminate graphics session.				*  MCQUIT  *
c								************
1000    if (isun.gt.0) then
	    write(6,'(''Enter <CR> to terminate graphics.'')')
	    call plstop
	    isun = 0
	end if
c
c	NOTE: The system call in psquit apparently requires that SunCore
c	graphics be terminated before attempting to print any PostScript file.
c
	call psquit(2)
c
c	Call to psquit does nothing if no PS device is open.
c
	if (mcdxnopen().gt.0) then
            if (iquit.eq.1) then
		write(6,'(a)')
     $      'Click right-hand mouse button in window to exit X.'
                call mcdxidle
            end if
        end if
c
        call mcdxquit
c
c       Remaining termination is done by devoff (label 200).
c
        iquit = 1
        go to 200
c
	end


	logical function gfxin()
        save
c                                                                ***********
c	Return TRUE iff graphin will work.                       *  GFXIN  *
c                                                                ***********
        character*80 device
        common /plot device/ device,aspect,idev
        common /dev details/ itek,ivers
c
	if (idev.eq.15.or.idev.eq.17.or.itek.eq.1
c     &          .or.(itek.eq.1.and.idev.ne.9.and.idev.ne.10)
     &          ) then
	    gfxin = .true.
	else
	    gfxin = .false.
	end if
c
	end


	integer function num_win(id)
        save
c                                                                ***********
c	Return the number of open windows associated with id.    * NUM_WIN *
c                                                                ***********
        character*80 device
        common /plot device/ device,aspect,idev
c
        logical ps_open
        external mcdxnopen	!$pragma C (mcdxnopen)
c
        num_win = 0
	if (id.gt.0) then
c
c           PostScript and X can be open without being active.
c
            if (id.eq.16) then
                if (ps_open()) num_win = 1
            else if (id.eq.17) then
                num_win = mcdxnopen()
            else
                if (id.eq.idev) num_win = 1
            end if
        end if
c
	end


	integer function curr_win()
        save
c                                                                ************
c	Return the current output window.                        * CURR_WIN *
c                                                                ************
        character*80 device
        common /plot device/ device,aspect,idev
c
        external mcdxcurrwin	!$pragma C (mcdxcurrwin)
c
	if (idev.le.0) then
	    curr_win = -1
        else
            curr_win = 0
            if (idev.eq.17) curr_win = mcdxcurrwin()
        end if
c
	end


	subroutine polyerase(r,s,n)
        save
c
c	Erase the specified polygonal area.
c
        external mcdxpolyf	!$pragma C (mcdxpolyf)
c
	dimension r(1),s(1)
        character*80 device
        common /plot device/ device,aspect,idev
        common /plot origin/ ro,so
        common /backgcolor/ ibackgcolor
        save /backgcolor/
c
	if(idev.eq.15)then
c
c           Switch to absolute coordinates:
c
	    do k=1,n
		r(k) = r(k) + ro
		s(k) = s(k) + so
            end do
c
            call setfillindex(0)
            call polygonabs2(r,s,n)
c
	    do k=1,n
		r(k) = r(k) - ro
		s(k) = s(k) - so
            end do
c
	else if(idev.eq.16)then
c
	    call ps stroke
	    call plot(r(n),s(n),3)
	    do i=1,n
                call plot(r(i),s(i),2)
            end do
c
c           PS background color is WHITE (grey 1.0)...
c
	    write(42,'(''1.0 setgray''/''fill''/f5.3,'' setgray'')')
     &		    ipen/256.
	    call ps color(ipen)
c
        else if (idev.eq.17) then
c
c           Switch to absolute coordinates:
c
	    do k=1,n
		r(k) = r(k) + ro
		s(k) = s(k) + so
            end do
c
            call mcdxpolyf(r,s,n,ibackgcolor)
c
	    do k=1,n
		r(k) = r(k) - ro
		s(k) = s(k) - so
            end do
c
	end if
c
	end
