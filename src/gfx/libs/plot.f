
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

        subroutine plot(rin,sin,npin)
        save
c								**********
c	Move/draw lines on the output device.			*  PLOT  *
c								**********
c
c	NOTE:  Very little of this is actually device/site-specific,
c	but we don't want to include any more subroutine calls than are
c	really necessary, as this is already a very low-level routine.
c
        character*80 device
c
        common /plot sizes/ xsize,ysize
        common /plot device/ device,aspect,idev
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /plain font/ wid
        common /plot origin/ ro,so
        common /last point/ rl,sl
        common /dev status/ idevon,idevpen,idevwt
        common /dev init/ init
c
c	PostScript info:
c	---------------
c
	common /ps enforced/ ibounds,ps rmax,ps smax
        common /ps strokes/ nstroke,nstrpage,nstroketot
	common /ps bounding box/ ixleft,iybot,ixright,iytop
c
	parameter (NSMAX = 100)
c
c       SparcPrinter has problem with short lines:
c       -----------------------------------------
c
        common /sparcbug/ isp
	parameter (SPARCTOL = 0.25)
c
c       For use with idev = 1:
c	---------------------
c
        character*3 ich,jch
c
c       This is for the benefit of NCAR:
c       -------------------------------
c
        common /sysplt/ mmajx  ,mmajy  ,mminx  ,mminy  ,mxlab  ,mylab  ,
     &                  mflg   ,mtype  ,mxa    ,mya    ,mxb    ,myb    ,
     &                  mx     ,my     ,mtypex ,mtypey ,xxa    ,yya    ,
     &                  xxb    ,yyb    ,xxc    ,yyc    ,xxd    ,yyd    ,
     &                  xfactr ,yfactr ,xadd   ,yadd   ,xx     ,yy     ,
     &                  mfmtx(3)       ,mfmty(3)       ,mumx   ,mumy   ,
     &                  msizx  ,msizy  ,mxdec  ,mydec  ,mxor   ,mop(19),
     &                  mname(19)      ,mxold  ,myold  ,mxmax  ,mymax  ,
     &                  mxfac  ,myfac  ,modef  ,mf2er  ,mshftx ,mshfty ,
     &                  mmgrx  ,mmgry  ,mmnrx  ,mmnry  ,mfrend ,mfrlst ,
     &                  mcrout ,mpair1 ,mpair2 ,msblen ,mflcnt ,mjxmin ,
     &                  mjymin ,mjxmax ,mjymax ,mnxsto ,mnysto ,mxxsto ,
     &                  mxysto ,mprint ,msybuf(360)    ,mncpw  ,minst  ,
     &                  mbufa  ,mbuflu ,mfwa(12)       ,mlwa(12)       ,
     &                  mipair ,mbprs(16)      ,mbufl  ,munit  ,small
c
c       For the HP plotter:
c	------------------
c
        character*30 hpout
        common /pen posn/ npo
        common /mline on/ imline
c
c       For the X interface:
c
        external mcdxmove !$pragma C (mcdxmove)
        external mcdxdraw !$pragma C (mcdxdraw)
c
c       For the Tektronix options:
c	-------------------------
c
        character*1 ctrl(0:31),
     &              null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        common /ctrlch/ ctrl,
     &                  null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
c
        character*1 vec(0:5),up(3),down(3)
        data up/' ','L','F'/down/' ','L','G'/
c
        data ro/0./so/0./npo/3/rl,sl/0.,0./imline/0/
c
c-----------------------------------------------------------------------------
c
        up(1) = esc
        down(1) = esc
        vec(0) = gs
c
        if (init.eq.0) then
            init = -1
	    call noclear
            call mcinit
            call devon
            call clear
        end if
c
        r = rin
        s = sin
        npen = npin
        iloop = 0
c
10      np = abs(npen)
        if (idevon.eq.0) call devon
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        if (idev.eq.15) then
c
c	    SunCore
c	    -------
c
            if (np.eq.2) then
                call lineabs2(r+ro,s+so)
            else
                call moveabs2(r+ro,s+so)
            end if
            go to 998
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	else if (idev.eq.16) then
c
c	    PostScript
c	    ----------
c
	    rr = min(ps rmax,max(0.,nx0+xfac*(r+ro)))
	    ss = min(999.9,max(0.,ny0+yfac*(s+so)))
	    if (ibounds.ne.0) ss = min(ps smax,ss)
c
            if (np.eq.2) then
c
                if (isp.ne.0) then
c
c                   The SparcPrinter can't draw short lines!
c
                    if (xfac*abs(r-rl) + yfac*abs(s-sl).lt.SPARCTOL)
     &                      write(42,30)rr,ss,'p'
                else
                    write(42,30)rr,ss,'l'
30                  format(f7.3,f8.3,1x,a1,'%')
c                   
c                   NOTE: The trailing "%" is for identification purposes
c                   when erasing output...
c                   
                end if
            else if (np.eq.3) then
                write(42,30)rr,ss,'m'
            end if
c
c	    if (np.eq.2.or.ixleft.eq.10000) then
	    if (np.eq.2) then
		irr = rr
		iss = ss
		ixleft = min(ixleft,irr)
		iybot = min(iybot,iss)
		ixright = max(ixright,irr+1)
		iytop = max(iytop,iss+1)
	    end if
c
	    if (nstroke.ge.NSMAX) then
		call ps stroke
		write(42,30)rr,ss,'m'
	    end if
	    nstroke = nstroke+1
	    nstrpage = nstrpage+1
c
c           Only count "real" strokes in the total.
c
            if (np.eq.2) nstroketot = nstroketot + 1
	    go to 998
c
        else if (idev.eq.17) then
c
c	    X
c	    -
c
            if (np.eq.2) then
                call mcdxdraw(r+ro,s+so)
            else
                call mcdxmove(r+ro,s+so)
            end if
            go to 998
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	Determine scalings for other devices
c	------------------------------------
c
        i = nx0+xfac*(r+ro)
        j = ny0+yfac*(s+so)
        if (idev.ne.2) then
            i = max(0,min(nxpix,i))
            j = max(0,min(nypix,j))
        end if
        if (idev.ge.7)go to 400
        go to (100,200,400,400,500,501), idev
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	Output to "plot file"
c	---------------------
c
100     ip = npen
        if (ip.lt.0)ip = ip+3
        call reduce(i,ich)
        call reduce(j,jch)
        write(60,110)ip,ich,jch
110     format(i1,2a3)
        go to 998
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	NCAR
c	----
c
200     mx = max0(0,min0(i,32767))
        my = max0(0,min0(j,32767))
        ipen = 3-np
        minst = max0(0,min0(1,ipen))
c
        call put42
        go to 998
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	Tektronix/Versaterm-PRO
c	-----------------------
c
400     if ((idev.eq.7.or.idev.eq.8.or.idev.eq.11.or.idev.eq.12)
     &      .and.idevwt.gt.1) then
            i = i-idevwt/2
            j = j+idevwt/2
        end if
c
        if (nxpix.gt.1023) then
            i4 = i/4
            ii = i-4*i4
            j4 = j/4
            jj = j-4*j4
            vec(2) = char(96+ii+4*jj)
            i = i4
            j = j4
        else
            vec(2) = '`'
        end if
        j32 = j/32
        vec(1) = char(32+j32)
        vec(3) = char(96+(j-32*j32))
        i32 = i/32
        vec(4) = char(32+i32)
        vec(5) = char(64+(i-32*i32))
c
c       (These assignments are ok because none of the values are >127)
c
        if (idev.lt.13) then
            if (np.eq.3) then
                call type string(vec(0),6)
            else
                call type string(vec(1),5)
            end if
        else
            if (np.eq.3) then
                call type string(up,3)
            else
                call type string(down,3)
            end if
            call type string(vec(1),5)
        end if
        go to 998
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	HP plotter
c	----------
c
500     continue
501     if (imline.eq.1.and.r.eq.rl.and.s.eq.sl.and.npen.eq.3)
     &      go to 999
c
c       (Quick fix for possible use of plotin in mline)
c
        if (np.ne.npo) then
            if (np.eq.2) then
                write(hpout(1:3),510)'PD'
510             format(a2,';')
            else
                write(hpout(1:3),510)'PU'
            end if
            iout = 3
        else
            iout = 0
        end if
        write(hpout(iout+1:iout+14),520)i,j
520     format('PA',i5,',',i5,';')
        write(6,525)(hpout(k:k),k=1,14+iout)
525     format(1x,30a1)
c
c-----------------------------------------------------------------------------
c
c	End of routine -- clean up.
c	--------------------------
c
998     if (npen.lt.0) then
            ro = ro+r
            so = so+s
            rl = 0.
            sl = 0.
            call getlhe(rr)
            call setlhe(rr-r)
            call getbot(ss)
            call setbot(ss-s)
        else
            rl = r
            sl = s
        end if
        npo = np
999     if (iloop.le.0) return
        go to 1000
c
c	Variation -- combined move and draw.
c
        entry segment(r1,s1,r2,s2)
        iloop = 1
        npen = 3
        r = r1
        s = s1
        go to 10
c
1000    iloop = iloop-1
        if (iloop.lt.0) return
        npen = 2
        r = r2
        s = s2
        go to 10
c
        end
