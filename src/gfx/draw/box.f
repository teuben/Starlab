
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

c------------------------------------------------------------------------
c
c       Box creation and resizing.
c
c       Contents:   setup
c                   newbox
c                   offbox
c                   popbox
c                   setaspects
c                   fitbox
c
c------------------------------------------------------------------------


        block data init box params
        save
c
        common /last aspect/ aspect0
        common /compress frames/ icompress
c
        data aspect0/1./icompress/0/
c
        end


        subroutine setup(device,id)
        save
c
c       Set up graphics and frame parameters.
c       
        character*(*) device,id
        common /plot sizes/ xsize,ysize
        common /draw params/ roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &                       idevset,jbox,iorig
c       
        common /initial/ roff0,soff0,hn0,hs0,hp0
        common /last aspect/ aspect0
        common /dev init/ init
c
        character*80 dev
        common /plot device/ dev,aspect,idev
        character*80 colormapfile
        common /mcdraw_colormap_file/ colormapfile
c
        common /prompt/ iprompt
        integer curr_win
c
        if (idevset.eq.0.and.init.eq.0) then
c
c           Establish some defaults and initialize the system.
c
            if (iprompt.eq.1) then
                if (id(1:1).ne.' ') then
                    write(6,'(a,a,a)')'Graphics setup triggered by "',
     $                        id,'" command'
                else
                    if (device(1:1).eq.' ') then
                        write(6,'(a)')'Initializing graphics with "de"'
                    else
                        write(6,'(a,a,a)')
     $                      'Initializing graphics with "de ',
     $                      device,'"'
                    end if
                end if
                write(6,'(a)')'Prior settings will be overridden'
            end if
c
            call options('-x 0.5 -b 255')
            call mcinita(device)
            call setpln
            call sethts(hs,hn)
            call nobounds
            idevset = 1
            colormapfile = 'default'

        else
c
c           Choose a new device.
c
c           NOTE: Context switching must be performed here if we want 
c                 different windows/devices to be independent.
c
            idevold = idev
            icurrwin = curr_win()
            call save_context(idev,icurrwin)
c
            call dev selecta(device)
            colormapfile = 'default'
c
            icurrwin = curr_win()
            call load_context(idev,icurrwin)
c
            idevset = 1
c
            if (idevold.eq.17.and.idev.ne.idevold
     $              .and.num_win(17).gt.0) then
c
c               Shouldn't be necessary to close X windows if we can
c               get input from them...
c
c                if (iprompt.ne.0) write(6,*)
c     $                  '(Iconifying inactive X-windows.)'
c                call iconify_all
            end if
        end if
c
        call getstatus(idum,idum,ic, iw)
c
c       The aspect ratio will have been set to the ratio for the
c       specified device.  Reset it here if aspect1 has been set.
c
        if (aspect1.gt.0.) aspect = aspect1
c
        call setaspects(aspect,ax,ay)
c
c       On return from setaspects, aspect will be the specified value or
c       the device default, and ax and ay will be set accordingly.
c
        if (iorig.eq.0) then
c
c           Initial setup of basic graphics parameters.
c
c           The initial definitions of roff0 and soff0 are in terms
c           of the device dimensions.  Because of the possibility of
c           horizontal y-labels, roff0 is somewhat larger than soff0.
c
            roff0 = xsize*roff0
            soff0 = ysize*soff0
c
c           NOTE: The height settings in effect at the first frame
c           initialization become the default.
c
            hs0 = hs
            hn0 = hn
            hp0 = hp
c
            iorig = 1
c
        end if
c
c       Establish the new offset of the box within the frame,
c       and other plot parameters.
c
        if (jbox.eq.0) then
            scale = 1.
            roff = roff0
            soff = soff0
            hs = hs0
            hn = hn0
            hp = hp0
        else
            scale = .5
            roff = roff0 / 1.5
            soff = soff0 / 1.25
            ss = sqrt(scale)
            hs = hs0 * ss
            hn = hn0 * ss
            hp = hp0 * ss
        end if
c
        call sethts(hs,hn)
c
c       Choose xlen and ylen, maintaining the correct aspect ratio,
c       and fitting into the space available.
c
        call fitbox(roff,soff,aspect,scale,xlen,ylen,xextra,yextra)
c
c       Establish the offset of the frame.  Note that the subdivision
c       splits the display into four equal pieces, unless the "compress"
c       flag has been sent (via common) to fitbox, and xextra and yextra
c       are nonzero.
c
        call setorigin(0.,0.)
        if (jbox.eq.2.or.jbox.eq.4) call plot(.5*xsize-xextra,0.,-3)
        if (jbox.eq.1.or.jbox.eq.2) call plot(0.,.5*ysize-yextra,-3)
c
c       Offset the box within the frame.
c
        call plot(roff,soff,-3)
        call setlhe(-roff)
        call setbot(-soff)
c
c       Save the current value of the aspect ratio.
c
        aspect0 = aspect
c       
        end
        
        
	subroutine newbox(input,nin,ierbox,*)
        save
c       
c	Choose a new output area.
c       The input arguments from mcdraw are in input(1:nin).
c       
	character*(*) input
	dimension ierbox(0:4)
c       
        common /plot sizes/ xsize,ysize
        common /draw params/ roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &          idevset,jbox,iorig
c       
	common /initial/ roff0,soff0,hn0,hs0,hp0
        common /last aspect/ aspect0
c       
        data inew/0/
c       
        jboxo = jbox
        idel = 1
c       
        if (nin.lt.2.or.input(2:2).eq.' '.or.input(2:3).eq.'- '
     &          .or.input(2:3).eq.'-0') then
            jbox = 0
            idel = 0
        else
            read(input(2:nin),'(i5)',iostat=io)jbox
            if (io.ne.0) then
                write(6,'(a,a,a)')
     $              'Error reading argument "',input(2:nin),'"'
                go to 1001
            end if
        end if
c       
        if (jbox.lt.0) then
            jbox = -jbox
            idel = 0
        end if
c       
        if (jbox.eq.5) then
            jbox = 0 
            idel = 0
        end if
c       
        jbox = max(0,min(4,jbox))
c       
        if (aspect1.gt.0.) aspect = aspect1
        call setaspects(aspect,ax,ay)
c       
        if (inew.eq.0.or.aspect.ne.aspect0.or.jbox.ne.jboxo) then
c           
            inew = 1
c           
            if (idevset.eq.0) then
                jb = jbox
                jbox = jboxo
                call setup(' ','newbox')
                jbox = jb
                aspect = aspect1
                call setaspects(aspect,ax,ay)
            end if
c           
            if (aspect.ne.aspect0.or.jbox*jboxo.eq.0) then
c               
c               Change plot size or shape.  Note that "=x" has the
c               effect of reinitializing the offsets and heights to
c               the default, or scaled default, values.
c               
                if (jbox.eq.0) then
                    scale = 1.
                    roff = roff0
                    soff = soff0
                    hs = hs0
                    hn = hn0
                    hp = hp0
                else
                    scale = .5
                    roff = roff0 / 1.5
                    soff = soff0 / 1.25
                    ss = sqrt(scale)
                    hs = hs0 * ss
                    hn = hn0 * ss
                    hp = hp0 * ss
                end if
c               
            end if
c           
            call sethts(hs,hn)
c           
        end if
c
c       Choose xlen and ylen, maintaining the correct aspect ratio,
c       and fitting into the space available.
c           
        call fitbox(roff,soff,aspect,scale,xlen,ylen,
     &              xextra,yextra)
c           
c       Reposition the bottom-left corner of the current frame.
c       Note that the subdivision splits the display into four
c       equal pieces unless the "compress" flag has been sent
c       (via common) to fitbox, and xextra and yextra are nonzero.
c           
        call setorigin(0.,0.)
        if (jbox.eq.2.or.jbox.eq.4)
     &          call plot(.5*xsize-xextra,0.,-3)
        if (jbox.eq.1.or.jbox.eq.2)
     &          call plot(0.,.5*ysize-yextra,-3)
c           
c       The old box will only be erased if both idel and ierbox are
c       set, or if a shape change has occurred.
c           
        if (idel.eq.1.and.(ierbox(jbox).eq.1
     &          .or.aspect.ne.aspect0)) then
            if (jbox.eq.0) then
                call clear
            else
c               
c               Clear any unused space.
c
                if (jbox.eq.1) then
                    call erase(0.,.5*xsize-xextra,0.,.5*ysize+yextra)
                else if (jbox.eq.2) then
                    call erase(0.,.5*xsize+xextra,0.,.5*ysize+yextra)
                else if (jbox.eq.3) then
                    call erase(0.,.5*xsize-xextra,0.,.5*ysize-yextra)
                else
                    call erase(0.,.5*xsize+xextra,0.,.5*ysize-yextra)
                end if
c
            end if
        end if
c       
        do k=0,4
            ierbox(k) = 1
        end do
        ierbox(jbox) = 0
c       
c       Establish the offset of the box within the frame.
c       
        call plot(roff,soff,-3)
        call setlhe(-roff)
        call setbot(-soff)
c       
c       Save the last-used aspect ratio.
c       
        aspect0 = aspect
c       
	return
1001	return 1
c       
	end


	subroutine offbox(input,nin,ierbox,xbox,ybox,bscale,nbox,*)
        save
c
c	Choose a new output area with arbitrary (absolute) offset and size.
c
	character*(*) input
	dimension ierbox(0:4)
	dimension xbox(1),ybox(1),bscale(1)
c
        common /plot sizes/ xsize,ysize
        common /draw params/ roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &			     idevset,jbox,iorig
c
	common /initial/ roff0,soff0,hn0,hs0,hp0
c
        character*40 arg(4)
c
	if (nin.le.0) go to 99999
	iflag = 0
c
	call gettokens(input(1:nin),arg,narg)
        if (narg.lt.2) go to 99999
c
        call readrq(input(1:nin),narg,
     &              x,y,scale,xiflag,*99999)
c
        if (narg.lt.3) scale = 1.
        if (narg.lt.4) xiflag = 0.
c
        iflag = nint(xiflag)
c
c	Apply some checks (don't restrict y to be less than ysize, note).
c
	if (x.lt.0..or.x.gt.xsize.or.y.lt.0.
     &	    .or.scale.le.0..or.scale.gt.1.) go to 99999
	ioff = 0
	go to 1
c
	entry offbox1(x1,y1,scale1,iflag1,ierbox,*)
c
	x = x1
	y = y1
	scale = scale1
	iflag = iflag1
	ioff = 1
c
c	Undo any standard subdivisions.
c
1	if (jbox.ne.0) call newbox('=-0',3,ierbox,*99999)
c
c	Set the new origin (absolute units).
c
	call setorigin(x,y)
c
c	Set up the new box parameters.
c
        aspect = aspect1
        call setaspects(aspect,ax,ay)
c
c       Establish a fictitious "offset" for the box, even though (x,y)
c       as specified refers to the botom-left corner of the box itself.
c
c       (Note that setup has already been called, so roff0 and soff0 are
c        in "standard" units, rather than relative to the display size.)
c
	ss = sqrt(scale)
c
	roff = ax*ss*roff0
	soff = ay*ss*soff0
c
        call setlhe(-roff)
        call setbot(-soff)
c
	hn = ss*hn0
	hs = ss*hs0
	hp = ss*hp0
	call sethts(hs,hn)
c
c       Choose xlen and ylen, maintaining the correct aspect ratio,
c       and fitting into the space available.
c
        call fitbox(roff,soff,aspect,scale,xlen,ylen,xextra,yextra)
c
c       Clear the appropriate region of the screen, if requested.
c       Note that, in this case, the erased area follows the box,
c       as opposed to being a fixed quadrant of the screen.
c
	if (iflag.ne.0) call erase(-1.1*roff,1.1*xlen,-1.1*soff,1.1*ylen)
c
	if (ioff.eq.0) then
c
c	    Save everything on the stack.
c
	    nbox = nbox + 1
	    xbox(nbox) = x
	    ybox(nbox) = y
	    bscale(nbox) = scale
	end if
c
	return
99999	return 1
c
	end


	subroutine popbox(input,nin,ierbox,xbox,ybox,bscale,nbox,*)
        save
c
c	Pop a previous output area from the box stack.
c
	character*(*) input
	dimension xbox(1),ybox(1),bscale(1)
        character*40 arg(4)
c
	common /initial/ roff0,soff0,hn0,hs0,hp0
c
c       Count the number of arguments.
c
	call gettokens(input(1:nin),arg,narg)
c
        if (narg.gt.0) call readiq(input(1:nin),narg,
     &                             ipop,iflag,idum,idum,*99999)
        if (narg.lt.1) ipop = 1
        if (narg.lt.2) iflag = 0
c
	if (ipop.gt.nbox) then
	    write(6,*)'Box stack underflow'
	    go to 99999
	end if
c
	nbox = nbox - ipop
	if (nbox.gt.0) then
	    call offbox1(xbox(nbox),ybox(nbox),bscale(nbox),iflag,
     &		         ierbox,*99999)
	else
	    call offbox1(roff0,soff0,1.,iflag,ierbox,*99999)
	end if
c
        return
99999	return 1
c
	end


        subroutine setaspects(aspect,ax,ay)
c
c       Standard interpretation of the aspect ratio.
c       On return, aspect is the specified aspect ratio, or the
c       device default, ax is the reduction factor in the x-
c       direction, ay is the reduction factor in the y-direction.
c
        character*80 device
        common /plot device/ device,devaspect,idev
c
c       Use the device aspect ratio if none is specified.
c
        if (aspect.le.0.) aspect = devaspect
c
        if (aspect.le.devaspect) then
            ax = 1.
            ay = aspect
        else
            ax = devaspect/aspect
            ay = 1.
        end if
c
        end


        subroutine fitbox(roff,soff,aspect,scale,xlen,ylen,
     &                    xextra,yextra)
c
c       Ensure that the box has the correct aspect ratio and will
c       fit in the space available.  The "excess" space at the
c       right or top is returned in xextra or yextra.
c
c       Reduce roff and soff, if necessary, to fit the plot in.
c
        common /plot sizes/ xsize,ysize
        common /compress frames/ icompress
c
        parameter (WHITESPACE = 0.7, OFFTOL = 0.4)
c
c       Allow for offset of box, and (scaled) spacing at right.
c
        dx = scale*(xsize - WHITESPACE)
        if (roff.gt.OFFTOL*dx) roff = OFFTOL*dx
        xlen = dx - roff
c
c       Choose ylen to get the correct aspect ratio.
c
        ylen = xlen*aspect
c
c       Make sure that the y-axis will fit, by shrinking the entire
c       plot if necessary.
c
        dy = scale*(ysize - WHITESPACE*min(aspect,1.))
        if (soff.gt.OFFTOL*dy) soff = OFFTOL*dy
        ylmax = dy - soff
c
        xextra = 0.
        yextra = 0.
        if (icompress.ne.0) yextra = max(0., ylmax - ylen)
c
        if (ylen.gt.ylmax) then
            xlen0 = xlen
            xlen = xlen*ylmax/ylen
            ylen = ylmax
            if (icompress.ne.0) xextra = xlen0 - xlen
        end if
c
        end
