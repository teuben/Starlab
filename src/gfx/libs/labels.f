
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
        
        subroutine labels(xctit,yctit)
        save
c       
        character*(*) xctit,yctit
c       
        common /scales/ xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
        common /fr hts/htl,htn
        common /fontc1/ dum2(5),rmax,rmin,smax,smin
        common /fr xnums/ xnumbot
        common /fr ylab pos/ slab
        common /fr rotn/ irot
        common /fr setax/ kax,lax
        common /dev details/ itek,ivers
        common /debug trace/ itrace
        common /fr conf/ scent,rnuml,rnumr,snumt,snumb,
     &                     dsnums,jrot,stopnum
c       
        parameter (kountmax = 1)
c       
        logical debug
        data debug/.false./
c       
c       Label axes with extended font set.  S/R fr spaces now looks for 
c       five spaces or a non-printing character to terminate the string.
c       
c       As specified, label sizes are both determined by htl.  Allow the
c       labels to be rescaled by embedding a '!scale!' at the start.
c
c       Note that this routine makes no permanent changes to pen/font size,
c       etc.  It is assumed that these have already been set up in advance.
c
c       First, the x-label.
c       ------------------
c       
        call fr spaces(xctit,nxtit)
c       
        if (nxtit.gt.0.and.kax.ne.0) then
c           
c           Plot the x-label
c           
            lx = len(xctit)
c           
c           Check for embedded rescaling.
c           
            call lscale(xctit,nxtit,factor,ifirst)
            htlsave = htl
            htl = htl*factor
            nxtit = nxtit - ifirst + 1
c           
c           Determine the extent of the label.
c           
            call simbol(0.,0.,-htl,xctit(ifirst:lx),0.,nxtit)
            siml = rmax-rmin
            simh = smax-smin
c           
c           Dimensions of the label are siml wide by simh high.
c           
            if (siml.gt.rlen) then
                htl = htl*rlen/siml
                simh = simh*rlen/siml
            end if
            xltop = xnumbot-.3*(htn+htl)
            simbot = xltop-sihm
c           
c           Draw the label.
c           
            call getbot(spbot)
            if (simbot.lt.spbot)htl = htl*(xltop-spbot)/simh
            call strpos(.5,1.)
            if (itrace.eq.1) write(2,*)'x call to simbol:',
     &              .5*rlen,xltop,htl,xctit(ifirst:ifirst+nxtit-1)
            call simbol(.5*rlen,xltop,htl,xctit(ifirst:lx),0.,nxtit)
            call clrstr
c           
            htl = htlsave
        end if
c       
c       Now for the y-label.
c       -------------------
c       
        if (lax.lt.0) return
c       
        if (debug) write(6,*)'y-label...'
        kount = 0
        ileadsp = 0
        inbl = 1
        ly = len(yctit)
        htlsave = htl
c       
c       Convention:  leading spaces in yctit
c        =  = > strip and plot the label VERTICALLY.
c       
        call fr spaces(yctit,nytit)
        do 100 inbl = 1,nytit
            if (yctit(inbl:inbl).gt.' ') go to 110
100     continue
        go to 500
c       
110     if (inbl.gt.1) then
            ileadsp = 1
            nytit = nytit-inbl+1
        end if
c       
c       Check for embedded rescaling.
c       
        call lscale(yctit(inbl:inbl+nytit-1),nytit,factor,ifirst)
        htl = htl*factor
        inbl = inbl + ifirst - 1
        nytit = nytit - ifirst + 1
c       
c       Determine the extent of the label.
c       
        call simbol(0.,0.,-htl,yctit(inbl:ly),0.,nytit)
        siml = rmax - rmin
        simh = smax - smin
c       
c       Dimensions of the label are siml wide by simh high.
c       Now determine its left and right edges.  Note that this uses
c       information returned by the last calls to the number-drawing
c       routines used by eframe:
c       
c       rnuml = the left-most extent of the numeric labels
c       rnumr = the right-most extent of the numeric labels
c       scent = the vertical center of the box
c       snumb = the level of the next numeric label below the center
c       snumt = the level of the next numeric label above the center
c       dsnums = the vertical label spacing
c       stopnum = the level of the top numeric label
c       jrot = 1 iff the numeric labels are drawn vertically
c       
        if (debug) write(6,*)'lax,ileadsp,siml,irot,nytit = ',
     &          lax,ileadsp,siml,irot,nytit
        
        if (lax.le.1) then
            rttlmax = rnuml-.5*htn
            call getlhe(rttlmin)
            rttlmin = max(rttlmin,-.3*rlen)
        else
            rttlmin = rnumr+.5*htn
            call getrhe(rttlmax)
            rttlmax = min(rttlmax,1.3*rlen)
        end if
c       
c       The basic dimensions of the label layout are now set.
c       
c       See if the user wants a horizontal label.
c       
        if (ileadsp.eq.0.and.irot.eq.0
     &          .and.(siml.le.7.5*simh.or.siml.le..25*slen)) go to 300
c       
c       Plot the y-label vertically
c       
200     slab = scent
        if (lax.le.1) then
            rlab = rttlmax-simh
            rlablhs = rttlmax-2.*simh
            rlabrhs = rlab
        else
            rlab = rttlmin+simh
            rlablhs = rlab
            rlabrhs = rttlmin+2.*simh
        end if
c       
c       Check label positioning.
c       
        if (debug) write(6,*)'Checking position for vertical label...'
        if (rlablhs.lt.rttlmin.or.rlabrhs.gt.rttlmax) then
            drwant = rttlmax-rttlmin
            if (drwant.le..25*simh.and.(itek.ne.1.or.ivers.eq.1)) then
                call display text('No room for the y-label.',24)
                go to 500
            end if
            if (debug) write(6,*)'simh,drwant,htn  = ',simh,drwant,htn
            if (simh-drwant.lt..25*htn) then
 	        if (lax.le.1) then
                    rlab = rttlmin+.5*simh
	        else
                    rlab = rttlmax-.5*simh
	        end if
            else
                htl = htl*drwant/simh
                siml = siml*drwant/simh
                rlab = .5*(rttlmin+rttlmax)
            end if
        end if
c       
        if (siml.gt.slen)htl = htl*slen/siml
        th = 90.
        if (lax.gt.1) th = -th
        if (debug) write(6,*)'th = ',th
c       
        go to 400
c       
c       Try to make the y-label horizontal.
c       ----------------------------------
c       
c       First, check vertical positioning
c       
300     slab = scent
        imovr = 0
        if (debug) write(6,*)'At 300: jrot = ',jrot
c       
        if (jrot.eq.0) then
c           
c           Horizontal numerical labels.
c           
            if (min(snumt-scent,scent-snumb).gt.1.5*(simh+htn))
     &              go to 350
c           
            slab = .5*(snumb+snumt)
            if (slab.lt.scent.and.slab+dsnums.lt.stopnum) then
                snumb = snumb+dsnums
                snumt = snumt+dsnums
                slab = slab+dsnums
            end if
            sclear = min(snumt-slab-.5*simh,slab-snumb-.5*simh)
            if (sclear.gt..75*simh) go to 350
        end if
c       
c       Not enough clearance to move right(left): put label's r(l)h edge
c       at l(r)h edge of numbers
c       
        rright = rttlmax
        rleft = rttlmin
        go to 360
c       
c       Label can be moved right(left).
c       
350     if (debug) write(6,*)'At 350: Move right...'
        if (lax.le.1) then
            rright = rttlmax
            if (rnuml.lt.-2.*htl) rright = .75*rnuml
            rleft = rright-siml
            if (rleft.gt.rnuml)rright = rnuml+siml
        else
            dr = rnumr-slen
            rleft = rttlmin
            if (dr.gt.2.*htl) rleft = rnumr-.25*dr
            rright = rleft + siml
            if (rright.lt.rnumr)rleft = rnumr-siml
        end if
c       
        imovr = 1
c       
c       Right(left)-hand edge of label is now set.  Check the left(right) end.
c       
360     if (lax.le.1) then
            room = rright-rttlmin
        else
            room = rttlmax-rleft
        end if
        if (debug) write(6,*)'room = ',room
c       
        if (room.lt.siml) then
            kount = kount+1
            if (kount.le.kountmax) then
                siml = siml*.95
                simh = simh*.95
                htl = htl*.95
                if (debug) write(6,*)'Reducing y-label size.'
                go to 300
            end if
            if (kount.eq.kountmax+1.and.(itek.ne.1.or.ivers.eq.1))
     &              call display text('Warning: problem plotting'//
     &              ' y-label horizontally.',47)
            if (room.lt..9*siml) then
                if (imovr.eq.1) then
c                   
c                   Try to move label farther right(left).
c                   
                    imovr = 2
                    if (lax.le.1) then
                        rright = min(rttlmin + siml,-htn)
                        rleft = rright - siml
                    else
		        rleft = rttlmax
		        rright = rleft + siml
                    end if
                    go to 360
                end if
c               
c               Give up--plot the label vertically!
c               
                siml = siml*htlsave*factor/htl
                simh = simh*htlsave*factor/htl
                htl = htlsave*factor
                if (itek.ne.1.or.ivers.eq.1)
     &                  call display text('Attempting to plot the'//
     &                  ' y-label vertically.',42)
                go to 200
            end if
            kount = kount+1
            fac = room/siml
            htl = htl*fac
            siml = siml*fac
            simh = simh*fac
            go to 350
        end if
        rlab = rright-.5*siml
        th = 0.
c       
c       Plot the y-label.
c       
400     if (debug) write(6,*)'Plotting the label !',rlab,slab
        call simbol(rlab,slab,htl,yctit(inbl:ly),th,-nytit)
        if (itrace.eq.1) write(2,*)'y call to simbol:',
     &                         rlab,slab,htl,yctit(1:nytit)
c
500     htl = htlsave
c
        end
