
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

        subroutine sim draw(xlhe,ylhe,chars,letts,bounds)
        save
        character*(*) chars
        integer*2 n,m,num,jl,jr,idic,long
        common/sim fc2/n,m,num(288),jl(288),jr(288),idic(288),
     +                 long(15610)
c
c	Font storage:
c
c	Each of the 288 characters is stored as a series of points (i,j) on
c	an integer grid. i=0 is at the (horizontal) center of the character
c	(specifically, it is on the axis of symmetry of an "a", "h", etc.).
c	j=0 is near the vertical center of a capital letter (3/7 of the way
c	up a capital "a", in fact). Capital letters are contained within
c	a 25 x 25 grid (i,j ranging from -12 to 12) and have their bases at
c	j=-9. In this routine, (xlhe,ylhe) refers to the bottom left-hand
c	corner of the string, so an offset of 9 is applied to all j-values.
c	For a given character j, jl(j) and jr(j) are the left- and right-
c	hand edges of the figure, possibly including some space, and define
c	horizontal width actually plotted. Hence, jl is also used as an
c	offset below. given the baseline, "height" is set up in "simbol"
c	so that a capital "a" is of the size input to that routine. Since
c	a capital "a" ranges from j=-9 to j=12, the "height" here is the
c	"hite" input to simbol divided by 21 (see implementation below).
c	num(j) is twice the number of points required to specify the
c	character (2 coordinates per point). Points are stored in long in
c	the form i1,j1,i2,j2,i3,j3, etc. idic(j) points to the i
c	coordinate of the first point of #j. Note that no coordinate
c	lies outside the range (-20,20), so long could be stored as a byte
c	array, if desired.
c
c	n, m are the lengths of the arrays: 288 and 15610, respectively.
c
        integer bounds,backs
        integer*4 nch
        common/sim ang/sint,cost/sim len/nch,height
        real*4 lastinc
        common/fontc1/offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
c
c	On return, offx and offy are the coordinates, in internal units,
c	of the bottom right-hand corner of the string (note: the actual
c	plotted angle is irrelevant). lastinc is the width of the last
c	character, in the same units. xp and yp are the "real" coordinates
c	(i.e. "inches") of the (offx,offy) point, i.e. (xp,yp) is where to
c	start if it is desired to continue with simbol. The entire string
c	is contained in the box defined by (xmin, xmax, ymin, ymax).
c	Note that the box is defined by the figure itself, not by the
c	box used in its construction.
c
        common/str limits/offxmin,offxmax,offymin,offymax
        common/framesize/nxpix,nx0,xfac,nypix,ny0,yfac
        common/debug trace/itrace
c
        character*80 device
        common /plot device/device,aspect,idev
c
        if(itrace.eq.1)write(2,*)'sim_draw:',
     +                           xlhe,ylhe,chars,letts,bounds

        joff=0
        backs=0
        talic=0.                                                              ! initially, italics are off.
        offx=0.
        offy=0.
        offxmin=1.e10
        offxmax=-1.e10
        offymin=1.e10
        offymax=-1.e10
        subs=1.
        hc=height*cost
        hs=height*sint
c
c	***** loop over characters to be plotted *****
c
        do 10 j10=1,nch
            kcur=ichar(chars(j10:j10))
            if(kcur.eq.64)then
                if(joff.eq.192)then
                    subs=subs*0.6                                                        ! %@ = decrease size
                    joff=-96
                end if
                joff=joff+96                                                             ! @ = second font
                if(joff-192)10,1,11                                                      ! @@ = backspace
1               backs=1
                joff=0
            else if (kcur.eq.37) then
                if(joff.eq.96)then
                    subs=subs/0.6                                                        ! @% = increase size
                    joff=-192
                end if
                joff=joff+192                                                            ! % = third font
                if(joff.gt.192)go to 11                                                  ! %% = end
            else if (kcur.eq.94) then
                subs=subs*0.6                                                            ! ^ = superscript
                offy=offy+21.*subs
            else if (kcur.eq.92) then
                subs=subs*0.6                                                            ! \ = subscript
                offy=offy-12.*subs
            else if (kcur.eq.126)then
                talic=1.-talic                                                           ! ~ = italics off/on
            else
                j=kcur-31+joff
                lefte=jl(j)                                                              ! offset of letter within box
                if(backs.eq.1.and.nback.eq.0)oldinc=lastinc
                lastinc=(jr(j)-jl(j))*subs                                               ! width of (letter+rh edge)
                if(backs.eq.1)then
                    diff=.5*(oldinc-lastinc)
                    offx=offx-oldinc+diff
                    lastinc=lastinc+diff
                end if
                np=num(j)/2                                                              ! # of pen strokes
                lastdrawn=0
                npp=np
                if(np.eq.0)then
                    npp=1
                    go to 10005
                end if
                lastdrawn=1
                indx=idic(j)                                                             ! where to start
                kc=j+31
c
c                tal=0.
c                if(kc.ge.65.and.kc.le.90)tal=talic                                        ! only italicize letters
c                if(kc.ge.97.and.kc.le.122)tal=talic
c
                 tal=talic                                                                ! italicize everything
c
c               Force reference point to be at a standard position
c               relative to device pixels.
c
c               This doesn't seem to work with the versaterm option (?).
c
                if(abs(cost).gt..7071)then
                    xref=xlhe+offx*hc
                    dxp=nx0+xref*xfac
                    ixp=dxp
                    dxp=dxp-ixp
                    if(dxp.ne.0..and.idev.ne.16)then
                        if(dxp.le..25)then
                            do=-dxp
                        else
                            do=1.-dxp
                        end if
                        offx=offx+do/(xfac*hc)
                    end if
                else
                    yref=ylhe+offx*hs
                    dyp=yref*yfac
                    iyp=dyp
                    dyp=dyp-iyp
                    if(dyp.ne.0..and.idev.ne.16)then
                        if(dyp.le..25)then
                            do=-dyp
                        else
                            do=1.-dyp
                        end if
                        offx=offx+do/(yfac*hs)
                    end if
                end if
c
10005           up=1.
c
c	        ***** loop over pen strokes *****
c
                do 5 l=1,npp
                    if(np.eq.0)then
                        ix=0
                        iy=0
                        subs=0.
                        go to 10105
                    end if
                    ix=long(indx)
                    iy=long(indx+1)
                    indx=indx+2
10105               if(ix.eq.31)then
                        up=1.
                        go to 5
                    else                                                                 !   off is the current offset
                        xl=offx+subs*(ix-lefte+0.20*tal*(iy+9))                          !  relative to the start of the
                        yl=offy+subs*(iy+9)                                              ! string (i.e. the bottom lh edge
                        if(bounds.eq.1)then                                              !  of a leading capital letter.)
                            if(xl.gt.offxmax)offxmax=xl
                            if(xl.lt.offxmin)offxmin=xl
                            if(yl.gt.offymax)offymax=yl
                            if(yl.lt.offymin)offymin=yl
                        end if
                        xp=xlhe+height*(xl*cost-yl*sint)
                        yp=ylhe+height*(yl*cost+xl*sint)
                    endif
                    kount bounds=0
                    if(up.eq.1.)then
                        if(letts.eq.1)call plot(xp,yp,3)
                        up=0.
                        xpo=xp
                        ypo=yp
                    else
                        if(letts.eq.1)call plot(xp,yp,2)
                        kount bounds=1
                    end if
                    if(bounds.eq.1.and.kount bounds.eq.1)then
                        xmax=max(xmax,xp,xpo)
                        xmin=min(xmin,xp,xpo)
                        ymax=max(ymax,yp,ypo)
                        ymin=min(ymin,yp,ypo)
                    end if
5               continue
6               joff=0
                nback=0
                if(backs.eq.1)nback=1
                backs=0
                subs=1.
                offy=0.
                offx=offx+lastinc
            endif
10      continue
11      continue
        xp=xlhe+height*(offx*cost)
        yp=ylhe+height*(offx*sint)
        if(lastdrawn.eq.0)then
c
c	    deal with trailing blanks.
c
            call plot(xp,yp,3)
            if(bounds.eq.1)then
                if(xp.gt.xmax)xmax=xp
                if(yp.gt.ymax)ymax=yp
                if(xp.lt.xmin)xmin=xp
                if(yp.lt.ymin)ymin=yp
            end if
        end if
c
        end
