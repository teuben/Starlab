
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

        subroutine symbl(rinput,sinput,ht,string,theta,numch)
        save
c
c       Use the local string drawer to mimic "simbol."
c	(Not guaranteed to work in all circumstances!)
c
        character*80 device
        common /plot device/ device,aspect,idev
        common /plot sizes/ xsize,ysize
        common /plot origin/ ro,so
        common /fontc1/ offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
        common /numsym int/ rs,ss,dx,dy,sint,cost
        common /str posn/ i pos set,frx,fry
c
        common /plain font/ wid
c		( = ratio of character width to height.)
c
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /ncar/ nxpix1,nypix1,nx01,ny01,xfac1,yfac1
        common /dev init/ init
        common /dev status/ idevon,idevpen,idevwt
        common /fr int/ iframe
        common /numbr on/ inumbr

        character*200 temp
c
	character*40 font
	common /ps font/ font,ifsize
	save /ps font/
c
        external mcdxtext !$pragma C (mcdxtext)
c
        character*(*) string
        integer*2 ichars(100)
c
	r = rinput
	s = sinput
	go to 1
c
	entry csymbl(ht,string,theta,numch)
	call simwhe(r,s)
c
1       if (init.eq.0) then
            init = -1
            call mcinit
            call devon
            call clear
        end if
c
        if (idevon.eq.0)call devon
        nn = abs(numch)
        if (nn.eq.0)nn = len(string)
	if (nn.gt.len(string)) then
	    do 5 i = 1,len(string)-1
		if (string(i:i).eq.'%') then
		    if (string(i+1:i+1).eq.'%') then
			nn = i-1
			go to 6
		    end if
		end if
5	    continue
6	end if
        th = .017453293*theta
        sint = sin(th)
        cost = cos(th)
c
c	Get/estimate size of string.
c
        if (idev.eq.15) then
c
c	    SunCore calls!
c	    -------------
c
            call setcharsize(abs(ht*wid),abs(ht))
            call inqtextextent2(string(1:nn),dx,dy)
c
	    dy = abs(ht)
	    width = dx/nn
	    sint = 0.
	    cost = 1.
c
c	    (The SUN raster font used has fixed orientation.)
c
        else
            dy = abs(ht)
            width = wid*dy
            dx = width*nn
        end if
c
        if (idev.eq.2) then
c
c	    NCAR output:
c	    -----------
c
            rs = r+ro
            ss = s+so
            if (numch.lt.0) then
                n = 1
                write(ichars(1),'(2a1)')char(-numch),' '
            else
                n = numch
            end if
            isiz = yfac1*dy*wid
            icent = -1
            rs1 = rs
            ss1 = ss
            if (ht.lt.0..or.(iframe.eq.1.and.i pos set.eq.0)) then
                icent = 0
            else if (i pos set.ne.0) then
                dxs = 0.
                if (frx.ge.0..and.frx.le.1.)dxs = dx*frx
                dys = 0.
                if (fry.ge.0..and.fry.le.1.)dys = dy*fry
                rs1 = rs1-(dxs*cost-dys*sint)
                ss1 = ss1-(dxs*sint+dys*cost)
            end if
            rs = rs1
            ss = ss1
            if (icent.eq.-1) then
                rs1 = rs1-.5*dy*sint
                ss1 = ss1+.5*dy*cost
            end if
            ir = nx01+xfac1*rs1
            js = ny01+yfac1*ss1
            ith = theta
            if (numch.ge.0) then
                n = min(n,200)
                do 209 i = 1,n
209             if (string(i:i).ge.'a'.and.string(i:i).le.'z')
     +          string(i:i) = char(ichar(string(i:i))-32)
                i1 = -1
                do 210 i = 1,n/2
                    i1 = i1+1
                    i2 = i1+1
                    write(ichars(i),'(2a1)')string(i1:i1),string(i2:i2)
210             continue
                if (i2.lt.n)write(ichars(n/2+1),'(2a1)')string(n:n),' '
            end if
c
c	    NCAR call!!!
c	    ------------
c
            if (n.ne.0) call pwrit(ir,js,ichars,n,isiz,ith,icent)
c
            rs = rs-ro
            ss = ss-so
c
        else
c
c	    Attempt to determine the offset of the string.
c
            fx = 0.
            fy = 0.
            if (i pos set.ne.0) then
                if (frx.gt.0..and.frx.le.1.)fx = frx
                if (fry.gt.0..and.fry.le.1.)fy = fry
            else if (numch.lt.0.or.(inumbr.eq.1.and.iframe.eq.1)) then
                fx = .5
                fy = .5
            end if
c
            dxs = dx*fx
            dys = dy*fy
       	    rs = r-(dxs*cost-dys*sint)
            ss = s-(dxs*sint+dys*cost)
c
	    if (idev.eq.15) then
		rs = rs-.5*dy*sint
		ss = ss+.5*dy*cost
	    end if
c
c           (rs,ss) is the bottom left corner of the output string.
c
            if (idev.eq.5.or.idev.eq.6) then
c
c		HP plotter:
c		----------
c
                l = index(string(1:nn),'%%')
                if (l.gt.0)nn = l-1
                dx = width*nn
                ntilde = 0
                if (nn.le.0)go to 700
                call plot(rs,ss,3)
                fac = 2.5e-4*nxpix
c
c               (Each pixel is 1/40 mm, x-width is 10 units, "si"
c                below wants string dimensions specified in cm!)
c
                write(6,610).6666667*fac*width,fac*dy,cost,sint
610             format(' SI',f8.4,',',f8.4,';DI',f6.4,',',f6.4)
                islant = 1
                i1 = 1
650             if (i1.gt.nn)go to 700
                i2 = index(string(i1:nn),'~')+i1-1
                if (i2.lt.i1)i2 = nn+1
                i3 = i2-1
                if (i2.lt.nn.and.string(i2+1:i2+1).eq.'~')i3 = i2
                if (i3.ge.i1)call hp sym out(string(i1:i3),islant)
                ntilde = ntilde+1
                islant = 1-islant
                i1 = i3+2
                go to 650
700             dx = dx-width*max(0,ntilde-1)
c
            else if (idev.eq.15) then
c
c		SunCore:
c		-------
c
                call plot(rs,ss,3)
                call text(string(1:nn))
c
	    else if (idev.eq.16) then
c
c		Postscript:
c		----------
c
c		We can do better than the above guesses of string size
c		and offset.  Recall that r and s still represent the
c		input string coordinates, and fx and fy determine the offset.
c
		call plot(r,s,3)
		write(42,*)
c		write(42,*)'gsave %%%% Begin symbl output'
		write(42,*)'%%%% Begin symbl output'
c
c		The leading factor seems to be necessary to make a
c		(say) 72-point font come out 1 inch high...
c
		fheight = 0.645*ifsize
		scale = abs(ht)*yfac/fheight
c
c		Determine the length of the string -- stringlen:
c
		write(42,*)'/stringlen'
		write(42,800)'(',(string(i:i),i = 1,nn),')'
800		format(500a1)
		write(42,*)'stringwidth pop def'
c
c		Determine the r and s offsets (relative to string):
c
		write(42,810)fx,scale
810		format('/roffset stringlen ',f9.4,' mul ',
     &			f9.3,' mul def')
c
		write(42,820)fheight,fy,scale
820		format('/soffset ',2f9.4,' mul ',f9.3,' mul def')
c
c		Determine the x and y offsets (rotated):
c
		write(42,830)-cost,' r',sint,' s'
		write(42,830)-sint,' r',-cost,' s'
830		format(f9.5,a,'offset mul ',f9.5,a,'offset mul add')
c
c		Move to the correct offset point:
c
		write(42,*)'rmoveto'
c
c		Show the string:
c
		write(42,800)'(',(string(i:i),i = 1,nn),')'
c
c		Scale the string:
c
		write(42,'(2f9.3,'' scale'')')scale,scale
c
c		Rotate the string:
c
		write(42,'(f9.3,'' rotate'')')theta
c
c		Finish up:
c
c		write(42,*)'show grestore %%%% End symbl output'
		write(42,*)'show %%%% End symbl output'
		write(42,*)
c
            else if (idev.eq.17) then
c
c		X-windows:
c		---------
c
                call plot(rs,ss,3)
                do i=1,nn
                    temp(i:i) = string(i:i)
                end do
                temp(nn+1:nn+1) = char(0)
                call mcdxtext(abs(ht),theta,
     &                        temp(1:nn+1))
c
            end if
        end if
c
        call mc sym lims
c
        end
