
c
c       Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
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

        subroutine ez2dplot(a,na,xlen,ylen,iprompt,itype,*)
        save
c	
c	Plot a 2-D contour representation of the array a.
c
	dimension a(1),xlevel(1000),temp(1000)
	character string*80,c*1,c2*1
	equivalence (string(1:1),c),(string(2:2),c2)
	character*50 list(20)
c
        common /local offset/ xoff,yoff
c
	save m,n,nl,xlevel,ipaint,invrt
	data m/0/n/0/nl/0/ipaint/0/invrt/0/
c
	a2(i,j) = a(i+m*(j-1))
c
	if (xlen.le.0..or.ylen.le.0.) return 1
c
	call devoff
	if (iprompt.eq.1) write(6,'(''*** Two-dimensional plot ***'')')
c
        call get_dimensions(m,n,na,iprompt,*99999)
c
	call minmax(a,m*n,amin,amax)
	if (iprompt.eq.1) write(6,'(a,1p,e12.4,a,e12.4)')
     &          'Array minimum = ',amin,', maximum = ',amax
	if (amin.ge.amax) return
c
	nl0 = nl
	ns = -1
c
	xoff = 0.
	yoff = 0.
c
c	Get the next sub-command (possibly prompting for input).
c
5	if (iprompt.eq.1) then
            call nextcmd(string,istart,ns,'Command(s) (h for help):')
        else
            call nextcmd(string,istart,ns,' ')
        end if
c
c	Parse the command string.
c
	if (c.eq.'q') then

	    if (nl.gt.0) then
		call sort(nl,xlevel)
		call devoff
		if (iprompt.eq.1) write(6,1000)(i,xlevel(i),i=1,nl)
1000		format('Levels: ',1p,4(i4,': ',e12.4):/
     &			(8x,4(i4,': ',e12.4)))
	    end if
	    return

	else if (c.eq.'h') then

	    call devoff
	    write(6,*)
	    write(6,*)'b: redraw box'
	    write(6,*)'co color: specify color'
	    write(6,*)'e: erase'
	    write(6,*)'d list: delete specified levels from the list'
	    write(6,*)'g: draw contour through graphics cursor'
	    write(6,*)'gd: delete contour through graphics cursor'
	    write(6,*)'l: list current contour levels'
	    write(6,*)'o: offset the contours'
	    write(6,*)'p: paint the box'
	    write(6,*)'p-: paint the box, inverted color map'
	    write(6,*)'q: back to main mcdraw prompt'
	    write(6,*)'r: repeat previously specified levels'
	    write(6,*)'s list: select specified levels from the list'
	    write(6,*)'t pattern: specify dashed lines'
	    write(6,*)'w weight: set line weight'
	    write(6,*)'x: delete contour list'
	    write(6,*)'<number>: draw contour at this level'
	    write(6,*)

	else if (c.eq.'b') then
c
	    call box(0.,xlen,0.,ylen)
c
	else if (c.eq.'c') then
c
c	    Set line color.
c
            call readiq(string(istart:ns),1,
     &                  icolor,idum,idum,idum,*5)
            call color(icolor)
c
	else if (c.eq.'d') then
c
c	    Delete specified levels, retain the rest.
c
	    if (nl.gt.0) then
	        call gettokens(string(istart:ns),list,nlist)
	        do 30 i=1,nlist
		    read(list(i),*,err=30,end=30)ii
		    call invert
		    call xcontor(a,m,n,xlevel(ii),1,0.,xlen,0.,ylen,0,
     &				 amin,amax,ipaint,invrt)
		    call invert
30	        continue
c
c		Remove the levels from the list.
c
		do 50 i=1,nlist
		    read(list(i),*,err=50,end=50)ii
		    do 40 j=ii+1,nl
                        xlevel(j-1) = xlevel(j)
40                  continue
		    if (nl.gt.0) nl = nl - 1
50		continue
	    end if

	else if (c.eq.'e') then
c
c	    Erase the screen.
c
	    call clear
	    call box(0.,xlen,0.,ylen)
	    ipaint = 0

	else if (c.eq.'g') then
	    call devoff
	    write(6,'(a)')'Entering graphics mode.'
     &		      //'  Exit by clicking outside box.'
	    call devon
c
c	    Graphic input.
c
60	    call getgfx(r,s)
	    if (r.le.0..or.s.le.0..or.r.ge.xlen.or.s.ge.ylen) go to 5
	    if (nl.eq.0.and.string(2:2).eq.'d') go to 5
c
	    x = 1. + r*(m-1.)/xlen
	    y = 1. + s*(n-1.)/ylen
	    ix = x
	    fx = x - ix
	    jy = y
	    fy = y - jy
	    xx = (1.-fy)*((1.-fx)*a2(ix,jy) + fx*a2(ix+1,jy))
     &		      + fy*((1.-fx)*a2(ix,jy+1) + fx*a2(ix+1,jy+1))
c
	    ii = 0
	    jtype = itype
	    if (string(2:2).eq.'d') then
c
c		Find the closest contour.
c
		dxmin = 1.e30
		jmin = 0
		do 65 j=1,nl
		    dx = abs(xx-xlevel(j))
		    if (dx.lt.dxmin) then
			dxmin = dx
			jmin = j
		    end if
65		continue
		if (jmin.gt.0) then
		    xx = xlevel(jmin)
		    do 66 j=jmin+1,nl
66		    xlevel(j-1) = xlevel(j)
		    if (nl.gt.0) nl = nl - 1
c
		    ii = -1
		    jtype = 0
		    call invert
		end if
	    else
		ii = iaddlevel(nl,xlevel,xx)
	    end if
c
c	    Always draw the contour (in case weight or color has changed).
c
	    call xcontor(a,m,n,xx,1,0.,xlen,0.,ylen,jtype,
     &                   amin,amax,ipaint,invrt)
	    if (ii.lt.0) call invert
	    go to 60
c
	else if (c.eq.'l') then
c
c	    Sort and list current contour levels.
c
	    call devoff
	    if (nl.gt.0) then
		call sort(nl,xlevel)
		write(6,1000)(i,xlevel(i),i=1,nl)
	    else
		write(6,'(a)')'No points.'
	    end if

	else if (c.eq.'o') then
c
c	    Set contour offset.
c
            call readrq(string(istart:ns),2,
     &                  xoff,yoff,dum,dum,*5)
c
	else if (c.eq.'p') then
c
c	    Paint the box, according to the values in a.
c
	    if (c2.ne.'-') then
		invrt = 0
	    else
		invrt = 1
	    end if
c
	    call color2d(xlen,ylen,a,m,n,amin,amax,invrt)
	    ipaint = 1
c
	else if (c.eq.'r') then
c
c	    Replay previously saved levels.
c
	    if (nl0.gt.0) then
		call xcontor(a,m,n,xlevel,nl0,0.,xlen,0.,ylen,itype,
     &				 amin,amax,ipaint,invrt)
		nl = nl0
	    end if
	else if (c.eq.'s') then
c
c	    Select specified levels, delete the rest.
c
	    if (nl.gt.0) then
	        call gettokens(string(istart:ns),list,nlist)
		nl1 = 0
	        do 70 i=1,nlist
		    read(list(i),*,err=70,end=70)ii
		    jj = iaddlevel(nl1,temp,xlevel(ii))
70	        continue
c
c		Redraw:
c
		call clear
		call box(0.,xlen,0.,ylen)
		nl = nl1
		do 80 i=1,nl
80		xlevel(i) = temp(i)
c
		call xcontor(a,m,n,xlevel,nl,0.,xlen,0.,ylen,itype,
     &				 amin,amax,ipaint,invrt)
	    end if
	else if (c.eq.'t') then
c
c	    Specify line type.
c
	    do 90 j=istart,ns
		if (string(j:j).gt.' ') go to 100
90	    continue
            itype=0
	    go to 5
100         call readiq(string(istart:ns),4,
     &                  i1,i2,i3,i4,*5)
            itype=1
            call setpat(i1,i2,i3,i4)
	else if (c.eq.'w') then
c
c	    Set line weight.
c
            call readiq(string(istart:ns),1,
     &                  iweight,idum,idum,idum,*5)
            call weight(iweight)
c
	else if (c.eq.'x') then
	    nl = 0
	else
	    call gettokens(string,list,nlist)
	    nl0 = nl
	    do 110 i=1,nlist
		if (list(i)(1:1).eq.'#') then
		    read(list(i)(2:len(list(i))),*,err=110,end=110)ix
		    if (ix.le.0.or.ix.gt.nl) go to 110
		    xx = xlevel(ix)
		else
                    call readrtoken(list(i),xx,xx)
		    ii = iaddlevel(nl,xlevel,xx)
		end if
c
c	        Always draw the contour (in case weight or color has changed).
c
	        call xcontor(a,m,n,xx,1,0.,xlen,0.,ylen,itype,
     &				 amin,amax,ipaint,invrt)
110	    continue
	end if
c
	nl0 = nl
	go to 5
c
9999	return
99999	return 1
	end


        subroutine get_dimensions(m,n,na,iprompt,*)
c
c       Read in and check array dimensions.
c
	character string*80
	character*50 list(20)
c
	mm = 0
	nn = 0
        if (iprompt.eq.1) then
            call getstring('Array sub-dimensions:',1,21,string)
        else
            call getstring(' ',1,1,string)
        end if
c
        call gettokens(string,list,nl)
        if (nl.eq.1) then
            call readitoken(list(1),mm,mm)
            if (mm.le.0) return 1
            nn = na/mm
            if (iprompt.eq.1) write(6,'(a,i6,a)')
     &              'Adopted ',nn,' as second array dimension.'
        else
            call readitoken(list(1),mm,mm)
            if (mm.le.0) return 1
            call readitoken(list(2),nn,nn)
            if (nn.le.0) return 1
        end if
c
	if (mm.le.0.or.nn.le.0) then
	    if (m.gt.0.and.n.gt.0) then
		if (iprompt.eq.1) write(6,'(a)')
     &                  'Using previous values '//
     &                  '-- may not replay correctly...'
	    else
		if (iprompt.eq.1) write(6,'(a)')'No previous values.'
		return
	    end if
	else
	    m = mm
	    n = nn
	end if
c
	if (m*n.gt.na) then
	    write(6,'(a)')'Too big!'
	    return 1
	end if
c
        end


	integer function iaddlevel(nl,xlevel,xx)
	save
c
c	Add the specified level to the contour list, if it isn't on it.
c
	dimension xlevel(1)
c
	iaddlevel = 0
	do i=1,nl
	    if (xx.eq.xlevel(i)) return
        end do
c
	iaddlevel = 1
	nl = nl + 1
	xlevel(nl) = xx
c
	end


	subroutine xcontor(a,m,n,v,nv,x1,x2,y1,y2,itype,
     &			   amin,amax,ipaint,invert)
        save
c
	dimension v(1)
        common /local offset/ xoff,yoff
c
c	Call contor or dcontor.
c
	if (ipaint.eq.0) then
	    if (itype.eq.0) then
		call contor(a,m,n,v,nv,x1+xoff,x2+xoff,y1+yoff,y2+yoff)
	    else
		call dcontor(a,m,n,v,nv,x1+xoff,x2+xoff,y1+yoff,y2+yoff)
	    end if
	else
	    call getcolor(icsave)
	    a2 = .5*(amin+amax)
c
	    do i=1,nv
c
		if ((v(i).le.a2.and.invert.eq.0)
     &			.or.(v(i).ge.a2.and.invert.eq.1))
     &		    call color(ncolors()-2)
c
		if (itype.eq.0) then
		    call contor(a,m,n,v(i),1,
     $                          x1+xoff,x2+xoff,y1+yoff,y2+yoff)
		else
		    call dcontor(a,m,n,v(i),1,
     $                           x1+xoff,x2+xoff,y1+yoff,y2+yoff)
		end if
	        call color(icsave)
            end do
	end if
c
	end


	subroutine color2d(xlen,ylen,array,nx,ny,amin0,amax0,invert)
        save
c
c	Color a rectangular region according to the values in array.
c
	real array(nx,ny)
	common /color limits/ icmin,icmax
c
	character*80 device
	real xp(4),yp(4)
c
	icmin = 0
	icmax = 0
c
c	Color Suns, X and PostScript only!
c
	call devicename(device)
	if (device(1:1).ne.'s'.and.device(1:1).ne.'p'
     $                        .and.device(1:1).ne.'x') return
	if (ncolors().le.2) return
c
c	Set up the ranges of the colors to be used.
c
	if (device(1:1).eq.'s') then
c
c	    Avoid strange colors at the ends of the SunCore color map!
c
	    icmin = 2
	    icmax = ncolors() - 2
	else if (device(1:1).eq.'x'.and.ncolors().gt.128) then
	    icmin = 1
	    icmax = ncolors() - 2
        else
	    icmin = 0
	    icmax = ncolors() - 1
	end if
	nc1 = icmax - icmin
c
	amin = amin0
	amax = amax0
	if (amin.ge.amax) call minmax(array,nx*ny,amin,amax)
	if (amin.ge.amax) return
c
c	Assume uniform spacing along both axes.
c
	dr = xlen/(nx-1.)
	ds = ylen/(ny-1.)
c
c	Color the box.
c
	do 100 i=1,nx-1
	do 100 j=1,ny-1
	    rr = (i-1)*dr
	    ss = (j-1)*ds
	    xp(1) = rr
	    yp(1) = ss
	    xp(2) = rr
	    yp(2) = ss + ds
	    xp(3) = rr + dr
	    yp(3) = ss + ds
	    xp(4) = rr + dr
	    yp(4) = ss
c
c	    Base color on average level in cell...
c
	    aa = .25*(array(i,j) + array(i+1,j)
     &				 + array(i,j+1) + array(i+1,j+1))
c
c	    Set up and clip the color, according to amin and amax:
c
	    ic = nint(icmin+nc1*(max(amin,aa)-amin)/(amax-amin))
	    ic = max(icmin,min(icmax,ic))
c
c	    Inverse color map:
c
	    if (invert.ne.0) ic = icmax + icmin - ic
c
	    call setfill(ic)
	    call polyfill(xp,yp,4)
c
100	continue
c
	end
