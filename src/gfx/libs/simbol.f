
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
c         *    Fancy font-set symbol/number drawers.                  *
c         *                                                           *
c         *************************************************************
c
        subroutine simbol(xcall,ycall,hite,chars,theta,numch)
        save
c
c	extended-font version of symbol
c	positioning info is returned in common/fontc1/...
c	numch.lt.0  : string centered on (xi,yi)
c	numch.ge.0  : string starts with (xi,yi) at lower lh corner,
c		      or the point defined by subroutine strpos.
c	hite.le.0.  : only positioning info returned: nothing drawn
c
        character*(*) chars
	character*10 first10
        integer*2 n,m,num,jl,jr,idic,long
        integer*4 numch,nch
        real*4 lastinc
c
        common /fontc1/ offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
        common /sim fc2/ n,m,num(288),jl(288),jr(288),idic(288),
     &		         long(20000)
        common /sim ang/ sint,cost/sim len/nch,height
        common /str posn/ i pos set,frx,fry
        common /str limits/ offxmin,offxmax,offymin,offymax
        common /debug trace/ itrace
c
        data flag/0./
c
        xi = xcall
        yi = ycall
        go to 10
c
        entry usimbol(xcall,ycall,hite,chars,theta,numch)
        entry usersimbol(xcall,ycall,hite,chars,theta,numch)
        call fr inches(xcall,ycall,xi,yi)
c
10      if (itrace.eq.1) write(2,*)'simbol:',
     &                           xi,yi,hite,chars(1:1),numch
	first10 = chars(1:min(10,abs(numch)))
	call routine id('simbol '//first10)
c
	if (flag.eq.2.) return
c
        if (flag.eq.0.) then
	    n = 0
	    call getfonts
            flag = 1.
	    if (n.eq.0) then
	        flag = 2.
		return
	    end if
        end if
c
        height = 0.04762*abs(hite)
c       (factor is 1./21.)
c
        xmax = -1.e10
        ymax = xmax
        xmin = -xmax
        ymin = xmin
        sint = sin(theta*.0174533)
        cost = cos(theta*.0174533)
        nch = iabs(numch)
        if(nch.gt.1000)stop 'error: >1000 characters sent to simbol.'
        if(nch.eq.0)nch = 1000
c
        if(numch.lt.0.or.i pos set.ne.0)then
            call sim draw(xi,yi,chars,0,1)
            if(numch.lt.0)then
c
c               n.b. numch < 0 takes precedence over strpos.
c
                dx = .5*(xmax+xmin)-xi
                dy = .5*(ymax+ymin)-yi
            else
                if(frx.ge.0..and.frx.le.1.)then
                    dxs = (offxmin+(offxmax-offxmin)*frx)*height
                else
                    dxs = 0.
                end if
                if(fry.ge.0..and.fry.le.1.)then
                    dys = (offymin+(offymax-offymin)*fry)*height
                else
                    dys = 0.
                end if
                dx = dxs*cost-dys*sint
                dy = dxs*sint+dys*cost
            end if
            xmax = xmax-dx
            xmin = xmin-dx
            xp = xp-dx
            ymax = ymax-dy
            ymin = ymin-dy
            yp = yp-dy
            if(hite.le.0.)return
            call sim draw(xi-dx,yi-dy,chars,1,0)
        else
            if(hite.le.0.)call sim draw(xi,yi,chars,0,1)
            if(hite.gt.0.)call sim draw(xi,yi,chars,1,1)
        end if
c
        end


	subroutine getfonts
        save
c
c	Read the SIMBOL fonts.
c
        integer*2 n,m,num,jl,jr,idic,long
        common /sim fc2/ n,m,num(288),jl(288),jr(288),idic(288),
     &		         long(20000)
c
	parameter (NDIR = 10)
	character*120 directory(NDIR)
	dimension ldir(NDIR)
c
        parameter (IOPTION = 2)
c
c	Specify possible locations for the font file.
c
        nd = NDIR
	call listdir(directory,ldir,nd,iunit)
c
	if (iunit.lt.0) then
	    write(6,*)'No free unit number!'
	    return
	end if
c
c	Read the fonts.
c
20	do 50 idir=1,nd
            n = ldir(idir)
	    if (n.le.0) go to 50
c
c-----------------------------------------------------------------------
c
            if (IOPTION.eq.1) then
c
c               Old version:
c               -----------
c
                open(iunit,status='old',form='unformatted',
     &               file=directory(idir)(1:n)//'/SIM.UNF',err=50)
                read(iunit)n,m,num,jl,jr,idic,(long(i),i=1,m)
c
            else if (IOPTION.eq.2) then
c
c               New version:
c               -----------
c
c                write(6,*)idir,'. Trying font file ',
c     &                    directory(idir)(1:n)//'/SIM.out'
c
                open(iunit,status='old',form='formatted',
     &               file=directory(idir)(1:n)//'/SIM.out',err=50)

                read(iunit,*)n,m
c       
                j = 0
                l = 1
100             read(iunit,*,err=101,end=101)
                read(iunit,*,err=101,end=101)
                read(iunit,*,err=101,end=101)i1,i2,i3
                j = j + 1
                num(j) = i1
                jl(j) = i2
                jr(j) = i3
c
                if (num(j).gt.0) then
                    idic(j) = l
                    read(iunit,*)(long(k),k=l,l+num(j)-1)
                    l = l + num(j)
                else
                    idic(j) = -1
                end if
c               
                go to 100
c
            else
                go to 50
            end if
c
c-----------------------------------------------------------------------
c
101	    close(iunit)
	    return
c
50	continue
c
	write(6,'(a)')'Can''t find the fonts!'
c
	end


	subroutine listdir(directory,ldir,ndir,iunit)
        save
c
c	Make a list of places to look for useful files, and return
c	a free unit number.  (Also used by fullhelp.)
c
c       On entry, ndir specifies the maximum allowable number of
c       list entries.  On exit, it is the actual length of the list.
c
	character*120 directory(ndir)
	dimension ldir(ndir)
	logical opened
c
	if (ndir.le.0) return
c
        do idir=1,min(6,ndir)
            directory(idir) = ' '
        end do
c
c	UNIX calls!
c	----------
c
        nd = 1
	call mygetenv('MCD_FONT_DIR',directory(nd))
c
c       Now list various likely locations.  Note the order of precedence.
c
	if (nd.lt.ndir) then
            nd = nd + 1
            directory(nd) = '.'
        end if
c
	if (nd.lt.ndir) then
            nd = nd + 1
            directory(nd) = '..'
        end if
c
	if (nd.lt.ndir) then
            nd = nd + 1
            nd = nd + 1
            call mygetenv('HOME',directory(nd))
            do n=120,1,-1
                if (directory(nd)(n:n).gt.' ') then
                    directory(nd)(n+1:n+4) = '/bin'
                    go to 20
                end if
            end do
20          continue
        end if
c
	if (nd.lt.ndir) then
            nd = nd + 1
            directory(nd) = '/usr/local/mcdraw'
        end if
c
	if (nd.lt.ndir) then
            nd = nd + 1
            directory(nd) = '/usr/local/mcdraw/libs'
        end if
c
	if (nd.lt.ndir) then
            nd = nd + 1
            directory(nd) = '/usr/local/mcdraw/draw'
        end if
c
        if (nd.lt.ndir) then
            nd = nd + 1
            call mygetenv('STARLAB_PATH',directory(nd))
            do n=120,1,-1
                if (directory(nd)(n:n).gt.' ') then

                    directory(nd)(n+1:n+13) = '/src/gfx/libs'
c
c     		    Add this for help documentation...
c
                    nd = nd + 1
                    directory(nd) = directory(nd-1)
                    directory(nd)(n+1:n+13) = '/src/gfx/draw'
                    go to 30
                end if
            end do
30          continue
        end if
c
        if (nd.lt.ndir) then
            nd = nd + 1
            call mygetenv('SIM_FONT_DIR',directory(nd))
        end if
c
        ndir = nd
c
c	Determine lengths of directory names.
c
	do idir=1,nd
	    do n=120,1,-1
		if (directory(idir)(n:n).gt.' ') then
		    ldir(idir) = n
		    go to 50
		end if
            end do
            ldir(idir) = 0
50	end do

c        do ii=1,nd
c           write(6,*)ii,directory(ii)(1:ldir(ii))
c        end do

c
c	Find a free unit number.
c
	do 100 iunit=10,99
	    inquire(iunit,opened=opened)
	    if (.not.opened) go to 120
100	continue
	iunit = -1
c
120	end


	block data sboxsetup
        save
	common /sboxdata/ iborder,ierase,fraction
	data iborder/0/ierase/0/fraction/.5/
	end


	subroutine sboxset(ib,ie,fr)
        save
	common /sboxdata/ iborder,ierase,fraction
c
	iborder = ib
	ierase = ie
	if (fr.ge.0.) fraction = fr
c
	end


	subroutine boxsim(x,y,ht,string,theta,nch)
        save
c
c	Same as "simbol," but optionally first erase the background
c	and/or add a bounding box.
c
	character*(*) string
c
	common /sboxdata/ iborder,ierase,fraction
        common/fontc1/offx,offy,lastinc,xp,yp,xmax,xmin,ymax,ymin
        dimension xbox(4),ybox(4)
c
	if (iborder.ne.0.or.ierase.ne.0) then
c
c	    Get string dimensions and location.
c
	    call simsize(ht,string,nch,dx,dy)
	    call simbol(x,y,-ht,string,theta,nch)
	    call simwhe(xp,yp)
c
c	    String "end" (xp) includes some space.  Correct by adjusting dx.
c
	    dx = dx + .5*ht
c
c	    Set up box containing the string.
c
	    cost = cos(3.14159*theta/180.)
	    sint = sin(3.14159*theta/180.)
c
	    xbox(1) = xp - dx*cost
	    xbox(2) = xp
	    xbox(3) = xp - dy*sint
	    xbox(4) = xp - dx*cost - dy*sint
	    ybox(1) = yp - dx*sint
	    ybox(2) = yp
	    ybox(3) = yp + dy*cost
	    ybox(4) = yp - dx*sint + dy*cost
c
c	    Enlarge the box (anisotropically).
c
c	    "Fraction" refers to the height direction.  Along the string,
c	    cap the enlargement at ht.
c
	    ddy = .5*fraction*dy
	    ddx = max(ddy,min(.5*ht,.5*fraction*dx))
c
	    xbox(1) = xbox(1) - ddx*cost + ddy*sint
	    xbox(2) = xbox(2) + ddx*cost + ddy*sint
	    xbox(3) = xbox(3) + ddx*cost - ddy*sint
	    xbox(4) = xbox(4) - ddx*cost - ddy*sint
c
	    ybox(1) = ybox(1) - ddx*sint - ddy*cost
	    ybox(2) = ybox(2) + ddx*sint - ddy*cost
	    ybox(3) = ybox(3) + ddx*sint + ddy*cost
	    ybox(4) = ybox(4) - ddx*sint + ddy*cost
c
c	    Add the embellishment(s).
c
	    if (ierase.ne.0) call polyerase(xbox,ybox,4)
	    if (iborder.ne.0) call polydraw(xbox,ybox,4)
c
	end if
c
c	Draw the string, as usual.
c
	call simbol(x,y,ht,string,theta,nch)
c
	end
