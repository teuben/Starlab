
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

	subroutine readcols(lunit,nhead,nread,ix,iy,x,y,z,n,nmax,
     &			    c1,info,iprompt,*)
        save
c
c	Read data from the specified column(s) of a file.  Skip the
c	first nhead lines, and read at most min(nread,nmax) lines.
c
	dimension x(nmax),y(nmax),z(nmax)
	character*1 c1
c
c       Read data in double precision, to allow offset if necessary.
c
	real*8 dummy(1000)
c
        character*4000 line
        character*40 column(200)
c
	common /data offset/ delx,dely,delz,facx,facy,facz
c
        common /file input/ inpmode
c
        character*1 comment
        common /file comment/ icomment,comment
c
        data icomment/0/comment/' '/
c
        n=0
        i1 = min(ix,iy)
        i2 = max(ix,iy)
c
c	Skip the header lines.
c
        rewind lunit
        do 10 i=1,nhead
10      read(lunit,*,err=120,end=120)
c
110     if (n.ge.nmax) then
            if (info.eq.1) then
                call devoff
                if (iprompt.eq.1) write(6,101)n
101             format('Warning: maximum # of points (',i6,
     &                 ') reached.')
            end if
            n = n + 1
        else
            n = n + 1
c
c           Note that the first line is read the hard way to check
c           that the number of columns is OK.
c
            if (n.gt.1.and.inpmode.eq.0) then
                if (i2.gt.i1) then
                    read(lunit,*,err=120,end=120)(dummy(i),i=1,i2)
		    x(n) = facx * (dummy(i1) - delx)
		    y(n) = facy * (dummy(i2) - dely)
                else
		    read(lunit,*,err=120,end=120)(dummy(i),i=1,i1)
                    if (c1.eq.'x') then
			x(n) = facx * (dummy(i1) - delx)
                    else if (c1.eq.'y') then
			y(n) = facy * (dummy(i1) - dely)
                    else if (c1.eq.'z') then
			z(n) = facz * (dummy(i1) - delz)
                    else
                        return 1
                    end if
                end if
            else
c
c               (MUCH) slower, but surer...
c
                read(lunit,'(a)',err=120,end=120)line
c
c               Locate the end of the string ('                    ').
c               This will fail if the intercolumn spacing is too great...
c
                nl = index(line,'                    ') - 1
c
                if (nl.le.0) then
                    n = n - 1
                    go to 118
                end if
c
c               Check for comments.
c
                if (icomment.ne.0.and.line(1:1).eq.comment) then
                    n = n - 1
                    go to 118
                end if
c
c               Split line into columns.
c
                call gettokens(line(1:nl),column,nc)
                if (nc.lt.max(i1,i2)) then
                    if (iprompt.eq.1) write(6,'(a,i3,a)')
     &                      'Data file has only ',nc,' columns'
                    return 1
                end if
c
                if (i2.gt.i1) then
                    call readdtoken(column(i1),dummy,dummy)
		    x(n) = dummy(1)
                    call readdtoken(column(i2),dummy,dummy)
		    y(n) = dummy(1)
		    x(n) = facx * (x(n) - delx)
		    y(n) = facy * (y(n) - dely)
                else
                    if (c1.eq.'x') then
                        call readrtoken(column(i1),x(n),x(n))
			x(n) = facx * (x(n) - delx)
		    else if (c1.eq.'y') then
                        call readrtoken(column(i1),y(n),y(n))
			y(n) = facy * (y(n) - dely)
                    else if (c1.eq.'z') then
                        call readrtoken(column(i1),z(n),z(n))
			z(n) = facz * (z(n) - delz)
                    else
                        return 1
                    end if
                end if
c
                if (iprompt.eq.1.and.100*(n/100).eq.n)
     &                  write(6,*)n,' lines read...'
c
            end if
c
118         if (n.lt.nread) go to 110
	    n = n + 1
c
	end if
c
120     n = n - 1
c
        if (info.eq.1) then
            call devoff
            if (iprompt.eq.1) write(6,125)n
125         format(i6,' points')
        end if
c
        if (ix.gt.iy.and.n.gt.0) call swap(x,y,n)
c
	end
