
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
	subroutine readall(lunit,string,i1,i2,x,n,nhead,nmax,
	1	           iprompt,*)
        save
c
c	Read the array x from the specified logical unit, without worrying
c	about columns of data.
c
	character*(*) string
	character*400 line
	character*50 list(100)
	dimension x(1)
	character*1 dummy
c
	if (i1.le.i2) then
	    nin = 0
            call readiq(string(i1:i2),1,
     &                  nin,idum,idum,idum,*10)
	    nin = min(nin,nmax)
	else
	    nin = 0
	end if
c
c	Skip a header, if necessary.
c
10	rewind lunit
	do 20 i=1,nhead
20	read(lunit,'(a)')dummy
c
	if (nin.gt.0) then
	    read(lunit,*,err=999,end=999)(x(i),i=1,min(nmax,nin))
	    if (nin.gt.nmax.and.iprompt.ne.0) write(6,'(a)')
	1	    'Maximum number of points reached.'
	    n = nin
	else
c
c           Note: specifying no number reads data the "smart but slow" way.
c
	    n = 0
30	    read(lunit,'(a)',err=50,end=50)line
	    call gettokens(line,list,nlist)
	    do i=1,nlist
		read(list(i),*,err=50,end=50)xx
		if (n.ge.nmax) then
		    if (iprompt.ne.0) write(6,'(a)')
	1		    'Maximum number of points reached.'
		    go to 50
		end if
		n = n + 1
		x(n) = xx
	    end do
	    go to 30
	end if
c
50	if (iprompt.ne.0) write(6,*)n,' points read'
	return
c
999	if (iprompt.ne.0) write(6,*)'Error reading input file'
	return 1
	end
