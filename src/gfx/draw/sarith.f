
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

	subroutine sarith(input,istart,nin,arr,nmax,narr,register,
     &		 	  iprompt,*)
        save
c
c	Perform scalar arithmetic (on arrays).
c
	character*(*) input
	dimension arr(nmax,3),narr(3)
c
	character*1 c2
	dimension iarg(3)
c
	if (nin.le.istart) go to 1001
	c2 = input(2:2)
c
        ibl=1
        narg=0
        do 450 i=istart,nin
            if (input(i:i).eq.' '.or.input(i:i).eq.',') then
                ibl=1
            else
                if (ibl.eq.1) then
                    narg=narg+1
                    if (narg.ne.2) then
                        call decode(input(i:i),iarg(narg),*1001)
                    else
                        read(input(i:nin),*,err=1001,end=1001)s
                    end if
                    if (narg.eq.3) go to 451
                    ibl=0
                end if
            end if
450     continue
c
        if (narg.lt.1) then
            go to 1001
	else if (narg.eq.1) then
	    s = register
        end if
c
        if (narg.le.2) iarg(3)=iarg(1)
c
451	n = narr(iarg(1))
c
	if (n.le.0) return
c
	narr(iarg(3)) = n
c
        if (c2.eq.'+') then
            do 460 i=1,n
460         arr(i,iarg(3))=arr(i,iarg(1))+s
        else if (c2.eq.'-') then
            do 470 i=1,n
470         arr(i,iarg(3))=arr(i,iarg(1))-s
        else if (c2.eq.'*') then
            do 480 i=1,n
480         arr(i,iarg(3))=arr(i,iarg(1))*s
        else if (c2.eq.'/') then
            if (s.eq.0.)go to 1001
            s=1./s
            do 490 i=1,n
490         arr(i,iarg(3))=arr(i,iarg(1))*s
        else if (c2.eq.'^') then
            do 500 i=1,n
500         if (arr(i,iarg(1)).gt.0.) arr(i,iarg(3))=arr(i,iarg(1))**s
	else if (c2.eq.'=') then
            do 510 i=1,n
510         arr(i,iarg(3))=s
	else
	    go to 1001
        end if
c
	return
1001	return 1
c
	end
