
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

	subroutine varith(input,istart,nin,arr,nmax,narr,iprompt,*)
        save
c
c	Vector arithmetic.
c
	character*(*) input
	dimension arr(nmax,3),narr(3)
c
	character*1 c2
	dimension iarg(3)
c
	if (nin.lt.istart) go to 1001
	c2 = input(2:2)
c
        ibl=1
        narg=0
c
c	Decode the argument list
c
        do 70 i=istart,nin
            if (input(i:i).eq.' '.or.input(i:i).eq.',') then
                ibl=1
            else
                if (ibl.eq.1) then
                    narg=narg+1
                    if (narg.ge.2.and.c2.eq.'c') then
                        call readrq(input(i:nin),2,
     &                              x1,x2,dum,dum,*69)
                        go to 71
c
69                      call readrq(input(i:nin),1,
     &                              x1,dum,dum,dum,*1001)
                        x2=1.e30
                        go to 71
c
		    else if (narg.ge.2.and.c2.eq.'o') then
			call readiq(input(i:nin),1,
     &                              noff,idum,idum,idum,*1001)
			go to 71
c
		    else if (narg.ge.2.and.c2.eq.'!') then
			call readiq(input(i:nin),1,
     &                              nskip,idum,idum,idum,*1001)
                        if (nskip.lt.0) go to 1001
			go to 71
                    else
                        call decode(input(i:i),iarg(narg),*1001)
                        if (narg.eq.3)go to 71
                    end if
                    ibl=0
                end if
            end if
70      continue
71      continue
c
        if (narg.lt.2) then
            if (c2.ne.'i'
	1	    .and.c2.ne.'a'
	2	    .and.c2.ne.'c'
	3	    .and.c2.ne.'o'
	4	    .and.c2.ne.'\\') go to 1001
        else if (narg.eq.2) then
            iarg(3)=iarg(2)
        end if
c
	n = narr(iarg(1))
c
c	Apply the appropriate operation:
c
        if (c2.eq.'+') then
            do 75 i=1,n
75          arr(i,iarg(3))=arr(i,iarg(1)) + arr(i,iarg(2))
	    narr(iarg(3)) = n
        else if (c2.eq.'-') then
            do 76 i=1,n
76          arr(i,iarg(3))=arr(i,iarg(1)) - arr(i,iarg(2))
	    narr(iarg(3)) = n
        else if (c2.eq.'*') then
            do 77 i=1,n
77          arr(i,iarg(3))=arr(i,iarg(1)) * arr(i,iarg(2))
	    narr(iarg(3)) = n
        else if (c2.eq.'/') then
            nerr=0
            do 78 i=1,n
                s=arr(i,iarg(2))
                if (s.eq.0.) then
                    nerr=nerr+1
                else
                    arr(i,iarg(3))=arr(i,iarg(1))/s
                end if
78          continue
	    narr(iarg(3)) = n
            if (nerr.gt.0) then
                call devoff
                if (iprompt.eq.1) write(6,'(i5,'' error(s)'')')nerr
            end if
        else if (c2.eq.'^') then
            nerr=0
            do 79 i=1,n
                if (arr(i,iarg(1)).le.0.) then
                    nerr=nerr+1
                else
                    arr(i,iarg(3))=arr(i,iarg(1))**arr(i,iarg(2))
                end if
79          continue
	    narr(iarg(3)) = n
            if (nerr.gt.0) then
                call devoff
                if (iprompt.eq.1) write(6,'(i5,'' error(s)'')')nerr
            end if
        else if (c2.eq.'i') then
            do 80 i=1,n
80          if (arr(i,iarg(1)).ne.0.) arr(i,iarg(1))=1./arr(i,iarg(1))
        else if (c2.eq.'a') then
            do 81 i=1,n
81          arr(i,iarg(1))=abs(arr(i,iarg(1)))
        else if (c2.eq.'c') then
            do 82 i=1,n
82          arr(i,iarg(1))=min(x2,max(x1,arr(i,iarg(1))))
        else if (c2.eq.'=') then
            do 83 i=1,n
83          arr(i,iarg(2))=arr(i,iarg(1))
	    narr(iarg(2)) = n
	else if (c2.eq.'o') then
	    if (noff.gt.0) then
		if (narr(iarg(1))+noff.gt.nmax) then
		    if (iprompt.eq.1) write(6,*)'Array overflow'
		    go to 1001
		else
		    do 84 i=narr(iarg(1)),1,-1
84		    arr(i+noff,iarg(1)) = arr(i,iarg(1))
		    do 85 i=1,noff
85		    arr(i,iarg(1)) = 0.
		end if
	    else if (noff.lt.0) then
		do 86 i=1-noff,narr(iarg(1))
86		arr(i+noff,iarg(1)) = arr(i,iarg(1))
	    end if
	    narr(iarg(1)) = narr(iarg(1)) + noff
        else if (c2.eq.'\\') then
            do 87 i=1,n/2
		temp = arr(i,iarg(1))
		arr(i,iarg(1)) = arr(n+1-i,iarg(1))
		arr(n+1-i,iarg(1)) = temp
87	    continue
        else if (c2.eq.'>') then
            do 88 i=1,n
88          arr(i,iarg(3))=max(arr(i,iarg(1)), arr(i,iarg(2)))
	    narr(iarg(3)) = n
        else if (c2.eq.'<') then
            do 89 i=1,n
89          arr(i,iarg(3))=max(arr(i,iarg(1)), arr(i,iarg(2)))
	    narr(iarg(3)) = n
	else if (c2.eq.'!') then
            if (nskip.gt.0) then
                n1 = 0
                do 90 i=1,n,nskip+1
                    n1 = n1 + 1
                    arr(n1,iarg(1)) = arr(i,iarg(1))
90              continue
                narr(iarg(1)) = n1
            end if
        else
            go to 1001
        end if
c
	return
1001	return 1
c
	end


	subroutine negate(n,a)
        save
c
c	Replace array a by -a.
c
	dimension a(1)
c
	do 10 i=1,n
10	a(i) = -a(i)
c
	end
