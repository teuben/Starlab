
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

	subroutine fitxy(input,istart,nin,x,y,z,w,n,iprompt,ier)
        save
c
c	Perform least-squares fit to y(x).
c       NOTE: It is assumed that x is sorted!
c
	character*(*) input
	dimension x(*),y(*),z(*),w(*)
c       *****  For now, w really is of dimension 1.  *****
c
        common/resid/resi
c
c	Equal weights for all points.
c
	nw = 1
	go to 10
c
	entry fitxyz(input,istart,nin,x,y,z,w,n,iprompt,ier)
c
c	Allow the weights to be specified in the z array...
c
	nw = n
c
10      if (istart.gt.nin) then
            i1=1
            nn=n
            call minmax(x,n,x1,x2)
        else
            call readrq(input(istart:nin),2,
     &                  x1,x2,dum,dum,*1001)
            dx = x(2) - x(1)
            i1 = 0
            do 175 i=1,n
c
                if (i.gt.1.and.(x(i)-x(i-1))*dx.lt.0.) then
c
c                   The x-array is not monotonic.
c
                    write(6,*)'The x array must be ordered.'//
     &                      '  Use "so2 x y" to sort x and reorder y.'
                    go to 1001
                end if
c
                if ((x2-x(i))*(x(i)-x1).ge.0.) then
                    if (i1.eq.0) i1 = i
                else if (i1.gt.0) then
                    go to 176
                end if
c
175         continue
            i=n+1
176         i1=max(1,i1)
            nn=i-i1
            if (nn.le.1) then
                write(6,*)'No points in specified range.'
                go to 1001
            end if
c
            x1 = x(i1)
            x2 = x(i1+nn-1)
        end if
c
        if (nw.gt.1) then
c
c           Make sure the weights are non-negative.
c
            do 177 i=1,n
                if (z(i).lt.0.) nw = 1
177         continue
        end if
c
	if (nw.eq.1) then
            call lsfit2(x(i1),y(i1),nn,w(1),1,aa,bb)
c           *****  For now, w really is of dimension 1.  *****
	else
            call lsfit2(x(i1),y(i1),nn,z(i1),nn,aa,bb)
	end if
c
c	Print some basic data on the fit.
c       --------------------------------
c
        if (iprompt.eq.1) then
            call devoff
            write(6,178)aa,bb,-bb/aa,resi
178         format(' slope =',1pe11.3,',  intercepts =',2e11.3,
     &             ',  residual =',e11.3)
        end if
c
        if (input(3:3).eq.'p') then
c
c	    Plot the fit.
c
            call fr inches(x1,aa*x1+bb,rr,ss)
            call devon
            call plotin(rr,ss,3)
            call fr inches(x2,aa*x2+bb,rr,ss)
            call plotin(rr,ss,2)
        end if
c
        ier = 0
	return
c
1001	ier = 1
        return
c
	end
