
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

	subroutine xyplot(input,istart,nin,x,y,z,n,
     &			  itype,jth,jsym,plot_symbol,hp,*)
        save
c
c	Plot y(x).
c
	character*(*) input
	dimension x(1),y(1),z(1)
	character*1 plot_symbol
c
	save i1plot,i2plot
c
	icplot = 0
	go to 1
c
	entry xyplotc(input,istart,nin,x,y,z,n,
     &                itype,jth,jsym,plot_symbol,hp,*)
	icplot = 1
c
1	if (input(2:2).ne.'x') then
            i1plot = 1
            i2plot = n
            call readiq(input(istart:nin),2,
     &                  i1plot,i2plot,idum,idum,*1001)
	else
	    x1 = -1.e30
	    x2 = 1.e30
            call readrq(input(istart:nin),2,
     &                  x1,x2,dum,dum,*1001)
            xx1 = min(x1,x2)
            x2 = max(x1,x2)
	    x1 = xx1
	    i1plot = 0
	    do i=1,n
	        if (i1plot.eq.0.and.x(i).ge.x1) i1plot = i
		if (x(i).gt.x2) then
		    i2plot = i-1
		    go to 11301
		end if
            end do
11301       continue
c            write(6,*)'x1, x2 = ',x1, x2
c            write(6,*)'i1, i2 = ',i1plot,i2plot
        end if
c
        i1plot = max(1,min(n,i1plot))
        i2plot = max(i1,min(n,i2plot))
c        write(6,*)'i1, i2 = ',i1plot,i2plot
c
	entry unplot(x,y,z,n,itype,jth,jsym,plot_symbol,hp)
c
c       Note that this entry inherits the previous i1plot, i2plot!
c
	call devon
        nnplot = i2plot - i1plot + 1
	if (nnplot.le.0) go to 1001
c
301     if (plot_symbol.ne.' ') then
c
c           Encode the symbol.
c
	    if (plot_symbol.eq.'#') then
		js = -1
	    else
		js = -ichar(plot_symbol)
	    end if
	else
	    js = abs(jsym)
	end if
c
	if (itype.eq.0) then
	    if (icplot.eq.0) then
                call mline(x(i1plot),y(i1plot),nnplot,jth,js,hp)
	    else
                call mlinec(x(i1plot),y(i1plot),z(i1plot),
     &				nnplot,jth,js,hp)
	    end if
        else
	    if (icplot.eq.0) then
                call dline(x(i1plot),y(i1plot),nnplot,jth,js,hp)
	    else
                call dlinec(x(i1plot),y(i1plot),z(i1plot),
     &				nnplot,jth,js,hp)
	    end if
        end if
c
	return
1001	return 1
c
	end
