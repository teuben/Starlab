
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
        subroutine xyzplot(input,istart,nin,x,y,z,n,jth,
     &                     jsym,plot_symbol,hp,*)
        save
c
c	Plot y(x), clipped according to z.  Place an ngon symbol
c       determined by jsym and plot_symbol at every |jth|-th point
c       (see mline).
c
	character*(*) input
	dimension x(1),y(1),z(1)
	character*1 plot_symbol,sim
c
	common /scales/ xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
c
        call readrq(input(istart:nin),2,
     &              z1,z2,dum,dum,*1001)
        call devon
c
        jj = abs(jth)
c
	in = 0
        do 290 i=1,n
            if (x(i).ge.xl.and.x(i).le.xr
     &            .and.y(i).ge.ybot.and.y(i).le.ytop
     &		    .and.z(i).ge.z1.and.z(i).le.z2) then
c
                call fr inches(x(i),y(i),rloc,sloc)
c
                if (jth.ne.0
     &                .and.(jsym.ne.0.or.plot_symbol.ne.' ')) then
                    if (jj*(i/jj).eq.i) then
			if (jsym.ne.0) then
			    call ngon(rloc,sloc,.5*hp,jsym,0.)
			else
			    if (plot_symbol.ne.'#') then
				sim = plot_symbol
			    else
c                         
c                               Encode the number of the point in sim.
c                         
				ii = i
				do while (ii.gt.61)
				    ii = ii - 61
				end do
				if (ii.le.9) then
				    sim = char(48+ii)
				else if (ii.le.35) then
				    sim = char(87+ii)
				else
				    sim = char(55+ii)
				end if
			    end if
			    call simbol(rloc,sloc,hp,sim,0.,-1)
			end if
		    end if
                end if
c
	        if (jth.ge.0.and.in.eq.1) then
		    call plotin(rloc0,sloc0,3)
		    call plotin(rloc,sloc,2)
		end if
		in = 1
		rloc0 = rloc
		sloc0 = sloc
	    else
		in = 0
            end if
290     continue
c
	return
1001	return 1
c
	end
