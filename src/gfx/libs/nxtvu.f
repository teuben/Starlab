
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

	block data setupvu
        save
	common /vuopt/ iopt
	data iopt/1/
	end


	subroutine setvu(jopt)
        save
	common /vuopt/ iopt
	iopt = jopt
	end


        subroutine nxtvu(ic,x,y,n,ier)
        save
c
c	Plot routine used in conjunction with plt3d - includes hidden
c	line removal.
c
c	This routine is quite choosy.  It MUST traverse the screen from
c	left to right, or it loses track of the hidden lines.
c
        dimension x(n),y(n)
c
        parameter (nn = 2000)
        common /nxtv1/ xx(nn),yy(nn),kk,ll
	common /vuopt/ iopt
	save /vuopt/
c
        if (ic.eq.0) then
            if (n.gt.nn) then
                ier = 1
                return
	    end if
c
            ll = nn-n+1
            i = ll
            xx(i) = x(1)
            yy(i) = y(1)
            call plot(xx(ll),yy(ll),3)
            do 10 j = 2,n
                i = i+1
                xx(i) = x(j)
                yy(i) = y(j)
                call plot(xx(i),yy(i),2)
 10         continue
            ier = 0
            return
	end if
c
        if (ier.ne.0) return
c
        ii = 1
        jj = ll
        kk = 0
        ya0 = y(1)
        yb0 = yy(ll)
        if (x(1).gt.xx(ll)) go to 70
        call plot(x(1),ya0,3)
c
 40     call outp(x(ii),y(ii),ier)
        if (ii.eq.n) go to 360
        ii = ii+1
        ya0 = y(ii)
        if (x(ii).le.xx(ll)) then
            call plot(x(ii),ya0,2)
            go to 40
	end if
c
        ii = ii-1
        xl = x(ii)
        yl = y(ii)
        ya0 = alin(x(ii),x(ii+1),y(ii),y(ii+1),xx(ll))
        x0 = xx(ll)
        if ((iopt.eq.1.and.ya0.gt.yb0)
     &		.or.(iopt.eq.2.and.ya0.lt.yb0)) then
	    iov0 = 1
	else
            call plot(x0,ya0,2)
            call outp(x0,ya0, ier)
            call outp(x0,yb0, ier)
	    iov0 = 0
	end if
	go to 120
c
 70     call outp(xx(jj),yy(jj),ier)
        if (jj.eq.nn) go to 380
        jj = jj+1
        yb0 = yy(jj)
        if (x(1).ge.xx(jj)) go to 70
c
        jj = jj-1
        yb0 = alin(xx(jj),xx(jj+1),yy(jj),yy(jj+1),x(1))
        x0 = x(1)
        if ((iopt.eq.1.and.ya0.le.yb0)
     &		.or.(iopt.eq.2.and.ya0.ge.yb0)) then
	    iov0 = 0
	else
            call outp (x0,yb0,ier)
            call outp(x0,ya0,ier)
            xl = x0
            yl = ya0
            iov0 = 1
	end if
c
 120    if (ii.eq.n) go to 300
        if (jj.eq.nn) go to 310
        if (x(ii+1).le.xx(jj+1)) then
            isw = +1
            ii = ii+1
            x1 = x(ii)
            ya1 = y(ii)
            yb1 = alin(xx(jj),xx(jj+1),yy(jj),yy(jj+1),x1)
        else
            if (xx(jj+1).ge.x(n)) go to 340
            isw = -1
            jj = jj+1
            x1 = xx(jj)
            ya1 = alin(x(ii),x(ii+1),y(ii),y(ii+1),x1)
            yb1 = yy(jj)
	end if
        if ((iopt.eq.1.and.ya1.le.yb1)
     &		.or.(iopt.eq.2.and.ya1.ge.yb1)) go to 160
        iov1 = 1
        if (iov0.eq.0) go to 170
 150    if (isw.eq.-1) go to 200
        call outp (x1,ya1,ier)
        call plot(xl,yl,3)
        call plot(x1,ya1,2)
        xl = x1
        yl = ya1
        go to 200
c
 160    iov1 = 0
        if (iov0.eq.0) go to 190
 170    frac = (yb0-ya0)/(ya1-yb1+yb0-ya0)
        xi = (x1-x0)*frac+x0
        yi = (ya1-ya0)*frac+ya0
        call outp(xi,yi,ier)
        if (iov0.eq.0) go to 180
        call plot(xl,yl,3)
        call plot(xi,yi,2)
        xl = xi
        yl = yi
        go to 190
c
 180    xl = xi
        yl = yi
        go to 150
c
 190    if (isw.eq.+1) go to 200
        call outp(xx(jj),yy(jj),ier)
 200    if (ier.ne.0) return
        x0 = x1
        ya0 = ya1
        yb0 = yb1
        iov0 = iov1
        go to 120
c
 310    x1 = xx(nn)
        ya1 = alin(x(ii),x(ii+1),y(ii),y(ii+1),x1)
        yb1 = yy(nn)
        if ((iopt.eq.1.and.ya1.gt.yb1)
     &		.or.(iopt.eq.2.and.ya1.lt.yb1)) go to 320
        call outp(x1,yb1,ier)
        call outp(x1,ya1,ier)
        call plot(x1,ya1,3)
        go to 330
c
 380    ii = 1
 320    call plot(x(ii),y(ii),3)
 330    if (ii.eq.n) go to 400
        ii = ii+1
        call outp(x(ii),y(ii),ier)
        call plot(x(ii),y(ii),2)
        go to 330
c
 300    if (jj.eq.nn) go to 400
 340    x1 = x(n)
        ya1 = y(n)
        yb1 = alin(xx(jj),xx(jj+1),yy(jj),yy(jj+1),x1)
        if ((iopt.eq.1.and.ya1.le.yb1)
     &		.or.(iopt.eq.2.and.ya1.ge.yb1)) go to 350
        call outp(x1,ya1,ier)
        call outp(x1,yb1,ier)
        call plot(xl,yl,3)
        call plot (x1,ya1,2)
350     if (jj.eq.nn) go to 400
        jj = jj+1
        call outp(xx(jj),yy(jj),ier)
        go to 350
c
 360    jj = 0
        go to 350
c
 400    ll = nn-kk+1
        i = ll
        do 410 j = 1,kk
            xx(i) = xx(j)
            yy(i) = yy(j)
            i = i+1
 410    continue
c
        end
