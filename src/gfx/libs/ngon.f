
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

        subroutine ngon(rcall,scall,h,n,thoff)
        save
c
c       Draws |n|-sided polygons  (n = 0 ==> "."
c                                      1     "+"
c                                      2     "x".
c       or "stars" (for n < 0)
c
        real*4 rr,ss,r,s,h,thoff
        integer*4 n
	common /findex/ index
	dimension rpoly(0:100),spoly(0:100),cth1(0:100),sth1(0:100)
c
	data n1/0/
c
        r=rcall
        s=scall
        go to 1
c
        entry ungon(rcall,scall,h,n,thoff)
        entry userngon(rcall,scall,h,n,thoff)
        call fr inches(rcall,scall,r,s)
c
1       if(h.lt.0..and.n.ne.0)then
            hh=-h
            rr=r+.7071*hh
            ss=s+.7071*hh
        else
            hh=h
            rr=r
            ss=s
        end if
c
        if(n.eq.0.or.h.eq.0.)then
            call point(rr,ss)
        else if(abs(n).eq.1)then
            call plot(rr-hh,ss,3)
            call plot(rr+hh,ss,2)
            call plot(rr,ss-hh,3)
            call plot(rr,ss+hh,2)
        else if(abs(n).eq.2)then
            hh=.7071*hh
            call plot(rr-hh,ss-hh,3)
            call plot(rr+hh,ss+hh,2)
            call plot(rr-hh,ss+hh,3)
            call plot(rr+hh,ss-hh,2)
        else
            nn=abs(n)
	    nn1=nn
	    if(n.lt.0)then
		nn2=nn/2
		ieven=0
	    	if(2*nn2.eq.nn)ieven=1
		if(ieven.eq.1)nn1=nn2
	    end if
            th0=thoff
            dth=360./float(nn)
            th0=th0-90.-.5*dth
            th0=.017453*th0
            dth=.017453*dth
            th=th0-dth
            ipen=3
            do 100 i=0,nn1
		if(nn.ne.n1)then
                    th=th+dth
		    cth1(i)=cos(th)
		    sth1(i)=sin(th)
		end if
		hc=hh*cth1(i)
		hs=hh*sth1(i)
                r1=rr+hc
                s1=ss+hs
                if(n.gt.0)then
		    if(index.lt.0)then
                        call plot(r1,s1,ipen)
                        ipen=2
		    else
			rpoly(i)=r1
			spoly(i)=s1
		    end if
                else
		    if(ieven.eq.0)then
                	call plot(rr,ss,3)
		    else
			call plot(rr-hc,ss-hs,3)
		    end if
                    call plot(r1,s1,2)
                end if
100         continue
	    if(index.ge.0.and.n.gt.0)call polyfill(rpoly,spoly,n+1)
	    n1=nn
        end if
c
        end
