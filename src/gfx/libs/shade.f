
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
c
c
c---	subroutine shade
c
c---	creates complete cross-hatch density map
c	note:  negative cross hatches are circles
c
c---	created  11.aug.80
c
        subroutine shade (a,m,n,sc,xmin,xmax,ymin,ymax,msiz,nsiz)
        save
c---	a	real array of data
c	m,n	subset of a to be plotted
c	sc	scale parameter:
c		sc < 0	then k = -sc*ln(a)
c		sc = 0	then k =  sc*a
c	xmin,xmax,ymin,ymax    delineates plotting area
c	msiz,nsiz   dimensions of a
c
        real*4 a(msiz,nsiz)
        da=0.3490658504
        dx=(xmax-xmin)/m
        dy=(ymax-ymin)/n
        if (sc) 10,30,20
10      scx=abs(sc)
        y=ymin-dy/2.
        yp=ymin
        ym=ymin-dy
        do 19 j=1,n
            y=y+dy
            yp=yp+dy
            ym=ym+dy
            x=xmin-dx/2.
            xp=xmin
            xm=xmin-dx
            do 19 i=1,m
                x=x+dx
                xp=xp+dx
                xm=xm+dx
                k=nint(scx*alog(a(i,j)))
                if (k) 15,19,11
11              ddx=dx/(k+1.)
                ddy=dy/(k+1.)
                do 14 kk=1,k
                    xx=xm+kk*ddx
                    yy=ym+kk*ddy
                    call plot(xm,yy,3)
                    call plot(xp,yy,2)
                    call plot(xx,ym,3)
                    call plot(xx,yp,2)
14              continue
                go to 19
15              ddx=.5*dx/(-k+1.)
                ddy=.5*dy/(-k+1.)
                do 17 kk=1,-k
                    call plot(x,y+kk*ddy,3)
                    do 17 kkk=1,18
                        xx=x+kk*ddx*sin(kkk*da)
                        yy=y+kk*ddy*cos(kkk*da)
                        call plot(xx,yy,2)
17              continue
19      continue
        go to 30
20      y=ymin-dy/2.
        yp=ymin
        ym=ymin-dy
        do 29 j=1,n
            y=y+dy
            yp=yp+dy
            ym=ym+dy
            x=xmin-dx/2.
            xp=xmin
            xm=xmin-dx
            do 29 i=1,m
                x=x+dx
                xp=xp+dx
                xm=xm+dx
                k=nint(sc*a(i,j))
                if (k) 25,29,21
21              ddx=dx/(k+1.)
                ddy=dy/(k+1.)
                do 24 kk=1,k
                    xx=xm+kk*ddx
                    yy=ym+kk*ddy
                    call plot(xm,yy,3)
                    call plot(xp,yy,2)
                    call plot(xx,ym,3)
                    call plot(xx,yp,2)
24              continue
                go to 29
25              ddx=.5*dx/(-k+1.)
                ddy=.5*dy/(-k+1.)
                do 27 kk=1,-k
                    call plot(x,y+kk*ddy,3)
                    do 27 kkk=1,18
                        xx=x+kk*ddx*sin(kkk*da)
                        yy=y+kk*ddy*cos(kkk*da)
                        call plot(xx,yy,2)
27              continue
29      continue
30      return
        end
