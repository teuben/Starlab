
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



        subroutine contor2 (a,x,y,m,n,v,iv,userplot)
        save
        external userplot
        dimension ix(50,250),xb(50,250),iy(50),yl(50),dx(250)
c
c       subroutine userplot converts from (x,y) coordinates to
c       the desired rectangular output grid. if the original grid
c       is already rectangular and correctly scaled (see mcinit),
c       just set userplot = plot in the call to contor.
c
        real v(iv),a(m,n),x(m),y(n)
        mm=m-1
        nn=n-1
        do 10 i=1,mm
            dx(i)=x(i+1)-x(i)
            do 10 k=1,iv
                ix(k,i)=0
10      continue
        do 60 j=1,nn
            dely=y(j+1)-y(j)
            do 15 k=1,iv
15          iy(k)=0
            do 50 i=1,mm
                delx=dx(i)
                v1=a(i,j)
                v2=a(i+1,j)
                v3=a(i,j+1)
                v4=a(i+1,j+1)
                do 40 k=1,iv
                    val=v(k)
                    if (val.lt.-1.e10) go to 40
                    icase=1
                    if (val.gt.v1) icase=icase+1
                    if (val.gt.v2) icase=icase+2
                    if (val.gt.v3) icase=icase+4
                    if (val.gt.v4) icase=9-icase
                    go to (40,24,25,26,27,28,29,30),icase
24                  x0=x(i)
                    if(iy(k).ne.0)then
                        y0=yl(k)
                    else
                        y0=y(j)+dely*(val-v1)/(v3-v1)
                    end if
                    if(ix(k,i).ne.0)then
                        x1=xb(k,i)
                    else
                        x1=x(i)+delx*(val-v1)/(v2-v1)
                    end if
                    y1=y(j)
                    go to 35
25                  if(ix(k,i).ne.0)then
                        x0=xb(k,i)
                    else
                        x0=x(i)+delx*(val-v1)/(v2-v1)
                    end if
                    y0=y(j)
                    x1=x(i+1)
                    y1=y(j)+dely*(val-v2)/(v4-v2)
                    go to 35
26                  x0=x(i)
                    if(iy(k).ne.0)then
                        y0=yl(k)
                    else
                        y0=y(j)+dely*(val-v1)/(v3-v1)
                    end if
                    x1=x(i+1)
                    y1=y(j)+dely*(val-v2)/(v4-v2)
                    go to 35
27                  x0=x(i)
                    if(iy(k).ne.0)then
                        y0=yl(k)
                    else
                        y0=y(j)+dely*(val-v1)/(v3-v1)
                    end if
                    x1=x(i)+delx*(val-v3)/(v4-v3)
                    y1=y(j+1)
                    go to 35
28                  if(ix(k,i).ne.0)then
                        x0=xb(k,i)
                    else
                        x0=x(i)+delx*(val-v1)/(v2-v1)
                    end if
                    y0=y(j)
                    x1=x(i)+delx*(val-v3)/(v4-v3)
                    y1=y(j+1)
                    go to 35
29                  x0=x(i)
                    if(iy(k).ne.0)then
                        y0=yl(k)
                    else
                        y0=y(j)+dely*(val-v1)/(v3-v1)
                    end if
                    if(ix(k,i).ne.0)then
                        x1=xb(k,i)
                    else
                        x1=x(i)+delx*(val-v1)/(v2-v1)
                    end if
                    y1=y(j)
                    call userplot(x0,y0,3)
                    call userplot(x1,y1,2)
30                  x0=x(i)+delx*(val-v3)/(v4-v3)
                    y0=y(j+1)
                    x1=x(i+1)
                    y1=y(j)+dely*(val-v2)/(v4-v2)
35                  call userplot(x0,y0,3)
                    call userplot(x1,y1,2)
                    if(x1.eq.x(i+1))then
                        iy(k)=1
                        yl(k)=y1
                    else
                        iy(k)=0
                    end if
                    if(y0.eq.y(j+1))then
                        ix(k,i)=1
                        xb(k,i)=x0
                    else if(y1.eq.y(j+1))then
                        ix(k,i)=1
                        xb(k,i)=x1
                    else
                        ix(k,i)=0
                    end if
40              continue
50          continue
60      continue
        end
