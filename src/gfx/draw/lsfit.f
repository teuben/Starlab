
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
        
c------------------------------------------------------------------------
c       
c       Least-squares fitting routines.
c       
c------------------------------------------------------------------------
        
        subroutine lsfit1(x,y,n,w,nw,a,b)
        implicit real (a-h,o-z)
        save
c
        dimension x(1),y(1),w(1)
        common/resid/resi
c       
c       Return least-squares best linear fit to y(x) (n points,
c       weighted by w) of the form y = a*x + b, where a is
c       specified and b is variable.
c       
        if(n.le.1)then
            b=0.
            resi=0.
            return
        end if
c       
        s1=0.
        s2=0.
        do i=1,n
c           
            if (nw.eq.1) then
                ww = w(1)
            else
                ww = w(i)
            end if
c           
            s1=s1+ww
            s2=s2+ww*(y(i)-a*x(i))
        end do
        b=s2/s1
c       
        resi=0.
        do i=1,n
c           
            if (nw.eq.1) then
                ww = w(1)
            else
                ww = w(i)
            end if
c           
            resi=resi+ww*(y(i)-a*x(i)-b)**2
        end do
        resi=sqrt(resi/s1)
c       
        dy=a*(x(1)-x(n))
        if(dy.ne.0.)then
            resi=resi/abs(dy)
        else
            resi=1.
        end if
c       
        end


        subroutine lsfit2(x,y,n,w,nw,a,b)
        implicit real (a-h,o-z)
        save
c
        dimension x(1),y(1),w(1)
        common/resid/resi
c       
c       Return least-squares best linear fit to y(x) (n points,
c       weighted by w) of the form y = a*x + b, where a and b
c       are both variable.
c       
        if(n.le.1)then
            a=0.
            b=0.
            resi=0.
            return
        end if
c       
c       We will used NORMALIZED numbers for the fitting, to avoid
c       problems with single-precision roundoff...
c
        xo = x(1)
        xfac = 1./(x(n) - xo)
        yo = y(1)
        yfac = 1./(y(n) - yo)
c
        s1=0.
        sx=0.
        sy=0.
        sxx=0.
        sxy=0.
        do i=1,n
c           
            if (nw.eq.1) then
                ww = w(1)
            else
                ww = w(i)
            end if
c
            xx = (x(i) - xo)*xfac
            yy = (y(i) - yo)*yfac
c
            s1=s1+ww
            sx=sx+ww*xx
            sy=sy+ww*yy
            sxx=sxx+ww*xx**2
            sxy=sxy+ww*xx*yy
        end do
        d=sxx*s1-sx**2
c
        a=0.
        b=0.
        if(d.ne.0.)then
            d=1./d
            a=d*(s1*sxy-sx*sy)
            b=d*(sxx*sy-sx*sxy)
        end if
c
c       Now undo the normalization...
c
        b = (b-a*xfac*xo)/yfac + yo
        a = a*xfac/yfac
c       
        resi=0.
        do i=1,n
c           
            if (nw.eq.1) then
                ww = w(1)
            else
                ww = w(i)
            end if
c           
            resi=resi+ww*(y(i)-a*x(i)-b)**2
        end do
        resi=sqrt(resi/s1)
c       
        dy=a*(x(1)-x(n))
        if(dy.ne.0.)then
            resi=resi/abs(dy)
        else
            resi=1.
        end if
c       
        end
