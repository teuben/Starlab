
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
c       Integration and differentiation of the input arrays.
c
c------------------------------------------------------------------------

        subroutine integrate(x,y,z,n,iprompt,*)
c
c       Integrate y with respect to x, using the trapezoid rule,
c       assuming that there are no pathologies in the data, placing
c       the result in z.
c
        dimension x(n),y(n),z(n)
c
        if (n.le.0) return
c
        if (n.eq.1) then
c
c           Result = 0 if only one point is supplied.
c
            z(1) = 0.
c
        else
            sum = 0.
            do i=2,n
                sum1 = sum + .5*(y(i)+y(i-1))*(x(i)-x(i-1))
                z(i-1) = sum
                sum = sum1
            end do
            z(n) = sum
        end if
c
        end


        subroutine differentiate(x,y,z,n,iprompt,*)
c
c       Differentiate y with respect to x, assuming that there are no
c       pathologies in the data, placing the result in z.  Note that
c       the array "z" may really be x or y.
c
        dimension x(n),y(n),z(n)
c
        if (n.le.0) return
c
        nerr = 0
        if (n.eq.1) then
c
c           Result = 0 if only one point is supplied.
c
            z(1) = 0.
c
        else
c
c           Don't (can't) center the first and last derivatives.
c           The rest are centered and hence second-order accurate.
c
            zz = (y(2) - y(1))/(x(2) - x(1))
c
            do i=2,n-1
                zz0 = zz
c
                dx = x(i+1) - x(i-1)
                dy = y(i+1) - y(i-1)
c
                if (dx.eq.0.) then
                    zz = 0.
                    nerr = nerr + 1
                else
                    zz = dy/dx
                end if
c
c               Set z(i-1) at the end of the i loop, since location
c               i-1 is not going to be used again.
c
                z(i-1) = zz0
            end do
c
            z(n) = (y(n) - y(n-1))/(x(n) - x(n-1))
            z(n-1) = zz
c
        end if
c
        if (nerr.gt.0) then
            if (iprompt.ne.0)then
                call devoff
                write(6,'(i6,a)'),nerr,' errors'
            end if
            return 1
        end if
c
        end
