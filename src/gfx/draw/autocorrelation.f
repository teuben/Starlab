
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

        subroutine autocorrel(a,n,w)
c
c       Replace a by its autocorrelation function.
c
        real a(n),w(0:n)
c
        if (n.le.1) return
c
c       First take out the mean signal.
c
        aver = 0.
        do i=1,n
            aver = aver + a(i)
        end do
c
        aver = aver/n
        do i = 1,n
            a(i) = a(i) - aver
        end do
c
        n2 = n/2
c
        do i=0,n2
            w(i) = 0.
            do j=1,n-i
                w(i) = w(i) + a(j)*a(j+i)
            end do
            w(i) = w(i) / (n-i)
        end do
c
c       Normalize by the i=0 term.
c
        do i=1,n2
            a(i) = w(i)/w(0)
        end do
c
        n = n2
c
        end
        
