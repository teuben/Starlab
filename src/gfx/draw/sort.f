
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
c       Sort routines, from "Numerical Recipes."
c       
c------------------------------------------------------------------------
        
        subroutine sort(n,a)
c       
c       Sort the real array a.
c       
        implicit real (a-h,o-z)
        dimension a(n)
c
        if (n.le.1) return
c       
        l=n/2+1
        ir=n
10      continue
        if(l.gt.1)then
            l=l-1
            aa=a(l)
        else
            aa=a(ir)
            a(ir)=a(1)
            ir=ir-1
            if(ir.eq.1)then
                a(1)=aa
                return
            end if
        end if
        i=l
        j=l+l
20      if(j.le.ir)then
            if(j.lt.ir)then
                if(a(j).lt.a(j+1))j=j+1
            end if
            if(aa.lt.a(j))then
                a(i)=a(j)
                i=j
                j=j+j
            else
                j=ir+1
            end if
            go to 20
        end if
        a(i)=aa
        go to 10
c       
        end


        subroutine sort2(n,a,b)
c       
c       Sort the real array a, carrying the real array b along
c       in the process.
c       
        implicit real (a-h,o-z)
        dimension a(n),b(n)
c       
        if (n.le.1) return
c       
        l=n/2+1
        ir=n
10      continue
        if(l.gt.1)then
            l=l-1
            aa=a(l)
            bb = b(l)
        else
            aa=a(ir)
            bb = b(ir)
            a(ir)=a(1)
            b(ir) = b(1)
            ir=ir-1
            if(ir.eq.1)then
                a(1)=aa
                b(1) = bb
                return
            end if
        end if
        i=l
        j=l+l
20      if(j.le.ir)then
            if(j.lt.ir)then
                if(a(j).lt.a(j+1))j=j+1
            end if
            if(aa.lt.a(j))then
                a(i)=a(j)
                b(i) = b(j)
                i=j
                j=j+j
            else
                j=ir+1
            end if
            go to 20
        end if
        a(i)=aa
        b(i) = bb
        go to 10
c       
        end


        subroutine sort3(n,a,b,c)
c       
c       Sort the real array a, carrying the real arrays b and c along
c       in the process.
c       
        implicit real (a-h,o-z)
        dimension a(n),b(n),c(n)
c       
        if (n.le.1) return
c       
        l=n/2+1
        ir=n
10      continue
        if(l.gt.1)then
            l = l-1
            aa = a(l)
            bb = b(l)
            cc = c(l)
        else
            aa = a(ir)
            bb = b(ir)
            cc = c(ir)
            a(ir) = a(1)
            b(ir) = b(1)
            c(ir) = c(1)
            ir = ir-1
            if(ir.eq.1)then
                a(1)=aa
                b(1) = bb
                c(1) = cc
                return
            end if
        end if
        i=l
        j=l+l
20      if(j.le.ir)then
            if(j.lt.ir)then
                if (a(j).lt.a(j+1)) j=j+1
            end if
            if(aa.lt.a(j))then
                a(i) = a(j)
                b(i) = b(j)
                c(i) = c(j)
                i=j
                j=j+j
            else
                j=ir+1
            end if
            go to 20
        end if
        a(i) = aa
        b(i) = bb
        c(i) = cc
        go to 10
c       
        end
