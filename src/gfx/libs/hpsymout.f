
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
        subroutine hp sym out(string,islant)
        save
        character*(*) string
        character*1 ss
        parameter (slant=.2)
c
        n=len(string)
        nn=n
        if(islant.eq.0)then
            write(6,10)0.,(string(i:i),i=1,n),char(3)
10          format(' SL',f3.1,';LB',80a1)
            return
        end if
        i1=1
15      do 20 i2=i1,n
            ss=string(i2:i2)
            if((ss.lt.'A'.and.ss.ne.' ').or.(ss.gt.'Z'.and.ss.lt.'a')
     +         .or.ss.gt.'z')go to 30
20      continue
30      i2=i2-1
        nn=i2-i1+1
        write(6,10)slant,string(i1:i2),char(3)
        i1=i2+1
        if(i1.gt.n)return
        do 40 i2=i1,n
            ss=string(i2:i2)
            if((ss.ge.'A'.and.ss.le.'Z')
     +         .or.(ss.ge.'a'.and.ss.le.'z'))go to 50
40      continue
50      i2=i2-1
        nn=i2-i1+1
        write(6,10)0.,string(i1:i2),char(3)
        i1=i2+1
        if(i1.gt.n)return
        go to 15
        end
