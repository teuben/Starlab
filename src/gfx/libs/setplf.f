
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

        subroutine setplf
        save
c
        common/fr plain/iplain/fr bare/ibare
        logical set,fset
c
        data set/.false./
c
        entry clrxtf
        entry setpln
        iplain = 1
        go to 100
c       
        entry setxtf
        entry clrplf
        entry clrpln
        iplain = 0
        ibare = 0
        go to 100
c       
        entry setbare
        ibare = 1
        iplain = 1
        go to 100
c       
        entry clrbare
        ibare = 0
        go to 100
c       
        entry setfnt(j1,j2)
        iplain = j1
        ibare = j2
        if (ibare.eq.1) iplain = 1
        if (iplain.eq.0) ibare = 0
        go to 100
c       
        entry getfnt(i1,i2)
        i1 = iplain
        i2 = ibare
c
100     set = .true.
        return
c       
        entry getfset(fset)
        fset = set
c
        end
