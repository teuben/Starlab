
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
        
        subroutine setpens(i1,i2,i3)
        save
c
        logical set,pset
c
        common/fr pens/ipens(3)
        data ipens/3*0/set/.false./
c
        ipens(1)=i1
        ipens(2)=i2
        ipens(3)=i3
        set = .true.
        return
c       
        entry getpens(j1,j2,j3)
        j1=ipens(1)
        j2=ipens(2)
        j3=ipens(3)
        return
c       
        entry getpset(pset)
        pset = set
c
        end
