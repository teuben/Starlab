
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
        
        subroutine setwts(i1,i2,i3,i4)
        save
c
        logical set,wset
c
        common/fr wts/iwts(4)
        data iwts/4*0/set/.false./
c
        iwts(1)=i1
        iwts(2)=i2
        iwts(3)=i3
        iwts(4)=i4
        set = .true.
        return
c       
        entry getwts(j1,j2,j3,j4)
        j1=iwts(1)
        j2=iwts(2)
        j3=iwts(3)
        j4=iwts(4)
        return
c       
        entry getwset(wset)
        wset = set
c
        end
