
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
        
        subroutine sethts(ht1,ht2)
        save
c
        logical set,hset
c
        common/fr hts/htl,htn
        data htl/.15/htn/.15/set/.false./
c       
        if(ht1.ge.0.)htl=ht1
        if(ht2.ge.0.)htn=ht2
        set = .true.
        return
c       
        entry gethts(hh1,hh2)
        hh1=htl
        hh2=htn
        return
c       
        entry gethset(hset)
        hset = set
c
        end
