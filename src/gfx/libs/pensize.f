
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
        subroutine pensize(nx,ny)
        save
c
c       Establish Versaterm pen size as an (nx by ny)-pixel rectangle.
c
        character*80 device
        common /plot device/ device,aspect,idev
        character*1 ctrl(0:31),
     +              null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        common /ctrlch/ ctrl,
     +                  null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
c
        if(idev.ne.7.and.idev.ne.8.and.idev.ne.11
     +     .and.idev.ne.12)return
        if(nx.eq.1.and.ny.eq.1)then
            call output char(cr)
            return
        end if
        call pixel(nx,ny,3)
        call output char(esc)
        if(idev.eq.11.or.idev.eq.12)then
            call output char('p')                                        ! thank you, Mr. Abelbeck!
        else
            call output char('i')
        end if
        end
