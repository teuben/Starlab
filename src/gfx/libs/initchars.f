
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

        subroutine init chars
        save
c
c	General routine to set up useful character strings.
c
        character*1 ctrl(0:31),
     &              null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        common /ctrlch/ ctrl,
     &                  null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
c
        character*1 tekc(4),hpc(3),hpon(3),hpoff(6)
        common /dev controls/ tekc,hpc,hpon,hpoff
        data hpc/'P','G',';'/
     &       hpon/' ','.','y'/hpoff/'P','U',';',' ','.','Z'/
c
        do 1 i=0,31
            ctrl(i) = char(i)
1       continue
c
        null = ctrl(0)
        ctrla = ctrl(1)
	tab = ctrl(9)
	lf = ctrl(10)
        ff = ctrl(12)
        cr = ctrl(13)
        ctrlx = ctrl(24)
        ctrlz = ctrl(26)
        esc = ctrl(27)
        gs = ctrl(29)
        del = char(127)
c
        tekc(1) = gs
        tekc(2) = esc
        tekc(3) = ff
        tekc(4) = gs
        hpon(1) = esc
        hpoff(4) = esc
c
        end
