
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

c       Routines to set and get the left-hand and bottom boundaries.
c       The boundaries are specified in standard units relative to
c       the CURRENT origin.
        
        subroutine setlhe(rl)
        save
c       
        common /plot sizes/ xsize,ysize
c
        common /lhlimit/ rlhe
        common /lowlimit/ sbot
c
        data rlhe/0./
        data sbot/0./
c       
        rlhe = rl
        return
c       
        entry getlhe(r1)
        r1 = rlhe
        return
c       
        entry getrhe(r2)
        r2 = rlhe + xsize
        return
c       
        entry setbot(sb)
        sbot = sb
        return
c       
        entry getbot(s1)
        s1 = sbot
        return
c
        entrygettop(s2)
        s2 = sbot + ysize
c       
        end
