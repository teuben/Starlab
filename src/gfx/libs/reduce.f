
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
c************************************************************************
c
c	Device-specific routines:
c	------------------------   
c
c
        subroutine reduce(i,ch)
        save
        character*3 ch
c
c       Decompose i (in [0,262143]) to base 64 (for use with idev = 1).
c
        i1=i/64
        ch(3:3)=char(32+i-64*i1)
        i2=i1/64
        ch(2:2)=char(32+i1-64*i2)
        i1=i2/64
        ch(3:3)=char(32+i2-64*i1)
        end
