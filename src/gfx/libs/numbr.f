
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
c-----------------------------------------------------------------------
c***********************************************************************
c
c       The following routines are non-essential. The calls to
c       external routines may be commented out, if desired.
c
        subroutine numbr(rr,ss,ht,fpn,theta,ndec)
        save
c
c       Use the local number drawer to mimic "nomber."
c
        character*80 string
        common /numbr on/ inumbr
c
	r=rr
	s=ss
	go to 10
c
	entry cnumbr(ht,fpn,theta,ndec)
	call simwhe(r,s)
c
10      call numsym(fpn,ndec,string,n)
        inumbr=1
        call symbl(r,s,ht,string,theta,n)
        inumbr=0
        end
