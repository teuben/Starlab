
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

	subroutine setaxes
        save
c
	common/fr setax/kax,lax
	data kax/1/lax/0/
c
c	kax:	1 ==> draw the x-axes
c		0 ==> omit the x-axes
c
c	lax:	0 ==> draw both y axes
c		1 ==> draw LH y axis only
c		2 ==> draw RH y axis only
c		< 0 ==> omit both y axes
c
	entry setxaxis(iax)
	kax = iax
	return
c
	entry setyaxis(jax)
	lax = jax
	return
c
	end
