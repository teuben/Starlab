
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
	program fonts
c
c	Convert from unformatted to formatted font file.
c
        integer*2 n,m,num,jl,jr,idic,long
        common /sim fc2/ n,m,num(288),jl(288),jr(288),idic(288),
     +		         long(20000)
c
	open(1,status='OLD',form='unformatted',file='SIM.UNF')
	read(1)n,m,num,jl,jr,idic,(long(i),i=1,m)
	close(1)
c
	open(1,status='NEW',form='formatted',file='SIM.FMT')
	write(1,*)n,m,num,jl,jr,idic,(long(i),i=1,m)
	close(1)
c
	end
