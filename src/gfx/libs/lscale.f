
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

      subroutine lscale(string,length,factor,ifirst)
        save
c
      character*(*) string
c
      i1 = 0
      i2 = 0
      do 100 i=1,length
	  if (string(i:i).eq.'!') then
	      if (i1.eq.0) then
		  i1 = i
	      else
		  if (i.gt.i1+1) then
		      i2 = i
		      go to 200
		  end if
	      end if
	  end if
100   continue
c
200   factor = 1.
      ifirst = 1
      if (i1.eq.0.or.i2.eq.0.or.i2.le.i1+1) return
c
      read(string(i1+1:i2-1),'(f6.3)',err=500,end=500)factor
      ifirst = i2 + 1
c
500   return
      end
