
c     
c     Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c     by Steve McMillan, Drexel University, Philadelphia, PA.
c     
c     All rights reserved.
c     
c     Redistribution and use in source and binary forms are permitted
c     provided that the above copyright notice and this paragraph are
c     duplicated in all such forms and that any documentation,
c     advertising materials, and other materials related to such
c     distribution and use acknowledge that the software was developed
c     by the author named above.
c     
c     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c     IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c     WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c     


      subroutine fr symbl(r,s,h,sym,t,n)
      save
c
      character*(*) sym
      common/fr plain/iplain
      integer*4 n
      common/debug trace/itrace
c     
      if (itrace.eq.1) write(2,*)'fr symbl:',r,s,h,sym
c
      write(6,*)'iplain = ',iplain,'  string = ',sym(1:n)
      if (iplain.eq.0) then
          call simbol(r,s,h,sym,t,n)
      else
          call symbl(r,s,h,sym,t,n)
      end if
c
      end
